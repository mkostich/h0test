## Flatten arguments destined for cat() into character scalars, one per argument.
##   cat() cannot handle language objects (e.g. formulas, calls) or non-atomic
##   objects (e.g. lists), which otherwise mask the message being logged with a
##   confusing cat() error; those are deparsed instead:

f.cat_args <- function(...) {
  args <- list(...)
  out <- lapply(args, function(arg) {
    if(is.null(arg) || length(arg) %in% 0) return(character(0))
    if(is.numeric(arg)) return(paste(format(arg, trim=T), collapse=" "))
    if(is.character(arg) || is.logical(arg) || is.factor(arg)) {
      return(paste(arg, collapse=" "))
    }
    paste(deparse(arg), collapse=" ")
  })
  return(unlist(out))
}

f.msg <- function(..., config) {
  if(is.null(config$log_file)) config$log_file <- ""
  cat(f.cat_args(...), "\n", file=config$log_file, append=T)
  utils::flush.console()
}

f.log <- function(..., config) {
  if(is.null(config$log_file)) config$log_file <- ""
  cat(f.cat_args(...), "at:", format(Sys.time(), format='%Y%m%d%H%M%S'), "\n",
    file=config$log_file, append=T)
  utils::flush.console()
}

f.log_block <- function(..., config) {
  f.msg("", config=config)              ## blank line
  f.log(..., config=config)
}

f.err <- function(..., config) {
  f.msg("", config=config)              ## blank line
  f.log("ERROR:", ..., config=config)
  stop("Stopping", call.=F)
}

f.log_obj <- function(obj, config) {
  
  if(is.null(config$log_file)) config$log_file <- ""
  
  if(config$log_file %in% "") {
    print(obj)
    utils::flush.console()
  } else {
    utils::capture.output(obj, file=config$log_file, append=T)
  }
}

f.save_tsv <- function(dat, file_out, config, row.names=T, col.names=T) {
  tryCatch(
    utils::write.table(dat, file=file_out, quote=F, sep="\t", 
      row.names=row.names, col.names=col.names),
      error=function(msg) f.err("write.table() error: writing to ", 
        file_out, ": ", msg$message, config=config),
      warning=function(msg) f.err("write.table() warning: writing to ", 
        file_out, ": ", msg$message, config=config)
  )
}

## TRUE iff x is a single formula, e.g. ~age+gender+age:gender or y~age.
##   NOTE: length() of a formula is the number of components of the underlying
##   call: 2 for a one-sided formula (`~` and rhs), 3 for a two-sided formula
##   (`~`, lhs, and rhs); it is never 1. Two-sided formulas are accepted, but
##   the dependent variable is dropped by f.parse_frm(), since the dependent is
##   always the expression values of one feature:

f.is_formula <- function(x) {
  return(inherits(x, "formula") && length(x) %in% c(2, 3))
}

## Operators supported in config$frm. Everything else ('|', '^', '/', '%in%',
##   parentheses, and any function call, like log() or I()) is rejected by
##   f.check_frm_expr(). '-' is supported only for removing the intercept:

f.frm_ops <- c("+", ":", "*", "-")

## Recursively check one expression from the right-hand side of a formula.
##   Throws an informative error on any unsupported construct; otherwise
##   returns TRUE invisibly:

f.check_frm_expr <- function(expr, config, top=NULL) {

  if(is.null(top)) top <- expr

  ## terminal symbol: a bare variable name; '.' would silently pull in every
  ##   column of the sample metadata, so is not supported:

  if(is.symbol(expr)) {
    if(as.character(expr) %in% ".") {
      f.err("f.check_frm_expr: '.' not supported in formula; frm rhs:",
        deparse(top), config=config)
    }
    return(invisible(TRUE))
  }

  ## terminal constant: only 0 and 1 are meaningful (intercept):

  if(is.numeric(expr) && length(expr) %in% 1) {
    if(!(expr %in% c(0, 1))) {
      f.err("f.check_frm_expr: constant", expr, "not supported in formula;",
        "only 0 and 1 (intercept) allowed; frm rhs:", deparse(top),
        config=config)
    }
    return(invisible(TRUE))
  }

  if(!is.call(expr)) {
    f.err("f.check_frm_expr: unsupported element", deparse(expr),
      "in formula; frm rhs:", deparse(top), config=config)
  }

  op <- as.character(expr[[1]])

  if(!(op %in% f.frm_ops)) {
    f.err("f.check_frm_expr: operator or function '", op,
      "' not supported in formula; supported operators:",
      paste(f.frm_ops, collapse=" "),
      "; variables must be bare column names (no transformations, like",
      "log() or I()); frm rhs:", deparse(top), config=config)
  }

  ## '-' only for dropping the intercept, so what is removed must be 1 or 0;
  ##   covers both binary (x-1) and unary (-1+x) minus:

  if(op %in% "-") {
    rhs <- expr[[length(expr)]]
    if(!(is.numeric(rhs) && length(rhs) %in% 1 && rhs %in% c(0, 1))) {
      f.err("f.check_frm_expr: '-' supported only for removing the intercept",
        "(e.g. ~x-1); cannot remove term", deparse(rhs), "; frm rhs:",
        deparse(top), config=config)
    }
  }

  for(idx in 2:length(expr)) f.check_frm_expr(expr[[idx]], config, top=top)

  return(invisible(TRUE))
}

## Sort the variables within each interaction term label, so that e.g.
##   'sex:age' and 'age:sex' both become 'age:sex'. Facilitates comparison of
##   terms and of a term with config$test_term:

f.canon_label <- function(labels) {
  if(length(labels) %in% 0) return(character(0))
  out <- vapply(
    labels,
    function(lbl) paste(sort(unlist(strsplit(lbl, ":", fixed=T))), collapse=":"),
    character(1)
  )
  return(unname(out))
}

## Parse and check a model formula. Single place where config$frm is
##   interpreted. Throws an informative error on unsupported constructs.
##   Returns a list:
##     $frm       one-sided formula, with any dependent variable dropped
##     $vars      character vector of variables (bare column names) in $frm
##     $labels    character vector of canonical term labels (see f.canon_label)
##     $factors   variables x terms incidence matrix from stats::terms()
##     $intercept 1 if model includes an intercept, else 0
##     $two_sided TRUE iff a dependent variable was given (and dropped)

f.parse_frm <- function(frm, config) {

  if(!f.is_formula(frm)) {
    f.err("f.parse_frm: config$frm not a formula; value:", frm, config=config)
  }

  ## dependent variable is always the expression values of one feature, so a
  ##   dependent given in frm carries no information; drop it, so that it is
  ##   not mistaken for a covariate downstream:

  two_sided <- length(frm) %in% 3
  if(two_sided) frm <- frm[-2]

  ## '|' would otherwise be reported as an unsupported '(', which is unhelpful,
  ##   since it is nearly always an attempt at a random effect:

  if(any(grepl("|", deparse(frm), fixed=T))) {
    f.err("f.parse_frm: '|' not supported in formula; mixed models (random",
      "effects) are not implemented; frm:", frm, config=config)
  }

  f.check_frm_expr(frm[[2]], config)

  trms <- stats::terms(frm)
  labels <- f.canon_label(attr(trms, "term.labels"))
  fmat <- attr(trms, "factors")

  if(length(labels) %in% 0) {
    fmat <- matrix(0, nrow=0, ncol=0)
  } else {
    fmat <- matrix(fmat, nrow=nrow(fmat), ncol=ncol(fmat),
      dimnames=list(rownames(fmat), labels))
  }

  if(any(duplicated(labels))) {
    f.err("f.parse_frm: duplicated term '", labels[duplicated(labels)][1],
      "' in frm:", frm, config=config)
  }

  out <- list(
    frm=frm,
    vars=all.vars(frm),
    labels=labels,
    factors=fmat,
    intercept=attr(trms, "intercept"),
    two_sided=two_sided
  )

  return(out)
}

## Rank of a design matrix, tolerant of the degenerate shapes that arise when a
##   feature was measured in too few observations, or when the reduced model has
##   no columns left:

f.design_rank <- function(mat) {
  if(nrow(mat) %in% 0 || ncol(mat) %in% 0) return(0L)
  return(qr(mat)$rank)
}

## The design matrix for config$frm, the columns of it that carry the test of
##   config$test_term, and the degrees of freedom of that test measured over all
##   observations. Single place where the full and reduced models are derived
##   from config, so that filter_features_by_estimability() and the hypothesis
##   tests cannot disagree about what is being tested.
##   A df_intend of zero means the reduced model spans the same column space as
##   the full model, so the two fits are indistinguishable and there is no
##   hypothesis to test; any p-value reported would be an artifact of how the
##   models were coded rather than a statement about the data. That is an error
##   rather than a warning, since there is no result to salvage. It is reached
##   whenever config$test_term contributes no column that the remaining terms do
##   not already imply, most easily by naming the intercept alongside a factor
##   that the reduced model codes to full rank: with frm ~grp and test_term '1',
##   'y ~ 1 + grp' and 'y ~ 0 + grp' are two codings of one model:

f.design_test_cols <- function(state, config) {

  ## throws an informative error if config$test_term does not fit config$frm:

  parsed <- f.parse_frm(config$frm, config)
  drops <- f.normalize_terms(config)$drop_terms

  X <- stats::model.matrix(parsed$frm, data=state$samples)

  if(nrow(X) != ncol(state$expression)) {
    f.err("f.design_test_cols: design matrix has", nrow(X),
      "rows, but state$expression has", ncol(state$expression), "columns;",
      "\n", "model.matrix() drops observations with missing covariate values",
      config=config)
  }

  ## attr(X, 'assign') indexes the term labels in order, with 0 for the
  ##   intercept; parsed$labels is those same labels, canonicalized:

  asgn <- attr(X, "assign")
  cols_test <- which(asgn %in% match(setdiff(drops, "1"), parsed$labels))
  if("1" %in% drops) cols_test <- sort(c(which(asgn %in% 0), cols_test))

  if(length(cols_test) %in% 0) {
    f.err("f.design_test_cols: config$test_term", config$test_term,
      "matches no columns of the design matrix;", "\n",
      "terms dropped to form the reduced model:", drops, "\n",
      "terms of config$frm:", parsed$labels, config=config)
  }

  rank_all <- f.design_rank(X)
  rank_red <- f.design_rank(X[, -cols_test, drop=F])
  df_intend <- rank_all - rank_red

  ## nothing left in the reduced model, which happens when config$frm suppresses
  ##   the intercept and config$test_term names every remaining term. The test is
  ##   then against zero rather than against a common mean: well defined, but on
  ##   log-scale abundances every feature rejects it, which is rarely the
  ##   question. Warn rather than stop, since the test asked for is the one
  ##   performed, and it is also what filter_features_by_estimability() screens:

  if(rank_red %in% 0) {
    f.msg("WARNING: f.design_test_cols: dropping config$test_term '",
      config$test_term, "' leaves a reduced model with no parameters, so the",
      "test is of whether the", config$test_term, "means are all zero, not of",
      "whether they differ from each other;", "\n",
      "config$frm:", deparse(parsed$frm), "; for the usual comparison among",
      "levels, keep the intercept in config$frm", config=config)
  }

  if(df_intend %in% 0) {
    f.err("f.design_test_cols: dropping config$test_term '", config$test_term,
      "' leaves a reduced model spanning the same space as the full model,",
      "so there is no hypothesis to test;", "\n",
      "config$frm:", deparse(parsed$frm), "; terms dropped to form the",
      "reduced model:", paste(drops, collapse=" "), "; design columns:",
      ncol(X), "; rank:", rank_all, "\n",
      "testing the intercept ('1') is only meaningful when no factor in",
      "config$frm is coded to full rank in the reduced model", config=config)
  }

  out <- list(
    parsed=parsed,
    drops=drops,
    X=X,
    cols_test=cols_test,
    rank_all=rank_all,
    df_intend=df_intend
  )

  return(out)
}

## How many design matrix columns a test method can test at once. The limma-family
##   methods take a vector of coefficients and return an F-test over all of them, so
##   they are unlimited. proDA is unlimited too, by a different route:
##   proDA::test_diff() takes either one contrast or a reduced model, and test_proda()
##   hands it the reduced model that f.design_test_cols() built whenever more than one
##   column carries the test, which is a likelihood ratio test over all of them. Two
##   engines can still only report one coefficient at a time. DEqMS moderates a
##   single coefficient's t-statistic: DEqMS::spectraCounteBayes() forms sca.t from
##   fit$coefficients[, coef_col] over fit$stdev.unscaled[, coef_col], and the
##   package has no F-analogue anywhere, so the limit is in the moderation itself.
##   DEqMS::outputResult() takes a single coef_col because there is nothing joint
##   for it to report, so it is the wrong place to look for the constraint.
##   msqrob2::hypothesisTest() loops over the columns of the contrast and returns
##   one table per column rather than a joint test over several; that one is an API
##   limit rather than a statistical one, since msqrob2 exports getCoef(),
##   getVcovUnscaled(), getVarPosterior() and getDfPosterior(), which is everything
##   a joint Wald test would need. Both are limitations of the engines, not of how
##   h0testr selects columns:

f.test_max_cols <- function(method) {
  if(method %in% c("deqms", "msqrob")) return(1L)
  return(Inf)
}

## How many terms of config$frm a test method can test at once; companion to
##   f.test_max_cols(), which counts design matrix columns. The two limits are
##   different things and a method can be bounded by either: prolfqua reports the
##   rows of an anova table, one per term, so a factor with several levels is
##   already a correct multi-df F-test there, but a joint test over several terms
##   is not something the table can express. So testing 'dose' in ~dose with a
##   three level dose runs under prolfqua (one term, 2 df) while deqms and msqrob
##   cannot (2 coefficients), and testing 'sex' in ~sex*batch is out of reach for
##   all three; proda is bounded by neither limit, testing several columns by
##   likelihood ratio against the reduced model:

f.test_max_terms <- function(method) {
  if(method %in% "prolfqua") return(1L)
  return(Inf)
}

## The design for config$frm and the columns of it carrying the test of
##   config$test_term, for a method that can test at most max_cols of them at once.
##   Wraps f.design_test_cols() so that a method which cannot express the test that
##   config$test_term implies says so, rather than quietly testing a narrower
##   hypothesis than the one that was asked for and that
##   filter_features_by_estimability() screened features against. The shortfall
##   arises from the marginality rule: testing a variable tests every term
##   containing it, so testing 'sex' in ~sex*batch is a joint test of sexM and
##   sexM:batchb2, and testing a factor with more than two levels is a joint test of
##   its contrasts:

f.design_test_cols_max <- function(state, config, caller, max_cols=Inf) {

  design <- f.design_test_cols(state, config)
  n_cols <- length(design$cols_test)

  if(n_cols > max_cols) {
    f.err(caller, ": testing config$test_term '", config$test_term, "' in",
      deparse(design$parsed$frm), "is a joint test of", n_cols, "coefficients (",
      paste(colnames(design$X)[design$cols_test], collapse=", "), "), but", caller,
      "can test at most", max_cols, "at a time;", "\n",
      "use test_method 'lm', 'trend', 'voom' or 'proda' for this test_term, or",
      "name a term of config$frm that resolves to a single coefficient",
      config=config)
  }

  return(design)
}

## Classify each variable in config$frm as "factor" or "numeric" (continuous).
##   Single place where covariate type is decided, so that value checking,
##   filtering, and hypothesis testing all agree. Rules, in order:
##     1. named in config$reference_levels: "factor"; a numeric column may be
##          declared there, which is how a numeric variable is made categorical;
##     2. already a factor (e.g. set by f.set_covariate_factor_levels()): "factor";
##     3. logical: "factor"; model.matrix() orders these FALSE, TRUE, which is
##          deterministic, so no declaration is needed;
##     4. numeric: "numeric" (continuous);
##     5. anything else (typically character): error. Deriving the reference
##          level by sorting the values is locale dependent, and silently sets
##          the meaning of the reported coefficients, so the reference level has
##          to be declared in config$reference_levels.
##   Returns a named character vector with one element per variable in
##   config$frm. config$covariate_types is reused when it covers exactly the
##   variables in config$frm; since the caller controls the calling order, a
##   cached value that does not match is ignored and recomputed:

f.covariate_types <- function(state, config) {

  vars <- sort(unique(f.parse_frm(config$frm, config)$vars))

  if(!all(vars %in% names(state$samples))) {
    f.err("f.covariate_types: !all(vars %in% names(state$samples)); vars:",
      vars, "; names(state$samples):", names(state$samples), config=config)
  }

  cached <- config$covariate_types
  if(!is.null(cached) && is.character(cached) && !is.null(names(cached)) &&
      setequal(names(cached), vars) && all(cached %in% c("factor", "numeric"))) {
    return(cached[vars])
  }

  out <- character(0)

  for(nom in vars) {

    v <- state$samples[[nom]]

    if(nom %in% names(config$reference_levels)) {
      out[nom] <- "factor"
    } else if(is.factor(v) || is.logical(v)) {
      out[nom] <- "factor"
    } else if(is.numeric(v)) {
      out[nom] <- "numeric"
    } else {
      f.err("f.covariate_types: covariate", nom, "in config$frm is of class",
        class(v), "and is not declared in config$reference_levels;", "\n",
        "add", nom, "to config$reference_levels to set its reference level;",
        "\n", "distinct values:",
        utils::head(sort(unique(as.character(v))), 10), config=config)
    }
  }

  return(out)
}

## Rebuild a factor covariate with the level ordering initialize() resolved into
##   config$factor_levels, with the declared reference level first. Handing a
##   downstream fit a character vector instead leaves it to re-derive the levels by
##   sorting, which renames the coefficients: with reference level "M", h0testr names
##   the tested column sexF while an alphabetical re-derivation names it sexM. Where
##   the fit reports one named contrast (test_msqrob()) that mismatch yields a table
##   of NAs with no error; where it reports a coefficient (test_prolfqua()) it changes
##   which contrast the reported numbers describe. Both are silent, hence this.
##   config$factor_levels is absent when initialize() has not run, in which case the
##   sorted levels are the best available and match what the fit would derive anyway:

f.relevel_covariate <- function(v, trm, config, caller) {

  lvls <- config$factor_levels[[trm]]
  if(is.null(lvls)) lvls <- levels(factor(v))

  out <- factor(as.character(v), levels=lvls)

  if(any(is.na(out))) {
    f.err(caller, ": covariate", trm, "has values that are not among the levels",
      "set by initialize();", "\n", "levels:", lvls, "\n", "unmatched values:",
      utils::head(sort(unique(as.character(v)[is.na(out)])), 10), config=config)
  }

  return(out)
}

## Check the values of the covariates in config$frm. Missing, non-finite, and
##   constant covariates are errors: they cannot be fit, and letting them
##   through means model.matrix() silently drops observations (so that the
##   number of observations differs between features) or yields a
##   rank-deficient design. A numeric covariate with few distinct values is
##   often a miscoded factor, so is warned about; threshold is
##   config$n_distinct_numeric_warn. Optional types from f.covariate_types():

f.check_covariate_values <- function(state, config, types=NULL) {

  if(is.null(types)) types <- f.covariate_types(state, config)

  ## observation labels, for reporting which observations are offending:

  obs <- NULL
  for(nom in c(config$obs_col, config$obs_id_col)) {
    if(length(nom) %in% 1 && nchar(nom) > 0 && nom %in% names(state$samples)) {
      obs <- as.character(state$samples[[nom]])
      break
    }
  }
  if(is.null(obs)) obs <- colnames(state$expression)
  if(length(obs) != nrow(state$samples)) obs <- as.character(1:nrow(state$samples))

  cutoff <- config$n_distinct_numeric_warn
  if(is.null(cutoff)) cutoff <- 5

  for(nom in names(types)) {

    v <- state$samples[[nom]]

    i <- is.na(v)
    if(any(i)) {
      f.err("f.check_covariate_values: covariate", nom, "has missing values,",
        "which are not supported;", "\n", "n missing:", sum(i),
        "; first offending observations:", utils::head(obs[i], 10),
        config=config)
    }

    if(types[nom] %in% "numeric" && !all(is.finite(v))) {
      i <- !is.finite(v)
      f.err("f.check_covariate_values: covariate", nom, "has non-finite",
        "values;", "\n", "n non-finite:", sum(i),
        "; first offending observations:", utils::head(obs[i], 10),
        config=config)
    }

    lvls <- unique(v)

    if(length(lvls) %in% 1) {
      f.err("f.check_covariate_values: covariate", nom, "is constant, with",
        "single distinct value:", utils::head(as.character(lvls), 1), ";", "\n",
        "it cannot be fit; drop it from config$frm", config=config)
    }

    if(types[nom] %in% "numeric" && length(lvls) <= cutoff) {
      f.msg("WARNING: numeric (continuous) covariate", nom, "has only",
        length(lvls), "distinct values:", sort(as.character(lvls)), "\n",
        "  if it is categorical, declare it in config$reference_levels;",
        "otherwise this warning can be ignored", config=config)
    }
  }

  return(invisible(TRUE))
}

## needs config$feat_col and config$obs_col:

f.check_state <- function(state, config) {

  if(is.null(config$feat_col) || config$feat_col %in% "") {
    f.err("f.check_state: config$feat_col unset", config=config) 
  }
  feats <- state$features[[config$feat_col]]
  if(!all(rownames(state$expression) == feats)) {
    f.err("f.check_state: state$features do not match rows of state$expression", 
      "config$feat_col: ", config$feat_col, config=config)
  }
  
  if(is.null(config$obs_col) || config$obs_col %in% "") {
    f.err("f.check_state: config$obs_col unset", config=config)
  }
  samps <- state$samples[[config$obs_col]]
  if(!all(colnames(state$expression) == samps)) {
    f.err("f.check_state: state$samples do not match columns of state$expression", 
      "config$obs_col: ", config$obs_col, config=config)
  }
  
  if(!is.matrix(state$expression)) {
    f.err("f.check_state: !is.matrix(state$expression)", config=config)
  }

  ## NA is the only indicator of a missing value; see f.zeros_to_na(). Where
  ##   normalize() was able to guarantee that an exact 0 cannot be a measurement,
  ##   it says so in config$log_from_raw, and one appearing afterwards can only
  ##   have been written by code using 0 to mean missing. That is easy to do by
  ##   accident: combining a replicate group with sum(na.rm=TRUE) returns 0 for a
  ##   group in which nothing was measured, and on a log scale 0 is not a neutral
  ##   value but the most extreme one in the matrix. Checked at every step
  ##   boundary because the value is indistinguishable from a real measurement
  ##   once anything downstream has read it. Unset means no guarantee, so no
  ##   check, which is also what a minimal config gets:

  if(isTRUE(config$log_from_raw)) {

    i_zero <- which(state$expression %in% 0)   ## NAs are not matched

    if(length(i_zero) > 0) {

      rnom <- rownames(state$expression)
      cnom <- colnames(state$expression)
      if(is.null(rnom)) rnom <- as.character(1:nrow(state$expression))
      if(is.null(cnom)) cnom <- as.character(1:ncol(state$expression))

      idxs <- utils::head(i_zero, 5)
      rr <- ((idxs - 1) %% nrow(state$expression)) + 1
      cc <- ((idxs - 1) %/% nrow(state$expression)) + 1

      f.err("f.check_state: exact zeros in state$expression;",
        "the data are log2(x + 1) of raw input whose zeros became NA, so a 0",
        "here is a missing value written as a measurement, not a measurement;",
        "\n", "zeros:", length(i_zero), "of", length(state$expression), ";",
        "first offenders (feature, observation):", "\n",
        paste(rnom[rr], cnom[cc], sep=", "), config=config)
    }
  }
}

## Resolve the scale of state$expression for a function that takes an
##   is_log_transformed argument. The argument wins when given, so that the
##   function can be called on its own with a config that says nothing about the
##   scale; otherwise config$is_log_transformed answers, having been set by
##   initialize() and updated by normalize(). Disagreement between the two is an
##   error rather than a silent preference, since it means the caller and the
##   workflow hold different beliefs about the data and only one of them can be
##   right:

f.is_log_transformed <- function(is_log_transformed, config, fn_name) {

  from_config <- config$is_log_transformed

  if(!is.null(from_config) &&
      !(is.logical(from_config) && length(from_config) %in% 1 &&
        !is.na(from_config))) {
    f.err(fn_name, ": config$is_log_transformed is not TRUE or FALSE; value:",
      from_config, config=config)
  }

  ## "" is accepted as unset for callers that pass an empty character:

  unset <- is.null(is_log_transformed) ||
    (is.character(is_log_transformed) && all(is_log_transformed %in% ""))

  if(unset) {
    if(is.null(from_config)) {
      f.err(fn_name, ": is_log_transformed and config$is_log_transformed both",
        "unset;", "\n",
        "set config$is_log_transformed to declare the scale of the data, or",
        "pass is_log_transformed to this function", config=config)
    }
    return(from_config)
  }

  if(!(is.logical(is_log_transformed) && length(is_log_transformed) %in% 1 &&
      !is.na(is_log_transformed))) {
    f.err(fn_name, ": is_log_transformed is not TRUE or FALSE; value:",
      is_log_transformed, "; typeof:", typeof(is_log_transformed),
      config=config)
  }

  if(!is.null(from_config) && !identical(is_log_transformed, from_config)) {
    f.err(fn_name, ": is_log_transformed argument is", is_log_transformed,
      "but config$is_log_transformed is", from_config, ";", "\n",
      "these describe the same data and cannot both be right; drop the",
      "argument to use the scale the workflow recorded, or correct",
      "config$is_log_transformed", config=config)
  }

  return(is_log_transformed)
}

## Drop the features that cannot contribute to a fit of missingness against
##   intensity, and refuse the fit outright if too little is left of it. A
##   feature measured nowhere has no observed intensity, so the f_mid summary of
##   it is NA. Substituting a stand-in value would enter it into the fit as a
##   genuine point, sitting at the extreme of the response axis where it carries
##   the most leverage over the slope, so it is dropped instead; it is still
##   imputed afterwards, from the curve the remaining features determine.
##   h0testr::filter() rejects features measured nowhere, so this normally only
##   arises when an imputer is called on its own. Degeneracy is judged on the
##   number of distinct intensities rather than the number of features, since
##   many features sharing one intensity still cannot identify a slope. The fit
##   is extrapolated across the whole intensity range and imputed values are
##   drawn from it, so a curve resting on a handful of points does not merely
##   estimate badly, it invents structure that then enters the data as
##   measurements:

f.drop_unfittable <- function(dat, config, fn_name, min_fit_pts) {

  if(!(is.numeric(min_fit_pts) && length(min_fit_pts) %in% 1 &&
      is.finite(min_fit_pts) && min_fit_pts >= 2)) {
    f.err(fn_name, ": min_fit_pts is not a finite numeric scalar >= 2; value:",
      min_fit_pts, config=config)
  }

  i_drop <- is.na(dat$m)

  if(any(i_drop)) {

    nom <- rownames(dat)[i_drop]
    if(is.null(nom)) nom <- as.character(which(i_drop))

    f.msg(fn_name, ": dropping", sum(i_drop), "of", nrow(dat), "features from",
      "the missingness fit; nothing was measured for them, so they have no",
      "intensity to fit against;", "\n",
      "  they are still imputed, from the fit the remaining features give;",
      "first:", paste(utils::head(nom, 5), collapse=", "), config=config)

    dat <- dat[!i_drop, , drop=F]
  }

  n_distinct <- length(unique(dat$m))

  if(n_distinct < min_fit_pts) {
    f.err(fn_name, ": too few distinct intensities to fit missingness against;",
      "distinct intensities:", n_distinct, "; required:", min_fit_pts, ";",
      "features retained:", nrow(dat), "\n",
      "  the fit is extrapolated over the whole intensity range and imputed",
      "values are drawn from it, so this few points would invent structure",
      "instead of estimating it;", "\n",
      "  use a different config$impute_method, or lower min_fit_pts if you",
      "understand the consequences", config=config)
  }

  return(dat)
}

## The lower bound of the interval the unif_ imputation methods draw from. A
##   missing value means the feature fell below detection in that sample, so the
##   bound belongs at the bottom of the scale the data are on. On the raw scale
##   that is zero abundance, a fixed point of the measurement that needs no
##   configuration. A log scale has no such point: zero there is a single count,
##   which after normalization usually sits near the top of the range rather than
##   the bottom, so the bound is placed relative to the dimmest value actually
##   measured and displaced by config$impute_floor_offset. That offset is also
##   what gives the interval any width when config$impute_quantile is 0, since
##   the upper bound is then the observed minimum itself:

f.impute_floor <- function(mat, config, fn_name, is_log_transformed=NULL) {

  is_log_transformed <- f.is_log_transformed(is_log_transformed, config, fn_name)

  if(!is_log_transformed) return(0)

  offset <- config$impute_floor_offset
  if(is.null(offset)) offset <- -1              ## the new_config() default

  if(!(is.numeric(offset) && length(offset) %in% 1 && is.finite(offset) &&
      offset <= 0)) {
    f.err(fn_name, ": config$impute_floor_offset is not a finite non-positive",
      "numeric scalar; value:", offset, config=config)
  }

  min_val <- suppressWarnings(min(mat, na.rm=T))

  if(!is.finite(min_val)) {
    f.err(fn_name, ": no finite observed value to place the imputation floor",
      "below, so there is nothing to impute from; min:", min_val, config=config)
  }

  return(min_val + offset)
}

f.report_state <- function(state, config) {

  f.msg("class(state$expression):", class(state$expression), config=config)

  f.msg("N features: ", nrow(state$expression), 
    "; N observations: ", ncol(state$expression), config=config)
    
  f.msg("signal distribution:", config=config)
  
  f.quantile(c(state$expression), config, digits=0) 
  
  f.msg("min(state$expression):", min(c(state$expression), na.rm=T), 
    "; mean(state$expression):", mean(c(state$expression), na.rm=T), 
    config=config)
    
  f.msg("num NAs: ", sum(is.na(c(state$expression))), config=config)
  f.msg("num non-NAs: ", sum(!is.na(c(state$expression))), config=config)
}

f.save_state <- function(state, config, prefix) {

  if(is.null(config$save_state) || !config$save_state) {
    f.msg("config$save_state not TRUE; no files saved", config=config)
    return(NULL)
  }
  
  file_out <- paste0(config$dir_out, "/", prefix, 
    config$data_mid_out, config$suffix_out)
  f.log("writing expression data to", file_out, config=config)
  f.save_tsv(state$expression, file_out, config)

  file_out <- paste0(config$dir_out, "/", prefix, 
    config$feature_mid_out, config$suffix_out)
  f.log("writing feature metadata to", file_out, config=config)
  f.save_tsv(state$features, file_out, config)

  file_out <- paste0(config$dir_out, "/", prefix, 
    config$sample_mid_out, config$suffix_out)
  f.log("writing sample metadata to", file_out, config=config)
  f.save_tsv(state$samples, file_out, config)
}

f.quantile <- function(v, config, probs=NULL, digits=3, na.rm=T) {

  if(is.null(probs)) probs <- config$probs
  if(is.null(probs)) probs <- c(0, 0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99, 1.0)
  if(is.null(config$log_file)) config$log_file <- ""
  
  if(config$log_file %in% "") {
    print(round(stats::quantile(v, probs=probs, na.rm=na.rm), digits=digits))
    utils::flush.console()
  } else {
    utils::capture.output(round(stats::quantile(v, probs=probs, na.rm=na.rm), digits=digits), 
      file=config$log_file, append=T)
  }
}

