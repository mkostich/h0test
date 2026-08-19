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

## TRUE iff config$contrast asks for a contrast. "" (the default) means none, the
##   same convention config$permute_var uses. Single place where that is decided,
##   so that the config checks, the estimability filter and every hypothesis test
##   agree about which of the two hypotheses a run is testing:

f.contrast_set <- function(config) {
  if(length(config$contrast) != 1) return(FALSE)
  if(is.na(config$contrast)) return(FALSE)
  return(nzchar(trimws(config$contrast)))
}

## What a run is testing, for progress messages: config$test_term, or the contrast
##   when config$contrast is set, the two being mutually exclusive:

f.test_label <- function(design, config) {
  if(!is.null(design$contrast)) return(paste("contrast:", trimws(config$contrast)))
  return(paste("test_term:", config$test_term))
}

## Operators supported in config$contrast. A contrast is arithmetic over
##   coefficient names rather than a model formula, so the supported set differs
##   from f.frm_ops: '(' and '/' are allowed here, since averaging a group of
##   levels ('(grpb + grpc)/2') is the usual way to write one side of a contrast,
##   while ':' and '*' are not interaction operators here. An interaction
##   coefficient is named by one backquoted name ('`sexM:batchb2`') rather than
##   built from its variables, since what a contrast weights is a column of the
##   design matrix that config$frm has already produced:

f.contrast_ops <- c("+", "-", "*", "/", "(")

## TRUE iff expr mentions any coefficient name; helper for f.check_contrast_expr(),
##   which needs to know whether an operand of '*' or '/' is a constant:

f.contrast_has_name <- function(expr) {
  if(is.symbol(expr)) return(TRUE)
  if(!is.call(expr) || length(expr) < 2) return(FALSE)
  for(idx in 2:length(expr)) {
    if(f.contrast_has_name(expr[[idx]])) return(TRUE)
  }
  return(FALSE)
}

## Recursively check one expression from config$contrast against cols, the
##   coefficient names of the design matrix for config$frm. Throws an informative
##   error on any unsupported construct; otherwise returns TRUE invisibly.
##   Same shape as f.check_frm_expr(), and for the same reason: the expression is
##   evaluated below, so what it may contain is decided here rather than by
##   whatever base::eval() would accept:

f.check_contrast_expr <- function(expr, cols, config, top=NULL) {

  if(is.null(top)) top <- expr

  ## terminal symbol: a coefficient name, which must be a column of the design.
  ##   Checking here rather than after evaluation is what makes a misspelled or
  ##   mis-cased level name say so, and say what was available instead:

  if(is.symbol(expr)) {

    nom <- as.character(expr)

    if(!(nom %in% cols)) {
      f.err("f.check_contrast_expr: '", nom, "' in config$contrast is not a",
        "coefficient of the design matrix for config$frm;", "\n",
        "config$contrast:", deparse(top), "\n",
        "coefficients available:", paste(cols, collapse=", "), "\n",
        "a coefficient whose name is not a syntactic name, such as an",
        "interaction, has to be backquoted, e.g. `sexM:batchb2`", config=config)
    }

    return(invisible(TRUE))
  }

  ## terminal constant: any finite number, unlike in a formula, where only 0 and 1
  ##   (the intercept) mean anything. Here a constant is a weight:

  if(is.numeric(expr) && length(expr) %in% 1) {
    if(!is.finite(expr)) {
      f.err("f.check_contrast_expr: constant", expr, "in config$contrast is not",
        "finite; config$contrast:", deparse(top), config=config)
    }
    return(invisible(TRUE))
  }

  if(!is.call(expr)) {
    f.err("f.check_contrast_expr: unsupported element", deparse(expr),
      "in config$contrast; config$contrast:", deparse(top), config=config)
  }

  op <- as.character(expr[[1]])

  if(!(op %in% f.contrast_ops)) {
    f.err("f.check_contrast_expr: operator or function '", op,
      "' not supported in config$contrast; supported operators:",
      paste(f.contrast_ops, collapse=" "),
      "; a contrast is a weighted sum of coefficients, so coefficient names must",
      "be bare (backquoted if not syntactic) and may be scaled by constants",
      "only; config$contrast:", deparse(top), config=config)
  }

  ## '*' and '/' scale a contrast by a constant. Between two coefficient names
  ##   they would multiply the two weight vectors elementwise, which is not a
  ##   linear combination of the coefficients at all, and would silently yield a
  ##   contrast of all zeros for two different coefficients:

  if(op %in% c("*", "/")) {

    if(length(expr) != 3) {
      f.err("f.check_contrast_expr: '", op, "' needs two operands in",
        "config$contrast; got", deparse(expr), "; config$contrast:",
        deparse(top), config=config)
    }

    lhs <- f.contrast_has_name(expr[[2]])
    rhs <- f.contrast_has_name(expr[[3]])

    if(op %in% "/" && rhs) {
      f.err("f.check_contrast_expr: the divisor of '/' in config$contrast must be",
        "a constant, not an expression in coefficient names; got",
        deparse(expr[[3]]), "; config$contrast:", deparse(top), config=config)
    }

    if(op %in% "*" && lhs && rhs) {
      f.err("f.check_contrast_expr: '*' in config$contrast multiplies two",
        "expressions in coefficient names, which is not a weighted sum of",
        "coefficients; one side must be a constant; got", deparse(expr),
        "; config$contrast:", deparse(top), config=config)
    }
  }

  for(idx in 2:length(expr)) {
    f.check_contrast_expr(expr[[idx]], cols, config, top=top)
  }

  return(invisible(TRUE))
}

## The weight vector L of config$contrast over cols, the coefficient names of the
##   design matrix for config$frm. The hypothesis tested is that the weighted sum
##   of the coefficients is zero, which is one degree of freedom whatever
##   config$frm looks like. Evaluated by binding each coefficient name to its own
##   indicator vector over cols and letting '+', '-', '*', '/' and '(' do ordinary
##   vector arithmetic, so that the weights are built by the same rules a reader
##   of the expression would apply. The environment's parent is baseenv() so that
##   those operators resolve and nothing else does:

f.contrast_vector <- function(config, cols, caller="f.contrast_vector") {

  txt <- trimws(config$contrast)
  expr <- try(str2lang(txt), silent=T)

  if(inherits(expr, "try-error")) {
    f.err(caller, ": config$contrast does not parse as an R expression;",
      "config$contrast:", txt, "\n", "parse error:",
      conditionMessage(attr(expr, "condition")), "\n",
      "expected a weighted sum of coefficient names, e.g. 'grpb - grpc' or",
      "'(grpb + grpc)/2 - grpd'", config=config)
  }

  f.check_contrast_expr(expr, cols, config)

  env <- new.env(parent=baseenv())

  for(idx in seq_along(cols)) {
    wts <- rep(0, length(cols))
    wts[idx] <- 1
    assign(cols[idx], wts, envir=env)
  }

  out <- try(eval(expr, envir=env), silent=T)

  if(inherits(out, "try-error") || !is.numeric(out) ||
      length(out) != length(cols) || any(!is.finite(out))) {
    f.err(caller, ": config$contrast did not evaluate to one weight per",
      "coefficient of the design matrix;", "config$contrast:", txt,
      "; coefficients:", length(cols), "; weights:",
      if(is.numeric(out)) length(out) else paste("evaluation failed:",
        if(inherits(out, "try-error")) conditionMessage(attr(out, "condition"))
        else class(out)), config=config)
  }

  names(out) <- cols

  ## all weights zero, most easily by writing a coefficient minus itself. There is
  ##   then no hypothesis, the same situation f.design_test_cols() reports as a
  ##   df_intend of zero for config$test_term:

  if(all(abs(out) < sqrt(.Machine$double.eps))) {
    f.err(caller, ": every coefficient weight in config$contrast is zero, so",
      "there is no hypothesis to test; config$contrast:", txt, config=config)
  }

  return(out)
}

## A basis for the null space of the contrast, as the columns of a matrix with one
##   column fewer than L has elements. Post-multiplying the design by it gives the
##   design of the model constrained so that the contrast is zero, which is nested
##   in the full model and one rank below it. That reduced design is what makes a
##   contrast reach the engines that test a full model against a reduced one
##   without any contrast-specific code of their own; see f.design_contrast().
##   The first column of the complete Q of L spans L, so the rest span its
##   orthogonal complement:

f.contrast_null_basis <- function(L) {
  qrl <- qr(matrix(L, ncol=1))
  return(qr.Q(qrl, complete=TRUE)[, -1, drop=F])
}

## The design matrix for config$frm over the observations of state, checked against
##   state$expression. Shared by the config$test_term and config$contrast branches
##   of f.design_test_cols():

f.design_X <- function(state, parsed, config, caller="f.design_test_cols") {

  X <- stats::model.matrix(parsed$frm, data=state$samples)

  if(nrow(X) != ncol(state$expression)) {
    f.err(caller, ": design matrix has", nrow(X),
      "rows, but state$expression has", ncol(state$expression), "columns;",
      "\n", "model.matrix() drops observations with missing covariate values",
      config=config)
  }

  return(X)
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

  ## config$contrast tests a weighted sum of coefficients within the full model
  ##   instead of a term of it, so it is a different hypothesis, derived below by
  ##   its own route. check_config() has already refused a config that sets both:

  if(f.contrast_set(config)) return(f.design_contrast(state, config))

  ## throws an informative error if config$test_term does not fit config$frm:

  parsed <- f.parse_frm(config$frm, config)
  drops <- f.normalize_terms(config)$drop_terms

  X <- f.design_X(state, parsed, config)

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

  ## nothing left in the reduced model. Two quite different causes, told apart here
  ##   because the advice differs and the wrong advice sends the reader to config$frm
  ##   when the formula is fine. If the full model has no rank either, there is nothing
  ##   to fit at all, which is what an upstream step that removed every observation looks
  ##   like from here. Otherwise config$frm suppresses the intercept and config$test_term
  ##   names every remaining term, so the test is against zero rather than against a
  ##   common mean: well defined, but on log-scale abundances every feature rejects it,
  ##   which is rarely the question. Warn rather than stop, since the test asked for is
  ##   the one performed, and it is also what filter_features_by_estimability() screens:

  if(rank_red %in% 0) {
    if(rank_all %in% 0) {
      f.msg("WARNING: f.design_test_cols: the design matrix has rank 0, so there is",
        "nothing to fit with or without config$test_term '", config$test_term, "';",
        "\n", "observations:", nrow(X), "; design columns:", ncol(X), "; config$frm:",
        deparse(parsed$frm), "; config$frm is not the likely problem here: an empty",
        "design is what dropping every observation looks like, so check the",
        "filtering steps and their thresholds", config=config)
    } else {
      f.msg("WARNING: f.design_test_cols: dropping config$test_term '",
        config$test_term, "' leaves a reduced model with no parameters, so the",
        "test is of whether the", config$test_term, "means are all zero, not of",
        "whether they differ from each other;", "\n",
        "config$frm:", deparse(parsed$frm), "; rank of the full model:", rank_all,
        "; for the usual comparison among levels, keep the intercept in config$frm",
        config=config)
    }
  }

  ## the intercept is the fitted value where every covariate is zero. For a factor
  ##   that is its reference level, which is an observed group; for a continuous
  ##   covariate it need not be anywhere near the observed range, and where zero falls
  ##   decides the answer rather than merely the labelling: testing the intercept in
  ##   ~age compares the space spanned by (1, age) against the space spanned by age
  ##   alone, and shifting age leaves the first unchanged while changing the second.
  ##   Centering the covariate in state$samples makes the intercept the fitted value at
  ##   its mean, which is usually the intended quantity:

  if("1" %in% drops) {
    types <- f.covariate_types(state, config)
    nums <- names(types)[types %in% "numeric"]
    if(length(nums)) {
      f.msg("WARNING: f.design_test_cols: config$test_term names the intercept",
        "('1'), and config$frm has the continuous covariate(s)",
        paste(paste(nums, collapse=", "), ","), "so the intercept is the fitted",
        "value where those are zero, which need not be near any observation;", "\n",
        "the test therefore depends on where zero falls for them: center the",
        "covariate(s) in state$samples to test at their mean instead", config=config)
    }
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

  ## the reduced model, returned explicitly rather than left to each engine to form
  ##   by dropping columns, so that the engines which compare a full model against a
  ##   reduced one need no contrast-specific code: see f.design_contrast(), where
  ##   this is not a subset of the columns of X:

  out <- list(
    parsed=parsed,
    drops=drops,
    X=X,
    X_red=X[, -cols_test, drop=F],
    cols_test=cols_test,
    rank_all=rank_all,
    df_intend=df_intend,
    contrast=NULL
  )

  return(out)
}

## The design matrix for config$frm and the weighted sum of its coefficients that
##   config$contrast tests, in the shape f.design_test_cols() returns for
##   config$test_term, so that the estimability filter and the hypothesis tests take
##   the two the same way. The test is of whether that weighted sum is zero, which
##   is one degree of freedom however many coefficients carry a non-zero weight. That
##   used to be the only way a capped engine could reach a multi-level factor; no engine
##   is capped any more, so a contrast is now a different hypothesis rather than a way
##   around one, and it stays one degree of freedom for the limma-family engines because
##   limma::contrasts.fit() leaves them a single coefficient.
##   $X_red is the design of the model constrained so that the contrast is zero,
##   which is nested in the full model and one rank below it, so comparing the two
##   is exactly the test of the contrast. It is a re-parameterization rather than a
##   subset of the columns of X, which is why f.design_test_cols() returns it
##   explicitly in both cases:

f.design_contrast <- function(state, config) {

  parsed <- f.parse_frm(config$frm, config)
  X <- f.design_X(state, parsed, config, caller="f.design_contrast")

  L <- f.contrast_vector(config, colnames(X), caller="f.design_contrast")
  cols_test <- which(abs(L) >= sqrt(.Machine$double.eps))

  ## the contrast has to be estimable over all observations, which it is exactly
  ##   when it lies in the row space of the design: otherwise no linear combination
  ##   of the fitted values estimates it, and every engine would report either a
  ##   silently aliased coefficient or an error of its own naming neither
  ##   config$contrast nor config$frm. Reached by weighting a coefficient that
  ##   config$frm codes as aliased with another. Per feature estimability is a
  ##   separate matter, screened by filter_features_by_estimability() from the same
  ##   pair of designs:

  rank_all <- f.design_rank(X)
  Lr <- L / sqrt(sum(L^2))

  if(f.design_rank(rbind(X, matrix(Lr, nrow=1))) != rank_all) {
    f.err("f.design_contrast: config$contrast", config$contrast, "is not",
      "estimable: it does not lie in the row space of the design matrix for",
      "config$frm, so no combination of the fitted values estimates it;", "\n",
      "config$frm:", deparse(parsed$frm), "; design columns:", ncol(X),
      "; rank:", rank_all, "\n", "coefficients weighted:",
      paste(names(L)[cols_test], collapse=", "), config=config)
  }

  X_red <- X %*% f.contrast_null_basis(L)
  colnames(X_red) <- paste0("h0red", seq_len(ncol(X_red)))
  df_intend <- rank_all - f.design_rank(X_red)

  ## the constrained model is one rank below the full model whenever the contrast is
  ##   estimable, which the check above has established, so this cannot fire; kept
  ##   because it is the assumption every engine below relies on, and a silent
  ##   failure of it would be reported as a test of the wrong number of degrees of
  ##   freedom rather than as an error:

  if(df_intend != 1) {
    f.err("f.design_contrast: the model constrained so that config$contrast",
      config$contrast, "is zero is", df_intend, "degrees of freedom below the",
      "full model rather than 1;", "\n", "config$frm:", deparse(parsed$frm),
      "; design rank:", rank_all, config=config)
  }

  f.warn_contrast_marginality(parsed, X, L, cols_test, state, config)

  out <- list(
    parsed=parsed,
    drops=character(0),
    X=X,
    X_red=X_red,
    cols_test=cols_test,
    rank_all=rank_all,
    df_intend=df_intend,
    contrast=L
  )

  return(out)
}

## Warn when config$contrast weights a coefficient of a term that a higher-order
##   term of config$frm also contains. Testing a term obeys marginality: dropping
##   'grp' from ~grp*sex drops the interaction too, so the test covers every term
##   containing grp. A contrast cannot, being a statement about named coefficients
##   inside one model, so 'grpb - grpc' in ~grp*sex compares those levels at the
##   reference level of sex alone rather than averaged over it. That is a legitimate
##   question and is the one asked for, so it is a warning rather than an error, but
##   it is rarely the question intended, and the coefficients involved do not say so
##   on their own. Same shape of problem, and the same remedy, as the intercept
##   warning in f.design_test_cols(): the answer depends on where the variables not
##   under test are held:

f.warn_contrast_marginality <- function(parsed, X, L, cols_test, state, config) {

  asgn <- attr(X, "assign")
  if(is.null(asgn)) return(invisible(NULL))

  ## the terms of config$frm carrying a weighted coefficient, and any term of
  ##   config$frm strictly containing one of them:

  trms <- unique(asgn[cols_test])
  labs <- parsed$labels[trms[trms > 0]]
  if(!length(labs)) return(invisible(NULL))

  higher <- character(0)

  for(lab in labs) {
    vars1 <- unlist(strsplit(lab, ":", fixed=T))
    for(lab2 in setdiff(parsed$labels, lab)) {
      vars2 <- unlist(strsplit(lab2, ":", fixed=T))
      if(all(vars1 %in% vars2) && length(vars2) > length(vars1)) {
        higher <- c(higher, lab2)
      }
    }
  }

  if(!length(higher)) return(invisible(NULL))

  higher <- unique(higher)
  vars_held <- setdiff(unlist(strsplit(higher, ":", fixed=T)),
    unlist(strsplit(labs, ":", fixed=T)))

  types <- f.covariate_types(state, config)
  held <- vapply(
    vars_held,
    function(nom) {
      if(!(nom %in% names(types))) return(nom)
      if(types[[nom]] %in% "factor") {
        lvl <- levels(f.relevel_covariate(state$samples[[nom]], nom, config,
          "f.warn_contrast_marginality"))[1]
        return(paste0(nom, " = ", lvl))
      }
      return(paste0(nom, " = 0"))
    },
    character(1)
  )

  f.msg("WARNING: f.design_contrast: config$contrast", config$contrast,
    "weights coefficients of the term(s)", paste(paste(labs, collapse=", "), ","),
    "which the higher-order term(s)", paste(paste(higher, collapse=", "), ""),
    "of config$frm also contain;", "\n",
    "the contrast therefore compares them with", paste(held, collapse=", "),
    "rather than averaged over", paste(paste(vars_held, collapse=", "), ","),
    "unlike config$test_term, which by marginality would test every term",
    "containing them;", "\n", "config$frm:", deparse(parsed$frm),
    "; to average instead, weight the higher-order coefficients too", config=config)

  return(invisible(NULL))
}

## f.test_max_cols() was here, returning how many design matrix columns a test method
##   could test at once, along with f.design_test_cols_max(), which wrapped
##   f.design_test_cols() so that a method which could not express the hypothesis
##   config$test_term implied said so rather than quietly testing a narrower one. Both
##   are gone, because no method is bounded that way any more.
##   The limma-family methods were always unlimited, taking a vector of coefficients and
##   returning an F over all of them, and proDA is unlimited by a different route:
##   proDA::test_diff() takes either one contrast or a reduced model, and test_proda()
##   hands it the reduced model f.design_test_cols() built whenever more than one column
##   carries the test, which is a likelihood ratio test over all of them.
##   The other two were bounded only at their APIs, not in what their fits could support.
##   msqrob2::hypothesisTest() loops over the columns of the contrast and returns one
##   table per column rather than a joint test over several, while the fit carries
##   everything a joint test needs; test_msqrob() computes that test from the fitted
##   models, see f.msqrob_wald(). DEqMS moderates one coefficient's t-statistic and has
##   no F-analogue anywhere, but the moderation is a variance prior: its
##   DEqMS::spectraCounteBayes() returns a per-gene posterior variance and a prior
##   degrees of freedom, neither of which mentions a coefficient, and coef_col enters
##   only where sca.t is formed from them. test_deqms() computes the joint test from
##   those, see f.deqms_moderated_f(). What was said here before, that "the limit is in
##   the moderation itself", was wrong on that second point.
##   So config$test_term and config$contrast now reach every method, and the only
##   remaining refusals are properties of an engine's model rather than of its API: see
##   the capability table in the tests directory.

## Whether a test method takes feature level input and reports gene level results,
##   whatever level state$expression is at. test_deqms() aggregates with
##   combine_features() and test_msqrob() with QFeatures::aggregateFeatures(), both
##   internally, so their results have one row per gene even when handed precursors;
##   "prolfqua_lmer" reaches the same place by modelling instead of aggregating, fitting
##   one mixed model per gene with a random effect for the feature, and "msqrob_agg" the
##   same way through msqrob2::msqrobAggregate(), which aggregates only to carry the
##   result and fits on the un-aggregated assay; every other method
##   is row-wise on state$expression and reports one row per row of it. Written once
##   here because four places need the same answer: test() (for the row ids and the
##   feature metadata it reports), f.feature_means() (for the level the average
##   expression is over), the f.format_*() functions (for which column of an engine's
##   table holds the id) and tune() (for whether to aggregate before testing):

f.gene_level_method <- function(method) {
  return(method %in% c("deqms", "msqrob", "prolfqua_lmer", "msqrob_agg"))
}

## Which column of state$features identifies the rows a test method returns.
##   config$feat_col is by definition the column matching rownames(state$expression)
##   at every point in the pipeline, so it is the answer for a row-wise method whether
##   state is at precursor level or has already been through combine_features(), which
##   sets config$feat_col to config$gene_id_col when it aggregates. The two gene level
##   methods aggregate for themselves, so they are keyed by the gene id either way:

f.test_id_col <- function(method, config) {
  if(f.gene_level_method(method)) return(config$gene_id_col)
  return(config$feat_col)
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

## Check the values of the covariates in config$frm. Missing, blank, non-finite,
##   and constant covariates are errors: they cannot be fit, and letting them
##   through means model.matrix() silently drops observations (so that the
##   number of observations differs between features) or yields a
##   rank-deficient design. A numeric covariate with few distinct values is
##   often a miscoded factor, so is warned about; threshold is
##   config$n_distinct_numeric_warn. Optional types from f.covariate_types().
##   caller names the function to blame in the messages, since this is called from
##   more than one entry point; warn_distinct=FALSE suppresses the distinct-value
##   warning for callers that run after initialize() has already issued it:

f.check_covariate_values <- function(state, config, types=NULL,
    caller="f.check_covariate_values", warn_distinct=TRUE) {

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
      f.err(caller, ": covariate", nom, "has missing values,",
        "which are not supported;", "\n", "n missing:", sum(i),
        "; first offending observations:", utils::head(obs[i], 10),
        config=config)
    }

    if(types[nom] %in% "numeric" && !all(is.finite(v))) {
      i <- !is.finite(v)
      f.err(caller, ": covariate", nom, "has non-finite",
        "values;", "\n", "n non-finite:", sum(i),
        "; first offending observations:", utils::head(obs[i], 10),
        config=config)
    }

    ## a blank is a missing value that does not look like one: utils::read.table(),
    ##   which read_data() uses, reads an empty field in a character column as "",
    ##   and only the strings in na.strings (default "NA") as NA. An empty cell in
    ##   the samples file therefore arrives here as a value rather than as a gap,
    ##   and would be carried into the design as a factor level of its own, silently
    ##   adding a group made of the observations whose annotation is missing. An
    ##   empty field in a numeric column does become NA, so this only applies to the
    ##   covariates that carry text:

    if(types[nom] %in% "factor") {
      i <- !nzchar(trimws(as.character(v)))
      if(any(i)) {
        f.err(caller, ": covariate", nom, "has blank values, which are treated",
          "as missing and are not supported;", "\n", "n blank:", sum(i),
          "; first offending observations:", utils::head(obs[i], 10), "\n",
          "an empty cell in the samples file is read as an empty string rather",
          "than as NA, so it would otherwise become a factor level of its own",
          config=config)
      }
    }

    lvls <- unique(v)

    if(length(lvls) %in% 1) {
      f.err(caller, ": covariate", nom, "is constant, with",
        "single distinct value:", utils::head(as.character(lvls), 1), ";", "\n",
        "it cannot be fit; drop it from config$frm", config=config)
    }

    if(warn_distinct && types[nom] %in% "numeric" && length(lvls) <= cutoff) {
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

## Resolve whether the variance prior of the moderation is fitted against mean feature
##   intensity, for a function that takes a trend argument. The argument wins when given,
##   config$test_trend answers otherwise, and FALSE when neither says anything, which is
##   what new_config() ships. Unlike f.is_log_transformed(), disagreement between the two
##   is not an error: the scale of the data is a fact and the two cannot both be right,
##   whereas this is a preference about how to fit a prior, so a caller passing one is
##   overriding the configuration on purpose. Which methods honor it, and what the ones
##   that cannot do instead, is test()'s business; see f.trend_methods():

f.is_trend <- function(trend, config, fn_name="f.is_trend") {

  from_config <- config$test_trend

  if(!is.null(from_config) &&
      !(is.logical(from_config) && length(from_config) %in% 1 &&
        !is.na(from_config))) {
    f.err(fn_name, ": config$test_trend is not TRUE or FALSE; value:",
      from_config, config=config)
  }

  ## "" is accepted as unset for callers that pass an empty character:

  unset <- is.null(trend) || (is.character(trend) && all(trend %in% ""))

  if(unset) {
    if(is.null(from_config)) return(FALSE)
    return(from_config)
  }

  if(!(is.logical(trend) && length(trend) %in% 1 && !is.na(trend))) {
    f.err(fn_name, ": trend is not TRUE or FALSE; value:", trend,
      "; typeof:", typeof(trend), config=config)
  }

  return(trend)
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

