## Flatten arguments destined for cat() into character scalars, one per argument:

f.cat_args <- function(...) {
  args <- list(...)
  out <- lapply(args, function(arg) {
    if(length(arg) %in% 0) return(paste(deparse(arg), collapse=" "))
    if(is.numeric(arg)) return(paste(format(arg, trim=T), collapse=" "))
    if(is.character(arg) || is.logical(arg) || is.factor(arg)) {
      return(paste(arg, collapse=" "))
    }
    paste(deparse(arg), collapse=" ")
  })
  return(unlist(out))
}

f.log_notice <- local({
  seen <- character(0)
  function(file) {
    if(file %in% seen) return(invisible(NULL))
    seen <<- c(seen, file)
    cat("WARNING: h0testr: cannot write to config$log_file '", file,
      "'; the messages below are going to the console instead\n", sep="")
    utils::flush.console()
    invisible(NULL)
  }
})

## Write one already-flattened message to config$log_file, falling back to 
##   console when file cannot be opened:

f.cat_log <- function(..., config) {

  if(is.null(config$log_file)) config$log_file <- ""

  if(config$log_file %in% "") {
    cat(...)
    utils::flush.console()
    return(invisible(NULL))
  }

  ok <- tryCatch({
      cat(..., file=config$log_file, append=T)
      TRUE
    },
    error=function(e) FALSE,
    warning=function(w) FALSE
  )

  if(!ok) {
    f.log_notice(config$log_file)
    cat(...)
  }

  utils::flush.console()
  invisible(NULL)
}

f.msg <- function(..., config) {
  f.cat_log(f.cat_args(...), "\n", config=config)
}

f.log <- function(..., config) {
  f.cat_log(f.cat_args(...), "at:", format(Sys.time(), format='%Y%m%d%H%M%S'),
    "\n", config=config)
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

f.pkg_install_cmd <- function(pkg) {

  bioc <- c("DEqMS", "edgeR", "impute", "limma", "MsCoreUtils", "msqrob2",
    "pcaMethods", "proDA", "QFeatures", "SummarizedExperiment", "vsn")
  gh <- c(prolfqua="wolski/prolfqua")

  if(pkg %in% names(gh)) {
    return(paste0('remotes::install_github("', gh[[pkg]], '")'))
  }
  if(pkg %in% bioc) return(paste0('BiocManager::install("', pkg, '")'))
  return(paste0('install.packages("', pkg, '")'))
}

f.need_pkgs <- function(pkgs, who, config) {

  if(length(pkgs) < 1) return(invisible(TRUE))
  i <- !vapply(pkgs, requireNamespace, logical(1), quietly=TRUE)
  if(!any(i)) return(invisible(TRUE))

  cmds <- vapply(pkgs[i], f.pkg_install_cmd, character(1))
  f.err(who, "needs the package", paste(pkgs[i], collapse=", "),
    "which could not be loaded; install with:", "\n  ",
    paste(cmds, collapse="\n   "), config=config)
}

f.log_obj <- function(obj, config) {

  if(is.null(config$log_file)) config$log_file <- ""

  if(config$log_file %in% "") {
    print(obj)
    utils::flush.console()
    return(invisible(NULL))
  }

  ## same fallback as f.cat_log() for printed objects (quantiles, tables):

  ok <- tryCatch({
      utils::capture.output(obj, file=config$log_file, append=T)
      TRUE
    },
    error=function(e) FALSE,
    warning=function(w) FALSE
  )

  if(!ok) {
    f.log_notice(config$log_file)
    print(obj)
  }

  utils::flush.console()
  invisible(NULL)
}

## Write one table to file_out as a tab delimited file:

f.save_tsv <- function(dat, file_out, config, row.names=T, col.names=T) {

  if(!(length(file_out) %in% 1 && is.character(file_out) && !is.na(file_out) &&
      nzchar(file_out))) {
    f.err("f.save_tsv: file_out has to be a single non-empty file name;", "\n",
      " file_out:", file_out, "; class:", class(file_out), "; length:",
      length(file_out), config=config)
  }

  file_tmp <- paste0(file_out, ".tmp", Sys.getpid())

  ## condition is returned rather than acted on inside handler:

  cond <- tryCatch({
      utils::write.table(dat, file=file_tmp, quote=F, sep="\t",
        row.names=row.names, col.names=col.names)
      NULL
    },
    error=function(msg) list(kind="error", msg=conditionMessage(msg)),
    warning=function(msg) list(kind="warning", msg=conditionMessage(msg))
  )

  if(is.null(cond) && !file.exists(file_tmp)) {
    cond <- list(kind="silent failure", msg="no file was written")
  }

  if(!is.null(cond)) {
    left <- file.exists(file_tmp) && !(unlink(file_tmp) %in% 0)
    f.err("f.save_tsv: write.table()", paste0(cond$kind, ":"), "writing to",
      paste0(file_out, ":"), cond$msg, "\n",
      " nothing was written to", paste0(file_out, ";"), "the table was being",
      "written to", paste0(file_tmp, ","), "which is moved into place only once",
      "it is complete, so any file already at", file_out, "is as it was",
      if(left) {
        paste("\n  the temporary file could not be removed and is still there:",
          file_tmp)
      } else "", config=config)
  }

  ## file.rename() replaces an existing destination on platforms this runs
  ##   on, but not on all of them. It cannot cross a file system:

  if(!file.rename(file_tmp, file_out)) {

    if(!file.copy(file_tmp, file_out, overwrite=TRUE)) {
      f.err("f.save_tsv: the table was written but could not be moved into",
        "place;", "\n", " from:", file_tmp, "\n", " to:", file_out, "\n",
        "the finished table is at the first of those; move it there by hand,",
        "or correct the destination and run the step again", config=config)
    }

    unlink(file_tmp)
  }

  return(invisible(file_out))
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

## Operators supported in config$frm:

f.frm_ops <- c("+", ":", "*", "-")

## Recursively check one expression from the right-hand side of a formula:

f.check_frm_expr <- function(expr, config, top=NULL) {

  if(is.null(top)) top <- expr

  if(is.symbol(expr)) {
    if(as.character(expr) %in% ".") {
      f.err("f.check_frm_expr: '.' not supported in formula; frm rhs:",
        deparse(top), config=config)
    }
    return(invisible(TRUE))
  }

  ## only 0 and 1 are meaningful (intercept):

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

  ## what is called is normally a symbol, an operator name among f.frm_ops:

  if(!is.symbol(expr[[1]])) {
    f.err("f.check_frm_expr: unsupported construct", deparse(expr),
      "in formula: what is called is itself an expression rather than a variable",
      "name or one of the supported operators:", paste(f.frm_ops, collapse=" "),
      "; frm rhs:", deparse(top), config=config)
  }

  op <- as.character(expr[[1]])

  if(!(op %in% f.frm_ops)) {
    f.err("f.check_frm_expr: operator or function '", op,
      "' not supported in formula; supported operators:",
      paste(f.frm_ops, collapse=" "),
      "; variables must be bare column names (no transformations, like",
      "log() or I()); frm rhs:", deparse(top), config=config)
  }

  ## operator called with no operands, e.g. ~`+`(), which is call of length 1:

  if(length(expr) < 2) {
    f.err("f.check_frm_expr: operator '", op, "' has no operands in formula;",
      "got", deparse(expr), "; frm rhs:", deparse(top), config=config)
  }

  ## '-' only for dropping intercept, so what is removed must be 1 or 0;
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

## Sort variables within each interaction term label:

f.canon_label <- function(labels) {
  if(length(labels) %in% 0) return(character(0))
  out <- vapply(
    labels,
    function(lbl) paste(sort(unlist(strsplit(lbl, ":", fixed=T))), collapse=":"),
    character(1)
  )
  return(unname(out))
}

## Parse and check model formula. Single place where config$frm is
##   interpreted. Throws an informative error on unsupported constructs.
##   Returns list:
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

  ## dependent variable always expression values of one feature:

  two_sided <- length(frm) %in% 3
  if(two_sided) frm <- frm[-2]

  ## '|' otherwise reported as unsupported:

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

## TRUE iff config$contrast asks for a contrast. "" (default) means none:

f.contrast_set <- function(config) {
  if(length(config$contrast) != 1) return(FALSE)
  if(is.na(config$contrast)) return(FALSE)
  return(nzchar(trimws(config$contrast)))
}

## What run is testing, for progress messages: 

f.test_label <- function(design, config) {
  if(!is.null(design$contrast)) return(paste("contrast:", trimws(config$contrast)))
  return(paste("test_term:", config$test_term))
}

## Operators supported in config$contrast:

f.contrast_ops <- c("+", "-", "*", "/", "(")

## TRUE iff expr mentions any coefficient name; helper for f.check_contrast_expr():

f.contrast_has_name <- function(expr) {
  if(is.symbol(expr)) return(TRUE)
  if(!is.call(expr) || length(expr) < 2) return(FALSE)
  for(idx in 2:length(expr)) {
    if(f.contrast_has_name(expr[[idx]])) return(TRUE)
  }
  return(FALSE)
}

## Recursively check one expression from config$contrast against cols, the
##   coefficient names of the design matrix for config$frm:

f.check_contrast_expr <- function(expr, cols, config, top=NULL) {

  if(is.null(top)) top <- expr

  ## terminal symbol: a coefficient name, which must be a column of the design:

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

  ## terminal constant: any finite number, unlike in a formula:

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

  ## same as in f.check_frm_expr():

  if(!is.symbol(expr[[1]])) {
    f.err("f.check_contrast_expr: unsupported construct", deparse(expr),
      "in config$contrast: what is called is itself an expression rather than a",
      "coefficient name or one of the supported operators:",
      paste(f.contrast_ops, collapse=" "), "; config$contrast:", deparse(top),
      config=config)
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

  ## operator called without operands e.g. '`+`()'; call of length 1:

  if(length(expr) < 2) {
    f.err("f.check_contrast_expr: operator '", op, "' has no operands in",
      "config$contrast; got", deparse(expr), "; config$contrast:", deparse(top),
      config=config)
  }

  ## '*' and '/' scale a contrast by a constant:

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
##   design matrix for config$frm:

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

  ## same bindings again, in environment where '+' and '-' add magnitudes
  ##   instead of cancelling:
  env_abs <- new.env(parent=baseenv())

  for(idx in seq_along(cols)) {
    wts <- rep(0, length(cols))
    wts[idx] <- 1
    assign(cols[idx], wts, envir=env)
    assign(cols[idx], wts, envir=env_abs)
  }

  assign("+", function(a, b) if(missing(b)) abs(a) else abs(a) + abs(b),
    envir=env_abs)
  assign("-", function(a, b) if(missing(b)) abs(a) else abs(a) + abs(b),
    envir=env_abs)
  assign("*", function(a, b) abs(a) * abs(b), envir=env_abs)
  assign("/", function(a, b) abs(a) / abs(b), envir=env_abs)

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

  ## whether there is a hypothesis at all. The hypothesis is that the weighted sum is zero, 
  ##   which is invariant under rescaling the weights, so '(grpb - grpc)/1e9' asks exactly
  ##   what 'grpb - grpc' asks while carrying weights of 1e-9. A contrast written
  ##   so that it cancels, on the other hand, such as
  ##   '(grpb + grpc)/3 - grpb/3 - grpc/3', is left with weights of about 1e-17
  ##   whose direction is rounding noise rather than a hypothesis. An absolute
  ##   threshold cannot separate those two:

  no_cancel <- try(eval(expr, envir=env_abs), silent=T)

  usable <- is.numeric(no_cancel) && length(no_cancel) %in% length(cols) &&
    all(is.finite(no_cancel))

  scl <- max(abs(out))
  floor_wt <- if(usable) max(no_cancel) * 128 * .Machine$double.eps else 0

  if(scl <= 0) {
    f.err(caller, ": every coefficient weight in config$contrast is zero, so",
      "there is no hypothesis to test; config$contrast:", txt, config=config)
  }

  if(scl <= floor_wt) {
    f.err(caller, ": the coefficient weights in config$contrast cancel to",
      "rounding error, so there is no hypothesis to test: their direction is",
      "what is left of the cancellation rather than a contrast;", "\n",
      "config$contrast:", txt, "; largest weight:", format(scl, digits=3),
      "; largest weight without the cancellation:",
      format(max(no_cancel), digits=3), "\n",
      "  a contrast that is a rescaling of another, such as '(grpb - grpc)/1e9',",
      "is not this: its weights are small but they do not cancel", config=config)
  }

  return(out)
}

## Basis for null space of contrast, as columns of a matrix with one
##   column fewer than L has elements:

f.contrast_null_basis <- function(L) {
  qrl <- qr(matrix(L, ncol=1))
  return(qr.Q(qrl, complete=TRUE)[, -1, drop=F])
}

## Design matrix for config$frm over observations of state, checked against
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

## A design with fewer independent columns than columns cannot support a test of
##   one of its terms: 

f.design_check_rank <- function(X, rank_all, parsed, config,
    caller="f.design_test_cols") {

  if(!(rank_all > 0 && rank_all < ncol(X))) return(invisible(NULL))

  ## qr() pivots columns it could not use to the end:

  qrx <- qr(X)
  aliased <- colnames(X)[qrx$pivot[-seq_len(qrx$rank)]]
  empty <- colnames(X)[colSums(X != 0) %in% 0]

  f.err(caller, ": the design matrix for config$frm is rank deficient: it has",
    ncol(X), "columns but only", rank_all, "of them are linearly independent, so",
    "the coefficient(s)", paste(aliased, collapse=" "), "cannot be estimated;",
    "\n", "config$frm:", deparse(parsed$frm), "; design columns:",
    paste(colnames(X), collapse=" "), "\n",
    if(length(empty)) {
      paste("column(s)", paste(empty, collapse=" "), "are all zeros, so no",
        "observation is at the level they code: drop the unused level(s) from",
        "the covariate in state$samples (droplevels()), or drop the term.",
        "h0testr::test_h0() and the other fitting entry points drop them for you;",
        "a design built directly from state$samples does not")
    } else {
      paste("no column is all zeros, so two terms of config$frm code the same",
        "grouping: drop the redundant one")
    }, config=config)
}

## Rank of a design matrix:

f.design_rank <- function(mat) {
  if(nrow(mat) %in% 0 || ncol(mat) %in% 0) return(0L)
  return(qr(mat)$rank)
}

## Helper for f.normalize_terms():

f.formula2terms <- function(config) {

  frm <- config$frm

  if(is.null(frm) || all(as.character(frm) %in% "")) {
    f.err("f.formula2terms: config$frm empty or undefined.", config=config)
  }

  parsed <- f.parse_frm(frm, config)

  ## make intercept explicit, so can be dropped or kept downstream:

  if(parsed$intercept %in% 1) {
    terms <- c("1", parsed$labels)
  } else {
    terms <- c("0", parsed$labels)
  }

  return(terms)
}

## Helper for f.normalize_terms():

f.test_term_drops <- function(parsed, test_term, config) {

  if(test_term %in% "1") {
    if(parsed$intercept %in% 0) {
      f.err("f.test_term_drops: test_term is intercept ('1'), but frm has no",
        "intercept; frm:", parsed$frm, config=config)
    }
    return("1")
  }

  vars <- rownames(parsed$factors)

  ## test_term names variable: drop every term it is found in:

  if(test_term %in% vars) {
    drops <- parsed$labels[parsed$factors[test_term, ] > 0]
    return(drops)
  }

  ## test_term names a term of model: drop it alone, but only if no
  ##   higher-order term contains it:

  if(test_term %in% parsed$labels) {
    vars_test <- vars[parsed$factors[, test_term] > 0]
    for(lbl in setdiff(parsed$labels, test_term)) {
      vars_lbl <- vars[parsed$factors[, lbl] > 0]
      if(all(vars_test %in% vars_lbl)) {
        f.err("f.test_term_drops: test_term '", test_term,
          "' is contained in higher-order term '", lbl,
          "' of frm, so cannot be tested on its own; test the variables of",
          "'", test_term, "' individually, or test '", lbl, "' instead",
          config=config)
      }
    }
    return(test_term)
  }

  f.err("f.test_term_drops: test_term '", test_term,
    "' is neither a variable nor a term of frm; variables:",
    paste(vars, collapse=" "), "; terms:",
    paste(parsed$labels, collapse=" "), config=config)
}

## Helper for f.design_test_cols():

f.normalize_terms <- function(config) {

  test_term <- config$test_term
  frm <- config$frm
  
  if(is.null(test_term) || test_term %in% "") {
    f.err("f.normalize_terms: config$test_term empty or undefined.", config=config)
  }
  
  if(is.null(frm) || all(as.character(frm) %in% "")) {
    f.err("f.normalize_terms: config$frm empty or undefined.", config=config)
  }
  
  if(length(test_term) != 1) {
    f.err("f.normalize_terms: length(test_term) != 1; test_term: '", 
      paste(test_term), "'", config=config)
  }
  
  if(grepl("[\\*\\-\\|\\(\\)\\^/]", test_term)) {
    f.err("f.normalize_terms: cannot handle '*', '-', '|', '^', '/', '(', or ')'",
      " in test_term: '", test_term, "'; test_term must name a single variable",
      " or a single interaction term, e.g. 'age' or 'age:sex'", config=config)
  }

  test_term <- gsub("[[:space:]]", "", test_term)

  if(test_term %in% "0") {
    f.err("f.normalize_terms: invalid test term: '", test_term, "'",
      config=config)
  }

  test_term <- f.canon_label(test_term)

  parsed <- f.parse_frm(frm, config)
  frm_terms <- f.formula2terms(config)    ## returns character vector

  ## terms of the full model to be dropped to form the reduced model; throws
  ##   error if test_term is not compatible with frm:

  drop_terms <- f.test_term_drops(parsed, test_term, config)

  return(list(test_term=test_term, frm_terms=frm_terms, drop_terms=drop_terms))
}

## Design matrix for config$frm, columns carry the test of
##   config$test_term, and the degrees of freedom of that test measured over
##   observations:

f.design_test_cols <- function(state, config) {

  ## config$contrast tests weighted sum of coefficients within full model
  ##   instead of a term:

  if(f.contrast_set(config)) return(f.design_contrast(state, config))

  ## throws informative error if config$test_term doesn't fit config$frm:

  parsed <- f.parse_frm(config$frm, config)
  drops <- f.normalize_terms(config)$drop_terms

  X <- f.design_X(state, parsed, config)

  ## attr(X, 'assign') indexes term labels in order; 0 for the
  ##   intercept; parsed$labels is same labels, but canonicalized:

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

  ## before df_intend arithmetic below:

  f.design_check_rank(X, rank_all, parsed, config, caller="f.design_test_cols")

  rank_red <- f.design_rank(X[, -cols_test, drop=F])
  df_intend <- rank_all - rank_red

  ## nothing left in reduced model:

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

  ## intercept is fitted value where every covariate is zero:

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
    hint <- if(rank_all %in% 0) {
        paste("the full model has no rank either, so there is nothing to fit",
          "with or without config$test_term: see the warning above, and check",
          "the filtering steps and the covariate values rather than config$frm")
      } else if("1" %in% drops) {
        paste("testing the intercept ('1') is only meaningful when no factor in",
          "config$frm is coded to full rank in the reduced model")
      } else {
        paste("the term(s) dropped add no rank to what the remaining terms of",
          "config$frm already span, so they code a grouping the rest imply:",
          "drop the redundant term, or test one that is not implied.")
      }
    f.err("f.design_test_cols: dropping config$test_term '", config$test_term,
      "' leaves a reduced model spanning the same space as the full model,",
      "so there is no hypothesis to test;", "\n",
      "config$frm:", deparse(parsed$frm), "; terms dropped to form the",
      "reduced model:", paste(drops, collapse=" "), "; design columns:",
      ncol(X), "; rank:", rank_all, "\n", hint, config=config)
  }

  ## reduced model:

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

## Design matrix for config$frm and the weighted sum of its coefficients that
##   config$contrast tests; in shape f.design_test_cols() returns for
##   config$test_term, so estimability filter and the hypothesis tests take
##   the two the same way. Test is whether weighted sum is zero:

f.design_contrast <- function(state, config) {

  parsed <- f.parse_frm(config$frm, config)
  X <- f.design_X(state, parsed, config, caller="f.design_contrast")

  L <- f.contrast_vector(config, colnames(X), caller="f.design_contrast")

  ## which coefficients contrast weights:
  cols_test <- which(abs(L) >= max(abs(L)) * sqrt(.Machine$double.eps))

  ## contrast must be estimable over all observations; true when it within 
  ##   row space of design: 

  ## normalized by  largest weight, rather than 2-norm, before rank
  ##   test and null space basis below: 

  rank_all <- f.design_rank(X)
  Lr <- L / max(abs(L))

  if(f.design_rank(rbind(X, matrix(Lr, nrow=1))) != rank_all) {
    f.err("f.design_contrast: config$contrast", config$contrast, "is not",
      "estimable: it does not lie in the row space of the design matrix for",
      "config$frm, so no combination of the fitted values estimates it;", "\n",
      "config$frm:", deparse(parsed$frm), "; design columns:", ncol(X),
      "; rank:", rank_all, "\n", "coefficients weighted:",
      paste(names(L)[cols_test], collapse=", "), config=config)
  }

  X_red <- X %*% f.contrast_null_basis(Lr)
  colnames(X_red) <- paste0("h0red", seq_len(ncol(X_red)))
  df_intend <- rank_all - f.design_rank(X_red)

  ## the constrained model is one rank below the full model whenever contrast is
  ##   estimable, which the check above has established; so this never trigger:

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
##   term of config$frm also contains:

f.warn_contrast_marginality <- function(parsed, X, L, cols_test, state, config) {

  asgn <- attr(X, "assign")
  if(is.null(asgn)) return(invisible(NULL))

  ## terms of config$frm carrying a weighted coefficient, and any term of
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

## Whether test method takes feature level input and reports gene level results:
f.gene_level_method <- function(method) {
  return(method %in% c("deqms", "msqrob", "prolfqua_lmer", "msqrob_agg"))
}

## Which column of state$features identifies rows a test method returns:

f.test_id_col <- function(method, config) {
  if(f.gene_level_method(method)) return(config$gene_id_col)
  return(config$feat_col)
}

## Classify each variable in config$frm as "factor" or "numeric" (continuous):

f.covariate_types <- function(state, config) {

  vars <- sort(unique(f.parse_frm(config$frm, config)$vars))

  if(!all(vars %in% names(state$samples))) {
    f.err("f.covariate_types: !all(vars %in% names(state$samples)); vars:",
      vars, "; names(state$samples):", names(state$samples), config=config)
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

  ## cache as cross-check, not override:

  cached <- config$covariate_types

  if(!is.null(cached) && is.character(cached) && !is.null(names(cached))) {

    shared <- intersect(names(cached), vars)
    ok <- !is.na(cached[shared]) & cached[shared] %in% c("factor", "numeric") &
      cached[shared] == out[shared]
    bad <- shared[!ok]

    if(length(bad)) {
      f.err("f.covariate_types: config$covariate_types disagrees with",
        "state$samples about the type of", paste(bad, collapse=" "), ":", "\n",
        " config$covariate_types:",
        paste(bad, cached[bad], sep="=", collapse=" "), "\n",
        " from state$samples and config$reference_levels:",
        paste(bad, out[bad], sep="=", collapse=" "), "\n",
        "config$covariate_types is set by init_state() for the state it was",
        "run on, so a disagreement means this config is being used with other",
        "data, or was edited by hand;", "\n",
        "set config$covariate_types to NULL to have it derived from",
        "state$samples, or, if the covariate is meant to be categorical, name",
        "it in config$reference_levels instead", config=config)
    }
  }

  return(out)
}

## Level ordering of one factor covariate, resolved from config:

f.covariate_levels <- function(v, trm, config, caller="f.covariate_levels") {

  lvls <- config$factor_levels[[trm]]
  if(!is.null(lvls)) return(lvls)

  ## NOTE: config$reference_levels is atomic vector:

  if(!(trm %in% names(config$reference_levels))) return(levels(factor(v)))
  ref1 <- config$reference_levels[[trm]]

  if(!(length(ref1) %in% 1 && is.character(ref1))) {
    f.err(caller, ": config$reference_levels[[trm]] not scalar character, for",
      "trm:", trm, "\n", "value:", ref1, "; class:", class(ref1), config=config)
  }

  vals <- sort(unique(as.character(v)))

  if(!(ref1 %in% vals)) {
    f.err(caller, ": reference level", ref1, "declared in",
      "config$reference_levels is not among the values of covariate", trm, ";",
      "\n", "its values:", vals, config=config)
  }

  return(c(ref1, setdiff(vals, ref1)))
}

## Rebuild factor covariate with level ordering f.covariate_levels()
##   resolves, with declared reference level first:

f.relevel_covariate <- function(v, trm, config, caller) {

  lvls <- f.covariate_levels(v, trm, config, caller=caller)

  out <- factor(as.character(v), levels=lvls)

  if(any(is.na(out))) {
    f.err(caller, ": covariate", trm, "has values that are not among its",
      "resolved levels;", "\n", "levels:", lvls, "\n", "unmatched values:",
      utils::head(sort(unique(as.character(v)[is.na(out)])), 10), config=config)
  }

  return(out)
}

## Rebuild covariates of config$frm in state$samples whose level ordering
##   config declares:

f.relevel_state_covariates <- function(state, config,
    caller="f.relevel_state_covariates") {

  vars <- sort(unique(f.parse_frm(config$frm, config)$vars))

  for(nom in vars) {

    declared <- !is.null(config$factor_levels[[nom]]) ||
      nom %in% names(config$reference_levels)

    if(!(nom %in% names(state$samples))) {
      if(!declared) next            ## left to the checks that report it as absent
      f.err(caller, ": covariate", nom, "is not a column of state$samples;",
        "\n", "columns present:", paste(names(state$samples), collapse=", "),
        config=config)
    }

    v <- state$samples[[nom]]

    if(declared) {
      v <- f.relevel_covariate(v, nom, config, caller)
    } else if(!is.factor(v)) {
      next                          ## character and numeric carry no empty level
    }

    empty <- setdiff(levels(v), unique(as.character(v)))

    if(length(empty)) {
      f.msg("WARNING: ", caller, ": covariate ", nom, " has level(s) that no",
        "observation is at:", empty, "\n", "dropped, since model.matrix() codes",
        "such a level as a column of zeros, which leaves the design rank",
        "deficient and unfittable; the levels kept are:",
        setdiff(levels(v), empty), config=config)
      v <- droplevels(v)
    }

    state$samples[[nom]] <- v
  }

  return(state)
}

## Check values of covariates in config$frm:

f.check_covariate_values <- function(state, config, types=NULL,
    caller="f.check_covariate_values", warn_distinct=TRUE) {

  if(is.null(types)) types <- f.covariate_types(state, config)

  ## every check below is over values of a covariate:

  if(nrow(state$samples) %in% 0) {
    f.err(caller, ": state$samples has no rows, so there are no covariate",
      "values to check and nothing to fit;", "\n",
      " covariates in config$frm:", names(types), "\n",
      "an earlier step kept no observations: check the prefilter and filter",
      "thresholds, and config$frm's covariates in the samples file",
      config=config)
  }

  ## observation labels, for reporting which observations are offending:

  obs <- NULL
  for(nom in c(config$obs_col, config$obs_id_col)) {
    if(length(nom) %in% 1 && nchar(nom) > 0 && nom %in% names(state$samples)) {
      obs <- as.character(state$samples[[nom]])
      break
    }
  }
  if(is.null(obs)) obs <- colnames(state$expression)
  if(length(obs) != nrow(state$samples)) {
    obs <- as.character(seq_len(nrow(state$samples)))
  }

  cutoff <- config$n_distinct_numeric_warn
  if(is.null(cutoff)) cutoff <- 5

  for(nom in names(types)) {

    ## covariate that is not a column of state$samples: 
    if(!(nom %in% names(state$samples))) {
      f.err(caller, ": covariate", nom, "is not a column of state$samples;",
        "\n", "columns present:", paste(names(state$samples), collapse=", "),
        config=config)
    }

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

    ## blank is a missing value: 

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

## helper for f.check_state() and f.check_parameters(): 

f.check_expr_matrix <- function(state, config, fn_name) {

  if(!is.matrix(state$expression)) {
    f.err(paste0(fn_name, ":"), "!is.matrix(state$expression)", "\n",
      "class(state$expression):", class(state$expression), config=config)
  }

  if(!is.numeric(state$expression)) {
    f.err(paste0(fn_name, ":"), "state$expression is not numeric;", "\n",
      "typeof(state$expression):", typeof(state$expression),
      "; class(state$expression):", class(state$expression), "\n",
      "  every step from here on does arithmetic on it, and missingness is",
      "settled by comparing against 0 and NA, which no other type can match",
      config=config)
  }

  return(invisible(TRUE))
}

## helper for f.check_state() and f.check_parameters(): 

f.check_dimnames <- function(state, config, key, fn_name="f.check_state") {

  rows <- key %in% c("feat_col", "feat_id_col")
  meta_lab <- if(rows) "state$features" else "state$samples"
  dim_lab <- if(rows) "rownames" else "colnames"
  unit <- if(rows) "features" else "observations"
  meta <- if(rows) state$features else state$samples
  nms <- if(rows) rownames(state$expression) else colnames(state$expression)
  n_dim <- if(rows) nrow(state$expression) else ncol(state$expression)
  nom <- config[[key]]
  me <- paste0(fn_name, ":")

  if(is.null(meta)) {
    f.err(me, meta_lab, "is missing, so the", unit, "of",
      "state$expression have no metadata to be checked against", config=config)
  }

  if(!(length(nom) %in% 1 && is.character(nom) && !is.na(nom))) {
    f.err(me, paste0("config$", key), "has to be a single column name of",
      paste0(meta_lab, ";"), "\n",
      " ", paste0("config$", key, ":"), nom, "\n",
      " ", "class:", class(nom), "; length:", length(nom), "\n",
      " ", paste0("names(", meta_lab, "):"), names(meta), config=config)
  }

  ids <- meta[[nom]]

  if(is.null(ids)) {
    f.err(me, paste0("config$", key), "names no column of",
      paste0(meta_lab, ","), "so the agreement of that column with the",
      dim_lab, "of state$expression cannot be checked;", "\n",
      " ", paste0("config$", key, ":"), nom, "\n",
      " ", paste0("names(", meta_lab, "):"), names(meta), "\n",
      "  to fix, name a column of", meta_lab, "there;",
      "h0testr::init_state() sets config$feat_col and config$obs_col from",
      "config$feat_id_col and config$obs_id_col", config=config)
  }

  if(length(ids) != n_dim) {
    f.err(me, "state$expression has", n_dim, unit, "but the",
      paste0("config$", key), "column of", meta_lab, "has", length(ids),
      "ids;", "\n",
      " ", paste0("config$", key, ":"), nom, "\n",
      "  they are matched by position downstream, so there has to be one",
      "metadata row per", if(rows) "feature" else "observation",
      config=config)
  }

  if(is.null(nms)) {
    f.err(me, "state$expression has no", dim_lab, ", so the",
      unit, "it holds cannot be matched to the", paste0("config$", key),
      "column of", paste0(meta_lab, ";"), "\n",
      "  they are matched by position downstream, so agreement of the ids is",
      "the only thing that shows the order is right", config=config)
  }

  ## an NA on either side compares to NA rather than to FALSE:

  if(any(is.na(nms)) || any(is.na(ids))) {
    i <- which(is.na(nms) | is.na(ids))
    f.err(me, "the", paste0("config$", key), "column of", meta_lab, "or the",
      dim_lab, "of state$expression has missing ids;", length(i), "of",
      length(nms), "are NA, first at index", paste0(i[1], ":"), "metadata",
      ids[i[1]], "vs matrix", nms[i[1]], ";", "\n",
      "  an NA identifies no", if(rows) "feature" else "observation", "and",
      "cannot be compared, so the agreement of the two cannot be checked",
      config=config)
  }

  if(!all(nms == ids)) {
    i <- which(nms != ids)
    f.err(me, "the", paste0("config$", key), "column of", meta_lab,
      "does not agree with the", dim_lab, "of state$expression;", length(i),
      "of", length(nms), "differ, first at index", paste0(i[1], ":"),
      "metadata", ids[i[1]], "vs matrix", nms[i[1]], ";", "\n",
      "  they are matched by position downstream, so the metadata would be",
      "attached to the wrong", unit, config=config)
  }

  return(invisible(TRUE))
}

## needs config$feat_col and config$obs_col:

f.check_state <- function(state, config) {

  ## before either check below reads its dimnames:

  f.check_expr_matrix(state, config, "f.check_state")

  if(is.null(config$feat_col) || config$feat_col %in% "") {
    f.err("f.check_state: config$feat_col unset", config=config)
  }
  f.check_dimnames(state, config, "feat_col")

  if(is.null(config$obs_col) || config$obs_col %in% "") {
    f.err("f.check_state: config$obs_col unset", config=config)
  }
  f.check_dimnames(state, config, "obs_col")

  ## NA is only indicator of a missing value:

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

## Resolve scale of state$expression for function that takes
##   is_log_transformed argument:

f.is_log_transformed <- function(is_log_transformed, config, fn_name) {

  from_config <- config$is_log_transformed

  if(!is.null(from_config) &&
      !(is.logical(from_config) && length(from_config) %in% 1 &&
        !is.na(from_config))) {
    f.err(fn_name, ": config$is_log_transformed is not TRUE or FALSE; value:",
      from_config, config=config)
  }

  ## "" is unset for callers that pass an empty character:

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

## Resolve whether the variance prior of moderation is fitted against mean feature
##   intensity, for function that takes trend argument:

f.is_trend <- function(trend, config, fn_name="f.is_trend") {

  from_config <- config$test_trend

  if(!is.null(from_config) &&
      !(is.logical(from_config) && length(from_config) %in% 1 &&
        !is.na(from_config))) {
    f.err(fn_name, ": config$test_trend is not TRUE or FALSE; value:",
      from_config, config=config)
  }

  ## "" is accepted as unset for callers that pass empty character:

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

## Drop features that cannot contribute to fit of missingness against
##   intensity, and refuse fit outright if too little is left of it:

f.drop_unfittable <- function(dat, config, fn_name, min_fit_pts) {

  if(!(is.numeric(min_fit_pts) && length(min_fit_pts) %in% 1 &&
      is.finite(min_fit_pts) && min_fit_pts >= 2)) {
    f.err(fn_name, ": min_fit_pts is not a finite numeric scalar >= 2; value:",
      min_fit_pts, config=config)
  }

  if(!(is.data.frame(dat) && "m" %in% names(dat) && is.numeric(dat[["m"]]))) {
    f.err(fn_name, ": the missingness fit needs a data.frame with a numeric",
      "column m, holding the intensity summary of each feature;", "\n",
      " class(dat):", class(dat), "; columns:",
      if(is.null(names(dat))) "none" else names(dat), "\n",
      " class(dat$m):", if("m" %in% names(dat)) class(dat[["m"]]) else "absent",
      config=config)
  }

  if(nrow(dat) %in% 0) {
    f.err(fn_name, ": the missingness fit was given no features at all;", "\n",
      " a frame of no rows is not too little data to fit, it is a step upstream",
      "that kept nothing: check the filtering steps rather than",
      "config$impute_method", config=config)
  }

  i_drop <- is.na(dat[["m"]])

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

  n_distinct <- length(unique(dat[["m"]]))

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

## Lower bound of interval unif_ imputation methods draw from.:

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

  v <- c(state$expression)
  n_ok <- sum(!is.na(v))

  if(n_ok %in% 0) {

    f.msg("signal distribution: none to report:", length(v), "value(s), none",
      "of them measured", config=config)

  } else {

    f.msg("signal distribution:", config=config)

    f.quantile(v, config, digits=0)

    f.msg("min(state$expression):", min(v, na.rm=T),
      "; mean(state$expression):", mean(v, na.rm=T), config=config)
  }

  f.msg("num NAs: ", sum(is.na(v)), config=config)
  f.msg("num non-NAs: ", n_ok, config=config)
}

f.save_state <- function(state, config, prefix) {

  if(is.null(config$save_state) || !config$save_state) {
    f.msg("config$save_state not TRUE; no files saved", config=config)
    return(NULL)
  }

  ## three paths below built by pasting config$dir_out:

  dir_out <- config$dir_out

  if(!(length(dir_out) %in% 1 && is.character(dir_out) && !is.na(dir_out) &&
      nzchar(dir_out))) {
    f.err("f.save_state: config$save_state is TRUE, so config$dir_out has to be",
      "a single non-empty directory name;", "\n", "config$dir_out:", dir_out,
      "; class:", class(dir_out), "; length:", length(dir_out), "\n",
      "for the working directory, set it to '.'", config=config)
  }

  if(!dir.exists(dir_out)) {
    f.err("f.save_state: config$dir_out names a directory that does not exist:",
      dir_out, "\n", "nothing was written; create the directory, correct",
      "config$dir_out, or set config$save_state to FALSE", config=config)
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

  ## through f.log_obj(), so that an unwritable config$log_file prints the quantiles:
  f.log_obj(round(stats::quantile(v, probs=probs, na.rm=na.rm), digits=digits), config)
}

