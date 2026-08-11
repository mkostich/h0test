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

