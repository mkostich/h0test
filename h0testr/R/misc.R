## Flatten arguments destined for cat() into character scalars, one per argument.
##   cat() cannot handle language objects (e.g. formulas, calls) or non-atomic
##   objects (e.g. lists), which otherwise mask the message being logged with a
##   confusing cat() error; those are deparsed instead. A zero-length argument is
##   deparsed too, rather than dropped: most such an argument is the value a
##   message exists to report, and dropping it left the message hanging on its
##   colon, e.g. 'f.parse_frm: config$frm not a formula; value:' for a config$frm
##   that is NULL. deparse() names the empty value ('NULL', 'character(0)'),
##   which is the fact the reader needs:

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

## Say once, per file, that config$log_file could not be written, then let the
##   messages themselves go to the console. Once rather than ahead of every
##   message, since a log file that cannot be opened cannot be opened for any of
##   them, and the notice repeated would bury the run's own output:

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

## Write one already-flattened message to config$log_file, falling back to the
##   console when that file cannot be opened. Every message in the package goes
##   through cat(file=config$log_file), which used to be taken on trust: a
##   log_file under a directory that is not there raised "cannot open the
##   connection" from cat() instead of printing, so a mistyped path replaced
##   every diagnostic in the package, including the error f.err() exists to
##   report, with that one message about a connection. A message that misses the
##   log is a nuisance; one that reaches nowhere loses the reason the run
##   stopped, so the console is used rather than the failure propagated:

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

f.log_obj <- function(obj, config) {

  if(is.null(config$log_file)) config$log_file <- ""

  if(config$log_file %in% "") {
    print(obj)
    utils::flush.console()
    return(invisible(NULL))
  }

  ## same fallback as f.cat_log(), for the printed objects (quantiles, tables):

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

## Write one table to file_out as a tab delimited file. The write goes to a
##   temporary file beside the destination and is moved into place only once it
##   has finished, so that a write that fails leaves whatever was already at
##   file_out as it was: write.table() streams its output, so a failure part way
##   through - a full disk, a volume that goes away - used to leave a truncated
##   file that looks like a complete one, under the name every later step reads.
##   Both error and warning stay fatal: an open that fails arrives as a warning
##   rather than an error (a read-only destination gives "cannot open file ...
##   Permission denied" that way), so demoting warnings would report a write
##   that never happened as a success. What they get instead is a message that
##   says which of the two it was and that the destination is untouched:

f.save_tsv <- function(dat, file_out, config, row.names=T, col.names=T) {

  if(!(length(file_out) %in% 1 && is.character(file_out) && !is.na(file_out) &&
      nzchar(file_out))) {
    f.err("f.save_tsv: file_out has to be a single non-empty file name;", "\n",
      " file_out:", file_out, "; class:", class(file_out), "; length:",
      length(file_out), config=config)
  }

  file_tmp <- paste0(file_out, ".tmp", Sys.getpid())

  ## the condition is returned rather than acted on inside the handler, so that
  ##   write.table()'s on.exit() has closed the connection before the temporary
  ##   file is removed; an open file cannot be removed on every platform:

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

  ## file.rename() replaces an existing destination on the platforms this runs
  ##   on, but not on all of them, and it cannot cross a file system; the copy
  ##   is the fallback for both, and the finished table is named if even that
  ##   fails, so that a completed write is never lost silently:

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

  ## what is called is normally a symbol, an operator name among f.frm_ops. It
  ##   can be a call itself, as in '(a)(b)', and as.character() of a call returns
  ##   one element per part, which the if() below cannot use: it failed with "the
  ##   condition has length > 1" rather than saying what was unsupported.
  ##   deparse() cannot stand in for as.character() here, since it backquotes an
  ##   operator name ('`+`'), which would then match nothing in f.frm_ops:

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

  ## an operator called with no operands at all, e.g. ~`+`(), which is a call of
  ##   length 1. Checked before the operands are looked at: the '-' check below
  ##   would otherwise take expr[[length(expr)]] to be the operator itself and
  ##   report it as a term that cannot be removed, and the loop at the end would
  ##   count down from 2 to 1 (R's empty range trap) and stop with 'subscript out
  ##   of bounds', naming neither the formula nor the construct:

  if(length(expr) < 2) {
    f.err("f.check_frm_expr: operator '", op, "' has no operands in formula;",
      "got", deparse(expr), "; frm rhs:", deparse(top), config=config)
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

  ## same as in f.check_frm_expr(), and more easily reached, config$contrast being
  ##   a string the caller writes: what is called can be a call itself, as in
  ##   '(grpb)(2)', and as.character() of a call returns one element per part,
  ##   which the if() below cannot use. deparse() cannot stand in for it here,
  ##   backquoting an operator name so that nothing in f.contrast_ops matches:

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

  ## an operator called with no operands at all, e.g. '`+`()', which is a call of
  ##   length 1. Checked before the operands are looked at, since the loop at the
  ##   end would otherwise count down from 2 to 1 (R's empty range trap) and stop
  ##   with 'subscript out of bounds', naming neither config$contrast nor the
  ##   construct. '*' and '/' have their own two-operand check below, which this
  ##   one precedes only so that every operator is covered by one message:

  if(length(expr) < 2) {
    f.err("f.check_contrast_expr: operator '", op, "' has no operands in",
      "config$contrast; got", deparse(expr), "; config$contrast:", deparse(top),
      config=config)
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

  ## the same bindings again, in an environment where '+' and '-' add magnitudes
  ##   instead of cancelling, so that evaluating the expression there gives the
  ##   scale the weights would have had had nothing cancelled; see below, where
  ##   that is what distinguishes a contrast of zero from a small one. '(' is left
  ##   to baseenv() in both, and the two are built together so that they cannot
  ##   come to hold different coefficients:

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

  ## whether there is a hypothesis at all, which the size of the weights alone
  ##   cannot answer. The hypothesis is that the weighted sum is zero, and that is
  ##   invariant under rescaling the weights, so '(grpb - grpc)/1e9' asks exactly
  ##   what 'grpb - grpc' asks while carrying weights of 1e-9. A contrast written
  ##   so that it cancels, on the other hand, such as
  ##   '(grpb + grpc)/3 - grpb/3 - grpc/3', is left with weights of about 1e-17
  ##   whose direction is rounding noise rather than a hypothesis. An absolute
  ##   threshold cannot separate those two, and the one used here,
  ##   sqrt(.Machine$double.eps) or about 1.5e-8, refused the first as readily as
  ##   the second. What separates them is whether cancellation happened, which is
  ##   measured rather than guessed at: env_abs above evaluates the same
  ##   expression with the magnitudes added instead of cancelled, so its result is
  ##   the scale the weights would have had without any cancellation, and a weight
  ##   negligible against that scale is what is left of a cancellation rather than
  ##   a small contrast. The weights are not rescaled to suit the comparison: the
  ##   hypothesis would be unchanged, but the effect size reported from them would
  ##   not, '2*grpb - 2*grpc' being twice the difference 'grpb - grpc' reports.
  ##   The two causes are reported separately, since the remedy differs: a
  ##   contrast that is zero by construction has to be rewritten, while one that
  ##   cancels was probably meant to be written some other way.
  ##   A no_cancel that does not come back usable leaves only the exact zero
  ##   refused, which is the weakest safe answer; it cannot happen for an
  ##   expression f.check_contrast_expr() has passed, every operator in it being
  ##   one of the four bound above or '(':

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

## A design with fewer independent columns than columns cannot support a test of
##   one of its terms: some coefficient is a linear combination of the others, so
##   the data do not decide which of them carries the effect. Refused rather than
##   left to the engines, which used to be handed such a design and each failed
##   their own way: limma stopped with "Coefficients not estimable: grpd" and then
##   "Subsetting to non-estimable coefficients is not allowed", naming neither the
##   covariate nor config, while the engines that fit feature by feature reported
##   NAs. Two routes reach one: a factor covariate carrying a level that nothing
##   was measured at, whose column model.matrix() codes as all zeros, and two
##   terms of config$frm that code the same grouping.
##   Only for the config$test_term branch. A config$contrast is a weighted sum of
##   coefficients, which can be estimable in a rank deficient design (the sum of
##   two aliased coefficients is, either one alone is not), and
##   f.design_contrast() tests exactly that by asking whether the contrast lies in
##   the row space, so the deficiency is not by itself a reason to refuse there.
##   Per feature estimability, where the design is full rank over all observations
##   but not over the observations one feature was measured in, is a third matter,
##   screened by filter_features_by_estimability() from this same design:

f.design_check_rank <- function(X, rank_all, parsed, config,
    caller="f.design_test_cols") {

  if(!(rank_all > 0 && rank_all < ncol(X))) return(invisible(NULL))

  ## qr() pivots the columns it could not use to the end, which names them:

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
        "h0testr::test() and the other fitting entry points drop them for you;",
        "a design built directly from state$samples does not")
    } else {
      paste("no column is all zeros, so two terms of config$frm code the same",
        "grouping: drop the redundant one")
    }, config=config)
}

## Rank of a design matrix, tolerant of the degenerate shapes that arise when a
##   feature was measured in too few observations, or when the reduced model has
##   no columns left:

f.design_rank <- function(mat) {
  if(nrow(mat) %in% 0 || ncol(mat) %in% 0) return(0L)
  return(qr(mat)$rank)
}

## Helper for f.normalize_terms(), which is a helper for f.design_test_cols(),
##   which serves test() and filter_features_by_estimability(). Converts formula
##   config$frm to character, makes intercept explicit (either '0' or '1'), sorts
##   variables in interaction terms (so e.g. 'sex:age' becomes 'age:sex'), then
##   returns formula representation as tokenized character vector. So would take
##   formula e.g. ~age + strain + strain:age, and return:
##   c("1", "age", "strain", "age:strain").

f.formula2terms <- function(config) {

  frm <- config$frm

  if(is.null(frm) || all(as.character(frm) %in% "")) {
    f.err("f.formula2terms: config$frm empty or undefined.", config=config)
  }

  parsed <- f.parse_frm(frm, config)

  ## make intercept explicit, so that it can be dropped or kept downstream;
  ##   '*', '^', and '/' expansion, as well as intercept removal, are handled
  ##   by stats::terms() within f.parse_frm(), so e.g. ~x1*x2 and
  ##   ~x1+x2+x1:x2 both give c("1", "x1", "x2", "x1:x2"):

  if(parsed$intercept %in% 1) {
    terms <- c("1", parsed$labels)
  } else {
    terms <- c("0", parsed$labels)
  }

  return(terms)
}

## Helper for f.normalize_terms(). Given a parsed formula (see f.parse_frm())
##   and a canonical config$test_term, returns the character vector of term
##   labels to be dropped from the full model to form the reduced model.
##   A test_term naming a variable drops that variable's term along with every
##   higher-order term containing it (so testing 'x1' in ~x1*x2 is a 2 df test
##   of 'x1' and 'x1:x2'); this keeps the reduced model hierarchical, so the
##   test does not depend on the contrast coding or on which level of x2 is the
##   reference. A test_term naming an interaction term drops just that term,
##   and is an error if some higher-order term in the formula contains it:

f.test_term_drops <- function(parsed, test_term, config) {

  if(test_term %in% "1") {
    if(parsed$intercept %in% 0) {
      f.err("f.test_term_drops: test_term is intercept ('1'), but frm has no",
        "intercept; frm:", parsed$frm, config=config)
    }
    return("1")
  }

  vars <- rownames(parsed$factors)

  ## test_term names a variable: drop every term it takes part in:

  if(test_term %in% vars) {
    drops <- parsed$labels[parsed$factors[test_term, ] > 0]
    return(drops)
  }

  ## test_term names a term of the model: drop it alone, but only if no
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

## Helper for f.design_test_cols(). Character scalar config$test_term, formula
##   config$frm; returns list with canonicalized character scalar $test_term,
##   tokenized character vector $frm_terms carrying either '1' for intercept or
##   '0' for no intercept, and character vector $drop_terms with the terms to be
##   dropped from the full model to form the reduced model. Interaction terms in
##   $frm_terms and $test_term are sorted alphabetically (so 'sex:age' becomes
##   'age:sex'), to facilitate formula/term comparison.

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
  ##   an error if test_term is not compatible with frm:

  drop_terms <- f.test_term_drops(parsed, test_term, config)

  return(list(test_term=test_term, frm_terms=frm_terms, drop_terms=drop_terms))
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

  ## before the df_intend arithmetic below, which cannot tell a design that is
  ##   short of rank from a term that is redundant, and used to hand both on:

  f.design_check_rank(X, rank_all, parsed, config, caller="f.design_test_cols")

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

  ## the advice depends on the cause, and the wrong advice sends the reader to the
  ##   intercept when config$test_term never named it, or to config$frm when the
  ##   formula is fine and the observations are gone. Three causes: a full model
  ##   with no rank, where there is nothing to fit either way; the intercept named
  ##   alongside a factor the reduced model codes to full rank; and a term that
  ##   adds no rank to what the remaining terms already span:

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
          "drop the redundant term, or test one that is not implied. A design",
          "whose columns are outright aliased is refused before this, by",
          "f.design_check_rank(), so the redundancy here is between terms rather",
          "than between columns")
      }
    f.err("f.design_test_cols: dropping config$test_term '", config$test_term,
      "' leaves a reduced model spanning the same space as the full model,",
      "so there is no hypothesis to test;", "\n",
      "config$frm:", deparse(parsed$frm), "; terms dropped to form the",
      "reduced model:", paste(drops, collapse=" "), "; design columns:",
      ncol(X), "; rank:", rank_all, "\n", hint, config=config)
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
  ## which coefficients the contrast weights, for the messages below and for the
  ##   effect size f.logfc_effect() forms from L[cols_test]. Relative to the
  ##   largest weight rather than absolute: the weights carry whatever scale the
  ##   caller wrote them in, so an absolute cut called every weight of
  ##   '(grpb - grpc)/1e9' unweighted, while a weight this far below the largest
  ##   moves the estimate by no more than rounding does. f.contrast_vector() has
  ##   already refused an L that is zero, so the largest weight is positive:

  cols_test <- which(abs(L) >= max(abs(L)) * sqrt(.Machine$double.eps))

  ## the contrast has to be estimable over all observations, which it is exactly
  ##   when it lies in the row space of the design: otherwise no linear combination
  ##   of the fitted values estimates it, and every engine would report either a
  ##   silently aliased coefficient or an error of its own naming neither
  ##   config$contrast nor config$frm. Reached by weighting a coefficient that
  ##   config$frm codes as aliased with another. Per feature estimability is a
  ##   separate matter, screened by filter_features_by_estimability() from the same
  ##   pair of designs:

  ## normalized by the largest weight, rather than by the 2-norm, before the rank
  ##   test and the null space basis below: the weights carry whatever scale the
  ##   caller wrote them in, and squaring them loses that scale at both ends of the
  ##   double range. sum(L^2) of '(grpb - grpc)/1e200' is exactly 0, so L/sqrt(...)
  ##   is NaN and Inf, and the rank test dies inside LAPACK ("NA/NaN/Inf in foreign
  ##   function call (arg 1)"), naming neither config$contrast nor config$frm;
  ##   sum(L^2) of the same contrast times 1e200 is Inf, so the normalized row is
  ##   all zeros, which the row space of every design contains, and an unestimable
  ##   contrast then passes the check below vacuously. Dividing by max(abs(L))
  ##   overflows and underflows for no L that f.contrast_vector() returns, and is
  ##   positive because f.contrast_vector() has already refused an L that is zero.
  ##   Both the rank test and the null space are scale invariant, so which of the
  ##   two normalizers is used changes nothing else. L itself is left as written,
  ##   since f.logfc_effect() forms the reported effect size from it:

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
##   config$frm. config$covariate_types records what this function decided, and
##   is checked against the columns rather than returned in place of them; see
##   below:

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

  ## the cache as a cross-check, not as an override. It used to be returned
  ##   whenever its names covered config$frm and its values were spelled right,
  ##   which let it name a type the column cannot have: a character or factor
  ##   covariate was then fit as continuous, or a numeric one split into levels,
  ##   and rule 5 above stopped refusing an undeclared character covariate, so
  ##   its reference level went back to being set by sorting. The type of a
  ##   column does not depend on which observations are left, so what the cache
  ##   says has to agree with the rules above; where it does not, the config was
  ##   built against other data or edited by hand, and neither answer can be
  ##   assumed to be the intended one. Deriving the types is cheap - one class()
  ##   per variable - so the cache saves nothing worth this. Variables the cache
  ##   covers and config$frm does not are ignored, config$frm being free to
  ##   change between calls:

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
        "config$covariate_types is set by initialize() for the state it was",
        "run on, so a disagreement means this config is being used with other",
        "data, or was edited by hand;", "\n",
        "set config$covariate_types to NULL to have it derived from",
        "state$samples, or, if the covariate is meant to be categorical, name",
        "it in config$reference_levels instead", config=config)
    }
  }

  return(out)
}

## The level ordering of one factor covariate, resolved from config. Single place
##   where that ordering is decided, so that it does not depend on which entry
##   point asked for it:
##     config$factor_levels[[trm]]      what initialize() resolved and logged;
##                                        preferred, since it is the ordering the
##                                        rest of the run has been reported under
##     config$reference_levels[[trm]]    the declared reference level first and
##                                        the remaining values sorted, which is
##                                        how initialize() resolves it; used when
##                                        initialize() has not run
##     levels(factor(v))                 neither declared: a covariate that is
##                                        already a factor, or is logical,
##                                        carries an ordering of its own that is
##                                        not locale dependent
##   config$reference_levels only ever names the first level, so the rest are
##   sorted; their order does not set the reference level, but it does set the
##   order in which the coefficients are reported. A declared reference level that
##   is not among the values is an error rather than an ordering, since it is a
##   declaration about data that are not there:

f.covariate_levels <- function(v, trm, config, caller="f.covariate_levels") {

  lvls <- config$factor_levels[[trm]]
  if(!is.null(lvls)) return(lvls)

  ## NOTE: config$reference_levels is an atomic vector, so [[ throws on a name
  ##   that is not present, rather than returning NULL:

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

## Rebuild a factor covariate with the level ordering f.covariate_levels()
##   resolves, with the declared reference level first. Handing a downstream fit a
##   character vector instead leaves it to re-derive the levels by sorting, which
##   renames the coefficients: with reference level "M", h0testr names the tested
##   column sexF while an alphabetical re-derivation names it sexM. Where the fit
##   reports one named contrast (test_msqrob()) that mismatch yields a table of
##   NAs with no error; where it reports a coefficient (test_prolfqua()) it changes
##   which contrast the reported numbers describe. Both are silent, hence this:

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

## Rebuild the covariates of config$frm in state$samples whose level ordering
##   config declares, so that the reference level a run reports is the declared one
##   whether or not initialize() has run. Without this, config$reference_levels
##   reaches the design only through initialize(), which is what turns it into
##   factor columns and into config$factor_levels: f.design_X() hands
##   state$samples to stats::model.matrix() as they are, and a character column
##   there is levelled by sorting, so a direct caller got coefficients named for a
##   reference level it did not ask for, with nothing said about it. Called once
##   per entry point that fits anything, rather than inside f.design_X(), because
##   the engines that build a model frame of their own (test_proda(),
##   test_msqrob(), test_prolfqua()) take it from state$samples rather than from
##   the design. The level ordering of a covariate config declares nothing about is
##   left exactly as it is, so what stats::model.matrix() does with it is unchanged.
##   Levels no observation has are dropped either way, since model.matrix() codes
##   such a level as a column of zeros and f.design_X() refuses the rank deficient
##   design that results. Two routes reach one: a factor column the caller built
##   with a level that was never measured, and config$factor_levels, which
##   initialize() records over all observations and which a later filtering step can
##   empty out. Dropped rather than refused, since the levels that remain are a
##   model that can be fit and are what the data support; said out loud, because
##   dropping a level renames no coefficient but removes one, so a config$contrast
##   naming it stops being parseable:

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

  ## every check below is over the values of a covariate, and a table of no rows
  ##   has none: is.na() of nothing is logical(0), unique() of nothing has length
  ##   0, so the missing, non-finite, blank, and constant checks all passed and
  ##   the covariates were reported as checked. The run then stopped at the first
  ##   design matrix, in stats::model.matrix(), with "contrasts can be applied
  ##   only to factors with 2 or more levels", which names neither the empty
  ##   table nor the step that emptied it:

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

    ## a covariate that is not a column of state$samples at all: every check
    ##   below would otherwise be run on the NULL that state$samples[[nom]]
    ##   returns, where is.na(NULL) is logical(0) and unique(NULL) has length 0,
    ##   so all of them pass vacuously and the covariate is reported as checked.
    ##   Reachable only from a caller that supplies types itself, since
    ##   f.covariate_types() checks membership, but silence is the wrong answer
    ##   for a check whose whole job is to refuse a covariate that cannot be fit:

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

## helper for f.check_state() and f.check_parameters(): state$expression has to
##   be a numeric matrix. Only is.matrix() was checked, which a character matrix
##   passes, and nothing downstream looks again: the exact zero check in
##   f.check_state() compares with %in% 0, which no string can match, and
##   f.zeros_to_na() compares with a numeric as well, so a matrix of text got as
##   far as the first arithmetic and failed there, naming neither the matrix nor
##   the step that should have refused it. read_data() checks the type of what it
##   reads off disk, but a state assembled in the session, which is how every
##   example builds one, arrives here unchecked. An integer matrix is numeric and
##   passes; a logical one does not, no normalization or imputation of TRUE and
##   FALSE being defined. fn_name only labels the messages:

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

## helper for f.check_state() and f.check_parameters(): the metadata rows and
##   state$expression are matched by position everywhere downstream, a column of
##   the matrix being the observation described by the corresponding row of
##   state$samples, so the metadata has to have one row per feature or
##   observation, and the ids it carries have to agree with the dimnames the
##   matrix carries. `==` returns logical(0) when either side is NULL, and
##   all(logical(0)) is TRUE, so comparing the two directly passed vacuously in
##   exactly the cases where the agreement could not be checked: a
##   state$expression with no dimnames, and a key naming no column of the
##   metadata. Unequal lengths were not caught either, `==` recycling the shorter
##   side, so metadata of twice the length compared against itself and passed.
##   All refused rather than assumed. Dimensions are checked before dimnames,
##   being answerable whether or not the matrix is named. key is "feat_col" or
##   "feat_id_col" for the rows, "obs_col" or "obs_id_col" for the columns;
##   f.check_state() has the first pair, f.check_parameters() runs before they
##   are set and has the second. fn_name only labels the messages:

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

  ## meta[[nom]] with a nom that is not a single name raises a raw R error from
  ##   the function whose whole purpose is to explain the mismatch: "attempt to
  ##   select less than one element in get1index" for a NULL key, and "subscript
  ##   out of bounds" for a key of two names, [[ reading those as recursive
  ##   indexing rather than as two columns. check_config() refuses a key of the
  ##   wrong shape, but f.check_parameters() checks only that the key names a
  ##   column, so a two-name key reached here whenever both happened to be
  ##   columns. The empty string is a name like any other here: it names no
  ##   column, and the message below says so:

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
      "h0testr::initialize() sets config$feat_col and config$obs_col from",
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

  ## an NA on either side compares to NA rather than to FALSE, so all() returned
  ##   NA and the if() below raised "missing value where TRUE/FALSE needed",
  ##   losing the message it was there to print. Reported here rather than made
  ##   NA-safe below, since a missing id is a different problem from a
  ##   disagreement: there is no id to compare, so no ordering can be right, and
  ##   the remedy is to supply one rather than to reorder. utils::read.table(),
  ##   which read_data() uses, reads the string NA as a missing value, so an id
  ##   of that name arrives as a gap; f.check_parameters() does not catch it
  ##   either, duplicated() not counting a single NA as a duplicate:

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

  ## dat$m used to be read on trust, and `$` partial matches on a data.frame:
  ##   a frame whose only m-like column was something else entirely (mean_int,
  ##   say) had that column fit as the intensity summary and was returned
  ##   unchanged, and a frame with no such column at all reached the degeneracy
  ##   check below as zero distinct intensities, which reported a caller's
  ##   mistake as too little data and advised changing config$impute_method.
  ##   Read by name from here on, and the shape stated before it is used:

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

  v <- c(state$expression)
  n_ok <- sum(!is.na(v))

  ## min() of nothing is Inf and mean() of nothing is NaN, both after a warning
  ##   that this function's own message then hid: an all-NA matrix, and a matrix
  ##   with no features at all, were reported as a signal distribution with a
  ##   minimum of Inf rather than as the absence of one. Both are states a
  ##   filtering step can produce, and this is where they are meant to be seen:

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

  ## the three paths below are built by pasting config$dir_out, which used to be
  ##   taken as given. A NULL or empty one made '/3.normalized.expression.tsv',
  ##   the root of the filesystem rather than a run directory, and one naming a
  ##   directory that is not there was reported only by f.save_tsv()'s warning
  ##   handler, after the log had said the data were being written. check_config()
  ##   catches a dir_out of the wrong shape, but only for a caller that runs it,
  ##   and it cannot know whether the directory exists at the time of the write;
  ##   checked here because this is the only function that writes state to disk.
  ##   Refused rather than created: an output directory that is not there usually
  ##   means the path is wrong, and creating it scatters run output over the
  ##   filesystem instead of saying so:

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

  ## through f.log_obj(), so that an unwritable config$log_file prints the
  ##   quantiles rather than stopping the run on the connection:

  f.log_obj(round(stats::quantile(v, probs=probs, na.rm=na.rm), digits=digits),
    config)
}

