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

## TRUE iff x is a single one-sided formula, e.g. ~age+gender+age:gender.
##   NOTE: length() of a formula is the number of components of the underlying
##   call: 2 for a one-sided formula (`~` and rhs), 3 for a two-sided formula
##   (`~`, lhs, and rhs); it is never 1. Two-sided formulas are rejected here
##   because the rest of the package assumes the dependent variable is implicit
##   (e.g. as.character(frm)[2] is taken to be the full set of terms):

f.is_formula <- function(x) {
  return(inherits(x, "formula") && length(x) %in% 2)
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

