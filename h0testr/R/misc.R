f.msg <- function(..., config) {
  if(is.null(config$log_file)) config$log_file <- ""
  cat(..., "\n", file=config$log_file, append=T)
  utils::flush.console()
}

f.log <- function(..., config) {
  if(is.null(config$log_file)) config$log_file <- ""
  cat(..., "at:", format(Sys.time(), format='%Y%m%d%H%M%S'), "\n", 
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

