f.msg <- function(..., config) {
  cat(..., "\n", file=config$log_file, append=T)
  utils::flush.console()
}

f.log <- function(..., config) {
  cat(..., "at:", format(Sys.time(), format='%Y%m%d%H%M%S'), "\n", 
    file=config$log_file, append=T)
  utils::flush.console()
}

f.log_block <- function(..., config) {
  f.msg("", config=config)
  f.log(..., config=config)
}

f.err <- function(..., config) {
  f.log("ERROR:", ..., config=config)
  stop("Stopping", call.=F)
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

f.check_state <- function(state, config) {

  feats <- state$features[, config$feat_id_col, drop=T]
  if(!all(rownames(state$expression) == feats)) {
    f.err("state$features do not match rows of state$expression", 
      config=config)
  }
  
  if(config$obs_col %in% "") config$obs_col <- config$obs_id_col
  samps <- state$samples[, config$obs_col, drop=T]
  if(!all(colnames(state$expression) == samps)) {
    f.err("state$samples do not match columns of state$expression", 
      config=config)
  }
}

f.report_state <- function(state, config) {

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

  file_out <- paste0(config$dir_out, "/", prefix, 
    config$data_mid_out, config$suffix_out)
  f.log("writing expression data to", file_out, config=config)
  f.save_tsv(state$expression, file_out, config)

  file_out <- paste0(config$dir_out, "/", prefix, 
    config$feature_mid_out, config$suffix_out)
  f.log("writing feature metadata to", file_out, config=config)
  f.save_tsv(state$feats, file_out, config)

  file_out <- paste0(config$dir_out, "/", prefix, 
    config$sample_mid_out, config$suffix_out)
  f.log("writing sample metadata to", file_out, config=config)
  f.save_tsv(state$samps, file_out, config)
}

f.quantile <- function(v, config, probs=NULL, digits=3, na.rm=T) {

  if(is.null(probs)) probs <- config$probs
  
  if(config$log_file %in% "") {
    print(round(stats::quantile(v, probs=probs, na.rm=na.rm), digits=digits))
  } else {
    utils::capture.output(round(stats::quantile(v, probs=probs, na.rm=na.rm), digits=digits), 
      file=config$log_file, append=T)
  }
}

