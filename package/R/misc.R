f.msg <- function(..., config) {
  cat(..., "\n", file=config$log_file, append=T)
  flush.console()
}

f.log <- function(..., config) {
  cat(..., "at:", format(Sys.time(), format='%Y%m%d%H%M%S'), "\n", 
    file=config$log_file, append=T)
  flush.console()
}

f.err <- function(..., config) {
  f.log("ERROR:", ..., config=config)
  stop("Stopping", call.=F)
}

f.save_tsv <- function(dat, file_out, config, row.names=T, col.names=T) {
  tryCatch(
    write.table(dat, file=file_out, quote=F, sep="\t", 
      row.names=row.names, col.names=col.names),
      error=function(msg) f.err("write.table() error: writing to ", 
        file_out, ": ", msg$message, config=config),
      warning=function(msg) f.err("write.table() warning: writing to ", 
        file_out, ": ", msg$message, config=config)
  )
}

f.check_state <- function(state, config) {
  if(!all(rownames(exprs) == feats[, config$feat_id_col, drop=T])) {
    f.err("!all(rownames(exprs) == feats[, config$feat_id_col])", 
      config=config)
  }
  if(!all(colnames(exprs) == samps[, config$obs_col, drop=T])) {
    f.err("!all(colnames(exprs) == samps[, config$obs_col])", 
      config=config)
  }
}

f.report_state <- function(state, config) {

  f.msg("N features: ", nrow(state$exprssion), 
    "; N observations: ", ncol(state$exprssion), config=config)
    
  f.msg("signal distribution:", config=config)
  
  f.quantile(c(state$exprssion), config, digits=0) 
  
  f.msg("min(state$exprssion):", min(c(state$exprssion), na.rm=T), 
    "; mean(state$exprssion):", mean(c(state$exprssion), na.rm=T), 
    config=config)
    
  f.msg("num NAs: ", sum(is.na(c(state$exprssion))), config=config)
}

f.save_state <- function(state, config, prefix) {

  file_out <- paste0(config$dir_out, "/", prefix, 
    config$data_mid_out, config$suffix_out)
  f.log("writing expression data to", file_out, config=config)
  f.save_tsv(state$exprs, file_out, config)

  file_out <- paste0(config$dir_out, "/", prefix, 
    config$feature_mid_out, config$suffix_out)
  f.log("writing feature metadata to", file_out, config=config)
  f.save_tsv(state$feats, file_out, config)

  file_out <- paste0(config$dir_out, "/", prefix, 
    config$sample_mid_out, config$suffix_out)
  f.log("writing sample metadata to", file_out, config=config)
  f.save_tsv(state$samps, file_out, config)
}

f.quantile <- function(v, config, digits=3, na.rm=T) {
  if(config$log_file %in% "") {
    print(round(quantile(v, probs=config$probs, na.rm=na.rm), digits=digits))
  } else {
    capture.output(round(quantile(v, probs=config$probs, na.rm=na.rm), digits=digits), 
      file=config$log_file, append=T)
  }
}

