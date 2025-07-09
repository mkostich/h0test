#' Combine technical replicates
#' @description
#'   Combine expression signals from technical replicates.
#' @details Combines signals for each feature across technical replicates by taking median.
#' @param state List with elements formatted like the list returned by \code{f.read_data()}:
#'   \tabular{ll}{
#'     \code{expression} \cr \tab Numeric matrix with non-negative expression values. \cr
#'     \code{features}   \cr \tab A data.frame with feature meta-data for rows of expression. \cr
#'     \code{samples}    \cr \tab A data.frame with observation meta-data for columns of expression. \cr
#'   } 
#' @param config List with configuration values. Uses the following keys:
#'   \tabular{ll}{
#'     \code{feat_col}      \cr \tab Name of column in feature metadata matching \code{rownames(state$expression)}. \cr
#'     \code{obs_col}       \cr \tab Name of column in observation metadata matching \code{colnames(state$expression)}. \cr
#'     \code{sample_id_col} \cr \tab Name of column in observation metadata with unique sample ids. \cr
#'   }
#' @return A list (the processed state) with the following elements:
#'   \tabular{ll}{
#'     \code{expression} \cr \tab Numeric matrix with non-negative expression values. \cr
#'     \code{features}   \cr \tab A data.frame with feature meta-data for rows of expression. \cr
#'     \code{samples}    \cr \tab A data.frame with observation meta-data for columns of expression. \cr
#'   } 
#' @examples
#' set.seed(101)
#' exprs <- h0testr::f.sim1(n_obs=6, n_feats=8)$mat
#' feats <- data.frame(feature_id=rownames(exprs))
#' samps <- data.frame(observation_id=colnames(exprs))
#' samps$sample_id=c("samp1", "samp2", "samp3", "samp1", "samp2", "samp3")
#' state <- list(expression=exprs, features=feats, samples=samps)
#' config <- list(feat_col="feature_id", obs_col="observation_id", sample_id_col="sample_id", save_state=FALSE)
#' out <- h0testr::f.combine_reps(state, config)
#' state
#' out$state

f.combine_reps <- function(state, config) {

  f.check_config(config)

  f <- function(v, s) tapply(v, s, stats::median, na.rm=T)
  sample_ids <- state$samples[, config$sample_id_col, drop=T]
  state$expression <- t(apply(state$expression, 1, f, sample_ids))
  state$samples <- state$samples[!duplicated(sample_ids), ]
  if(config$obs_col != config$sample_id_col) {
    state$samples[, config$obs_col] <- NULL
  }
  sample_ids <- state$samples[, config$sample_id_col, drop=T]
  if(!all(sample_ids %in% colnames(state$expression))) {
    f.err("f.normalize: !all(samples[, config$sample_id_col] %in% colnames(expression))", 
      config=config)
  }  
  state$expression <- state$expression[, sample_ids]
  config$obs_col <- config$sample_id_col

  f.check_state(state, config)
  f.report_state(state, config)
  i <- config$run_order %in% "combine_reps"
  prfx <- paste0(which(i)[1] + 2, ".combined")
  f.save_state(state, config, prefix=prfx)
  
  return(list(state=state, config=config))
}

f.combine_feats <- function(state, config) {

  f.check_config(config)
  
  return(list(state=state, config=config))
}
