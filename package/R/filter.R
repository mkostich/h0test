#' Filter features based on number of samples
#' @description
#' `f.filter_features` filters features based on number of samples expressing 
#'    feature.
#' @details Feature considered expressed if exprs > 0; NAs count as no 
#'   expression.
#' @param state A list with elements like that returned by `f.read_data()`:
#' \itemize{
#'   \item{expression}{Numeric matrix with feature rows and observation columns.}
#'   \item{features}{A data.frame with feature rows corresponding to rows of returned expression.}
#'   \item{samples}{A data.frame with observation rows corresponding to columns of returned expression.}
#' }
#' @param config List with configuration settings.
#' @return An updated `state` list with the following elements:
#' \itemize{
#'   \item{expression}{Numeric matrix with feature rows and observation columns.}
#'   \item{features}{A data.frame with feature rows corresponding to rows of returned expression.}
#'   \item{samples}{A data.frame with observation rows corresponding to columns of returned expression.}
#' }
#' @examples
#' state <- list(expression=exprs, features=feats, samples=samps)
#' config <- list(log_file="", n_samples_min=2)
#' state2 <- f.filter_features(state, config)
#' exprs2 <- state2$exprs
#' samps2 <- state2$samps
#' feats2 <- state2$feats
#'
#' config <- list(log_file="")
#' state2 <- f.filter_features(state, config, n_samples_min=2)
#' exprs2 <- state2$exprs
#' samps2 <- state2$samps
#' feats2 <- state2$feats

f.filter_features <- function(state, config, n_samples_min=NULL) {

  if(!is.matrix(state$expression)) {
    f.err("f.filter_features: !is.matrix(state$expression)", config=config)
  }
  
  if(is.null(n_samples_min)) n_samples_min <- config$n_samples_min

  f <- function(v) {
    i <- sum(v > 0, na.rm=T) >= n_samples_min
    i[is.na(i)] <- F
    return(i)
  }
  i <- apply(state$expression, 1, f)

  f.msg("filtering", sum(!i), "features, keeping", sum(i), config=config)
  state$expression <- state$expression[i, ]
  state$features <- state$features[i, ]
  f.msg("features filtered: nrow(state$expression):", nrow(state$expression), 
    "; ncol(state$features):", ncol(state$features), config=config)

  return(state)
}

#' Filter samples based on number of features
#' @description
#' `f.filter_samples` filters samples based on number of features with expression > 0
#' @details Feature considered expressed if exprs > 0; NAs count as no expression.
#' @param state A list with elements like that returned by `f.read_data()`:
#' \itemize{
#'   \item{expression}{Numeric matrix with feature rows and observation columns.}
#'   \item{features}{A data.frame with feature rows corresponding to rows of returned expression.}
#'   \item{samples}{A data.frame with observation rows corresponding to columns of returned expression.}
#' }
#' @param config List with configuration settings.
#' @return An updated `state` list with the following elements:
#' \itemize{
#'   \item{expression}{Numeric matrix with feature rows and observation columns.}
#'   \item{features}{A data.frame with feature rows corresponding to rows of returned expression.}
#'   \item{samples}{A data.frame with observation rows corresponding to columns of returned expression.}
#' }
#' @examples
#' state <- list(expression=exprs, features=feats, samples=samps)
#' config <- list(log_file="", n_features_min=1000)
#' state2 <- f.filter_features(state, config)
#' exprs2 <- state2$exprs
#' samps2 <- state2$samps
#' feats2 <- state2$feats
#'
#' config <- list(log_file="")
#' state2 <- f.filter_features(state, config, n_features_min=1000)
#' exprs2 <- state2$exprs
#' samps2 <- state2$samps
#' feats2 <- state2$feats

f.filter_samples <- function(state, config, n_features_min=NULL) {

  if(!is.matrix(state$expression)) {
    f.err("f.filter_samples: !is.matrix(state$expression)", config=config)
  }
  
  if(is.null(n_features_min)) n_features_min <- config$n_features_min

  f <- function(v) {
    i <- sum(v > 0, na.rm=T) >= n_features_min
    i[is.na(i)] <- F
    return(i)
  }
  i <- apply(state$expression, 2, f)

  f.msg("filtering", sum(!i), "samples, keeping", sum(i), config=config)
  state$expression <- state$expression[, i]
  state$samples <- state$samples[i, ]
  f.msg("samples filtered: nrow(exprs):", nrow(exprs), "; ncol(exprs):", 
    ncol(exprs), config=config)

  return(state)
}

#' Number of samples expressing each feature
#' @description
#' `f.samples_per_feature` calculates the number of samples having 
#' @details Feature considered expressed if exprs > 0; NAs count as no expression.
#' @param state A list with elements like that returned by `f.read_data()`:
#' \itemize{
#'   \item{expression}{Numeric matrix with feature rows and observation columns.}
#'   \item{features}{A data.frame with feature rows corresponding to rows of returned expression.}
#'   \item{samples}{A data.frame with observation rows corresponding to columns of returned expression.}
#' }
#' @param config List with configuration settings.
#' @return A numeric vector of length nrow(exprs) with non-negative sample 
#'   counts for each feature
#' @examples
#' state <- list(expression=exprs, features=feats, samples=samps)
#' config <- list(log_file="")
#' n_samples <- f.samples_per_feature(state, config)

f.samples_per_feature <- function(state, config) {

  if(!is.matrix(state$expression)) {
    f.err("f.samples_per_feature: !is.matrix(state$expression)", config=config)
  }

  f <- function(v) {
    n <- sum(v > 0, na.rm=T)
    n[is.na(n)] <- 0      ## only if all(is.na(v))
    return(n)
  }
  n <- apply(state$expression, 1, f)
  
  return(n)
}

#' Median expression of each feature in each expressing sample
#' @description
#' `f.feature_median_expression` calculates the median expression of each 
#'   feature in each expressing sample. 
#' @details 
#' Sample considered to express feature if exprs > 0; NAs count as no expression.
#'   Note that NAs and values <= 0 do not count toward the median.
#' @param state A list with elements like that returned by `f.read_data()`:
#' \itemize{
#'   \item{expression}{Numeric matrix with feature rows and observation columns.}
#'   \item{features}{A data.frame with feature rows corresponding to rows of returned expression.}
#'   \item{samples}{A data.frame with observation rows corresponding to columns of returned expression.}
#' }
#' @param config List with configuration settings.
#' @return A numeric vector of length nrow(exprs) with median expression in 
#'   each expressing sample
#' @examples
#' state <- list(expression=exprs, features=feats, samples=samps)
#' config <- list(log_file="")
#' median_expressions <- f.feature_median_expression(state, config)

f.feature_median_expression <- function(state, config) {

  if(!is.matrix(state$expression)) {
    f.err("f.feature_median_expression: !is.matrix(state$expression)", config=config)
  }
  
  f <- function(v) {
    m <- median(v, na.rm=T)
    m[is.na(m)] <- 0      ## only if all(is.na(v))
    return(m)
  }
  m <- apply(state$expression, 1, f)
  
  return(m)
}

#' Numer of expressed features per sample
#' @description
#' `f.features_per_sample` calculates the number of expressed features in 
#'    each sample. 
#' @details 
#' Features are considered to be expressed if exprs > 0; NAs count as no expression.
#' @param state A list with elements like that returned by `f.read_data()`:
#' \itemize{
#'   \item{expression}{Numeric matrix with feature rows and observation columns.}
#'   \item{features}{A data.frame with feature rows corresponding to rows of returned expression.}
#'   \item{samples}{A data.frame with observation rows corresponding to columns of returned expression.}
#' }
#' @param config List with configuration settings.
#' @return A numeric vector of length ncol(exprs) with number of features .
#'   expressed in each sample.
#' @examples
#' state <- list(expression=exprs, features=feats, samples=samps)
#' config <- list(log_file="")
#' n_features <- f.features_per_sample(state, config)

f.features_per_sample <- function(state, config) {

  if(!is.matrix(state$expression)) {
    f.err("f.features_per_sample: !is.matrix(state$expression)", config=config)
  }
  
  f <- function(v) {
    n <- sum(v > 0, na.rm=T)
    n[is.na(n)] <- 0      ## only if all(is.na(v))
    return(n)
  }
  n <- apply(state$expression, 2, f)
  
  return(n)
}

