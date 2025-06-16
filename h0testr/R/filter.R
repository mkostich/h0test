#' Filter features based on number of samples
#' @description
#'   Filter features based on number of samples expressing feature.
#' @details Feature considered expressed if \code{state$expression > 0}; 
#'   \code{NA}s count as no expression.
#' @param state A list with elements like that returned by `f.read_data()`:
#'   \tabular{ll}{
#'     \code{expression} \cr \tab Numeric matrix with non-negative expression values. \cr
#'     \code{features}   \cr \tab A data.frame with feature meta-data for rows of expression. \cr
#'     \code{samples}    \cr \tab A data.frame with observation meta-data for columns of expression. \cr
#'   } 
#' @param config List with configuration settings.
#' @param n_samples_min Minimum number of samples expressing feature. Non-negative numeric.
#' @return An updated \code{state} list with the following elements:
#'   \tabular{ll}{
#'     \code{expression} \cr \tab Numeric matrix with non-negative expression values. \cr
#'     \code{features}   \cr \tab A data.frame with feature meta-data for rows of expression. \cr
#'     \code{samples}    \cr \tab A data.frame with observation meta-data for columns of expression. \cr
#'   } 
#' @examples
#'   \dontrun{
#'     state <- list(expression=exprs, features=feats, samples=samps)
#'     config <- list(log_file="", n_samples_min=2)
#'     state2 <- f.filter_features(state, config)
#'     exprs2 <- state2$exprs
#'     samps2 <- state2$samps
#'     feats2 <- state2$feats
#'
#'     config <- list(log_file="")
#'     state2 <- f.filter_features(state, config, n_samples_min=2)
#'     exprs2 <- state2$exprs
#'     samps2 <- state2$samps
#'     feats2 <- state2$feats
#'   } 

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
#'   Filters samples based on number of features with \code{state$expression > 0}.
#' @details 
#'   Feature considered expressed if \code{state$expression > 0}; \code{NA}s 
#'     count as no expression.
#' @param state A list with elements like that returned by \code{f.read_data()}:
#'   \tabular{ll}{
#'     \code{expression} \cr \tab Numeric matrix with non-negative expression values. \cr
#'     \code{features}   \cr \tab A data.frame with feature meta-data for rows of expression. \cr
#'     \code{samples}    \cr \tab A data.frame with observation meta-data for columns of expression. \cr
#'   } 
#' @param config List with configuration settings.
#' @param n_features_min Minimum number of features expressed per sample. Non-negative numeric.
#' @return An updated \code{state} list with the following elements:
#'   \tabular{ll}{
#'     \code{expression} \cr \tab Numeric matrix with non-negative expression values. \cr
#'     \code{features}   \cr \tab A data.frame with feature meta-data for rows of expression. \cr
#'     \code{samples}    \cr \tab A data.frame with observation meta-data for columns of expression. \cr
#'   } 
#' @examples
#'   \dontrun{
#'     state <- list(expression=exprs, features=feats, samples=samps)
#'     config <- list(log_file="", n_features_min=1000)
#'     state2 <- f.filter_features(state, config)
#'     exprs2 <- state2$exprs
#'     samps2 <- state2$samps
#'     feats2 <- state2$feats
#'
#'     config <- list(log_file="")
#'     state2 <- f.filter_features(state, config, n_features_min=1000)
#'     exprs2 <- state2$exprs
#'     samps2 <- state2$samps
#'     feats2 <- state2$feats
#'   }

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
  f.msg(
    "n_features_min:", n_features_min, 
    "; samples filtered: nrow(state$expression):", nrow(state$expression), 
    "; ncol(state$expression):", ncol(state$expression), 
    config=config
  )

  return(state)
}

#' Number of samples expressing each feature
#' @description
#'   Calculates the number of samples expressing each feature
#' @details Feature considered expressed if \code{state$expression > 0}; 
#'   \code{NA}s count as no expression.
#' @param state A list with elements like that returned by \code{f.read_data()}:
#'   \tabular{ll}{
#'     \code{expression} \cr \tab Numeric matrix with non-negative expression values. \cr
#'     \code{features}   \cr \tab A data.frame with feature meta-data for rows of expression. \cr
#'     \code{samples}    \cr \tab A data.frame with observation meta-data for columns of expression. \cr
#'   } 
#' @param config List with configuration settings.
#' @return A numeric vector of length \code{nrow(state$expression)} with non-negative sample 
#'   counts for each feature.
#' @examples
#'   \dontrun{
#'     state <- list(expression=exprs, features=feats, samples=samps)
#'     config <- list(log_file="")
#'     n_samples <- f.samples_per_feature(state, config)
#'   }

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
#'   Calculates the median expression of each feature in each expressing sample. 
#' @details 
#'   Sample considered to express feature if \code{state$expression > 0}; 
#'     \code{NA}s count as no expression.
#'     Note that \code{NA}s and values less than or equal to zero do not count 
#'       toward the median.
#' @param state A list with elements like that returned by \code{f.read_data()}:
#'   \tabular{ll}{
#'     \code{expression} \cr \tab Numeric matrix with non-negative expression values. \cr
#'     \code{features}   \cr \tab A data.frame with feature meta-data for rows of expression. \cr
#'     \code{samples}    \cr \tab A data.frame with observation meta-data for columns of expression. \cr
#'   } 
#' @param config List with configuration settings.
#' @return A numeric vector of length \code{nrow(state$expression)} with median expression in 
#'   each expressing sample.
#' @examples
#'   \dontrun{
#'     state <- list(expression=exprs, features=feats, samples=samps)
#'     config <- list(log_file="")
#'     median_expressions <- f.feature_median_expression(state, config)
#'   }

f.feature_median_expression <- function(state, config) {

  if(!is.matrix(state$expression)) {
    f.err("f.feature_median_expression: !is.matrix(state$expression)", 
      config=config)
  }
  
  f <- function(v) {
    m <- stats::median(v, na.rm=T)
    m[is.na(m)] <- 0      ## only if all(is.na(v))
    return(m)
  }
  m <- apply(state$expression, 1, f)
  
  return(m)
}

#' Numer of expressed features per sample
#' @description
#'   Calculates the number of expressed features in each sample. 
#' @details 
#'   Features are considered to be expressed if \code{state$expression > 0}; 
#'     \code{NA}s count as no expression.
#' @param state A list with elements like that returned by `f.read_data()`:
#'   \tabular{ll}{
#'     \code{expression} \cr \tab Numeric matrix with non-negative expression values. \cr
#'     \code{features}   \cr \tab A data.frame with feature meta-data for rows of expression. \cr
#'     \code{samples}    \cr \tab A data.frame with observation meta-data for columns of expression. \cr
#'   }
#' @param config List with configuration settings.
#' @return A numeric vector of length \code{ncol(state$expression)} with 
#'   number of features expressed in each sample.
#' @examples
#'   \dontrun{
#'     state <- list(expression=exprs, features=feats, samples=samps)
#'     config <- list(log_file="")
#'     n_features <- f.features_per_sample(state, config)
#'   }

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

#' Filter features and samples
#' @description
#'   Filter features and samples based on expression. 
#' @details 
#'   Filters out features with too few expressing samples, and filters out 
#'     samples with too few expressed features. Features are considered to be 
#'     expressed if \code{state$expression > 0}; \code{NA}s count as no expression.
#' @param state A list with elements like that returned by `f.read_data()`:
#'   \tabular{ll}{
#'     \code{expression} \cr \tab Numeric matrix with non-negative expression values. \cr
#'     \code{features}   \cr \tab A data.frame with feature meta-data for rows of expression. \cr
#'     \code{samples}    \cr \tab A data.frame with observation meta-data for columns of expression. \cr
#'   }
#' @param config List with configuration settings.
#' @return A list with elements like that returned by `f.read_data()`:
#'   \tabular{ll}{
#'     \code{expression} \cr \tab Numeric matrix with non-negative expression values. \cr
#'     \code{features}   \cr \tab A data.frame with feature meta-data for rows of expression. \cr
#'     \code{samples}    \cr \tab A data.frame with observation meta-data for columns of expression. \cr
#'   }
#' @examples
#'   \dontrun{
#'     state <- list(expression=exprs, features=feats, samples=samps)
#'     config <- list(log_file="")
#'     out <- f.filter(state, config)
#'     config2 <- out$config
#'     state2 <- out$state
#'     exprs2 <- state2$expression
#'     feats2 <- state2$features
#'     samps2 <- state2$samples
#'   }

f.filter <- function(state, config) {

  state <- f.filter_features(state, config)
  state <- f.filter_samples(state, config)

  n <- f.samples_per_feature(state, config)
  state$features[, config$n_samples_expr_col] <- n
  
  m <- f.feature_median_expression(state, config)
  state$features[, config$median_raw_col] <- m
  
  n <- f.features_per_sample(state, config)
  state$samples[, config$n_features_expr_col] <- n
  
  f.check_state(state, config)
  f.report_state(state, config)
  f.save_state(state, config, prefix="5.filtered")
  
  return(list(state=state, config=config))
}
