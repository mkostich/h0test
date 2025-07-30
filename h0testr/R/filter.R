#' Filter features based on number of samples
#' @description
#'   Filter features based on number of samples expressing feature.
#' @details
#'   Feature considered expressed if \code{state$expression > 0}; 
#'     \code{NA}s count as no expression.
#'   Feature constant if \code{length(unique(expression_values)) \%in\% 1}.
#'   See documentation for \code{h0testr::new_config()} 
#'     for more detailed description of configuration parameters. 
#' @param state A list with elements like that returned by \code{read_data()}:
#'   \tabular{ll}{
#'     \code{expression} \cr \tab Numeric matrix with non-negative expression values. \cr
#'     \code{features}   \cr \tab A data.frame with feature meta-data for rows of expression. \cr
#'     \code{samples}    \cr \tab A data.frame with observation meta-data for columns of expression. \cr
#'   } 
#' @param config List with configuration values. Uses the following keys:
#'   \tabular{ll}{
#'     \code{n_samples_min} \cr \tab Minimum number (non-negative numeric) of samples expressing feature to keep feature. \cr
#'   }
#' @param n_samples_min Minimum number of samples expressing feature. Non-negative numeric.
#' @param remove_constant Logical scalar: if constant features of \code{state$expression} should be removed.
#' @return An updated \code{state} list with the following elements:
#'   \tabular{ll}{
#'     \code{expression} \cr \tab Numeric matrix with non-negative expression values. \cr
#'     \code{features}   \cr \tab A data.frame with feature meta-data for rows of expression. \cr
#'     \code{samples}    \cr \tab A data.frame with observation meta-data for columns of expression. \cr
#'   } 
#' @examples
#' set.seed(101)
#' exprs <- h0testr::sim1(n_obs=6, n_feats=12, mcar_p=0.5)$mat
#' feats <- data.frame(feature_id=rownames(exprs))
#' samps <- data.frame(observation_id=colnames(exprs))
#' state <- list(expression=exprs, features=feats, samples=samps)
#' config <- list(n_samples_min=3)
#' state2 <- h0testr::filter_features(state, config)
#' print(state)
#' print(state2)

filter_features <- function(state, config, 
    n_samples_min=NULL, remove_constant=TRUE) {

  check_config(config)
  
  if(!is.matrix(state$expression)) {
    f.err("filter_features: !is.matrix(state$expression)", config=config)
  }
  
  if(is.null(n_samples_min)) n_samples_min <- config$n_samples_min
  if(is.null(n_samples_min)) {
    f.err("filter_features: n_samples_min unset", config=config)
  }
  
  f <- function(v) {
    i <- sum(v > 0, na.rm=T) >= n_samples_min
    i[is.na(i)] <- F
    return(i)
  }
  i <- apply(state$expression, 1, f)
  
  f.msg("filtering", sum(!i), "features, keeping", sum(i), config=config)
  state$expression <- state$expression[i, , drop=F]
  state$features <- state$features[i, , drop=F]
  f.msg("features filtered: nrow(state$expression):", nrow(state$expression), 
    "; ncol(state$features):", ncol(state$features), config=config)

  if(remove_constant) {
    f <- function(v) {
      i_na <- is.na(v)
      if(any(!i_na)) {
        return(length(unique(v[!i_na])) >= 2)
      } else {
        return(FALSE)
      }
    }
    i <- apply(state$expression, 1, f)
    i[is.na(i)] <- F
    f.msg("filtering", sum(!i), "constant features; keeping", sum(i), config=config)
    state$expression <- state$expression[i, , drop=F]
    state$features <- state$features[i, , drop=F]
    f.msg("constant features removed:", sum(!i), config=config)
    f.msg("non-constant features left:", nrow(state$features), config=config)
  }
  
  return(state)
}

#' Filter samples based on number of features
#' @description
#'   Filter samples based on number of features with \code{state$expression > 0}.
#' @details 
#'   Feature considered expressed if \code{state$expression > 0}; \code{NA}s 
#'     count as no expression.
#'   Sample constant if \code{length(unique(expression_values)) \%in\% 1)}.
#'   See documentation for \code{h0testr::new_config()} 
#'     for more detailed description of configuration parameters. 
#' @param state A list with elements like that returned by \code{read_data()}:
#'   \tabular{ll}{
#'     \code{expression} \cr \tab Numeric matrix with non-negative expression values. \cr
#'     \code{features}   \cr \tab A data.frame with feature meta-data for rows of expression. \cr
#'     \code{samples}    \cr \tab A data.frame with observation meta-data for columns of expression. \cr
#'   } 
#' @param config List with configuration values. Uses the following keys:
#'   \tabular{ll}{
#'     \code{n_features_min} \cr \tab Minimum number (non-negative numeric) of features expressed in observation to keep observation. \cr
#'   }
#' @param n_features_min Minimum number of features expressed per sample. Non-negative numeric.
#' @param remove_constant Logical scalar: if constant observations of \code{state$expression} should be removed.
#' @return An updated \code{state} list with the following elements:
#'   \tabular{ll}{
#'     \code{expression} \cr \tab Numeric matrix with non-negative expression values. \cr
#'     \code{features}   \cr \tab A data.frame with feature meta-data for rows of expression. \cr
#'     \code{samples}    \cr \tab A data.frame with observation meta-data for columns of expression. \cr
#'   } 
#' @examples
#' set.seed(101)
#' exprs <- h0testr::sim1(n_obs=6, n_feats=12, mcar_p=0.5)$mat
#' feats <- data.frame(feature_id=rownames(exprs))
#' samps <- data.frame(observation_id=colnames(exprs))
#' state <- list(expression=exprs, features=feats, samples=samps)
#' config <- list(n_features_min=4)
#' state2 <- h0testr::filter_observations(state, config)
#' print(state)
#' print(state2)

filter_observations <- function(state, config, 
    n_features_min=NULL, remove_constant=TRUE) {

  check_config(config)

  if(!is.matrix(state$expression)) {
    f.err("filter_observations: !is.matrix(state$expression)", config=config)
  }
  
  if(is.null(n_features_min)) n_features_min <- config$n_features_min
  if(is.null(n_features_min)) {
    f.err("filter_observations: n_features_min unset", config=config)
  }
  
  f <- function(v) {
    i <- sum(v > 0, na.rm=T) >= n_features_min
    i[is.na(i)] <- F
    return(i)
  }
  i <- apply(state$expression, 2, f)
  f.msg("filtering", sum(!i), "observations, keeping", sum(i), config=config)
  state$expression <- state$expression[, i, drop=F]
  state$samples <- state$samples[i, , drop=F]
  f.msg("n_features_min:", n_features_min, config=config)
  f.msg("observations filtered:", sum(!i), config=config)
  f.msg("observations left:", ncol(state$expression), config=config)
  
  if(remove_constant) {
    f <- function(v) {
      i_na <- is.na(v)
      if(any(!i_na)) {
        return(length(unique(v[!i_na])) >= 2)
      } else {
        return(FALSE)
      }
    }
    i <- apply(state$expression, 2, f)
    f.msg("filtering", sum(!i), "constant observations, keeping", sum(i), config=config)
    state$expression <- state$expression[, i, drop=F]
    state$samples <- state$samples[i, , drop=F]
    f.msg("constant observations removed:", sum(!i), config=config)
    f.msg("non-constant observations left:", nrow(state$expression), config=config)
  }

  return(state)
}

#' Number of samples expressing each feature
#' @description
#'   Calculates the number of samples expressing each feature
#' @details Feature considered expressed if \code{state$expression > 0}; 
#'   \code{NA}s count as no expression.
#' @param state A list with elements like that returned by \code{read_data()}:
#'   \tabular{ll}{
#'     \code{expression} \cr \tab Numeric matrix with non-negative expression values. \cr
#'     \code{features}   \cr \tab A data.frame with feature meta-data for rows of expression. \cr
#'     \code{samples}    \cr \tab A data.frame with observation meta-data for columns of expression. \cr
#'   } 
#' @param config List with configuration values. Does not use any params, so can pass empty list.
#' @return A numeric vector of length \code{nrow(state$expression)} with non-negative sample 
#'   counts for each feature.
#' @examples
#' set.seed(101)
#' exprs <- h0testr::sim1(n_obs=6, n_feats=12, mcar_p=0.5)$mat
#' feats <- data.frame(feature_id=rownames(exprs))
#' samps <- data.frame(observation_id=colnames(exprs))
#' state <- list(expression=exprs, features=feats, samples=samps)
#' config <- list()
#' h0testr::samples_per_feature(state, config)

samples_per_feature <- function(state, config) {

  if(!is.matrix(state$expression)) {
    f.err("samples_per_feature: !is.matrix(state$expression)", config=config)
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
#' @param state A list with elements like that returned by \code{read_data()}:
#'   \tabular{ll}{
#'     \code{expression} \cr \tab Numeric matrix with non-negative expression values. \cr
#'     \code{features}   \cr \tab A data.frame with feature meta-data for rows of expression. \cr
#'     \code{samples}    \cr \tab A data.frame with observation meta-data for columns of expression. \cr
#'   } 
#' @param config List with configuration values. Does not use any params, so can pass empty list.
#' @return A numeric vector of length \code{nrow(state$expression)} with median expression in 
#'   each expressing sample.
#' @examples
#' set.seed(101)
#' exprs <- h0testr::sim1(n_obs=6, n_feats=12)$mat
#' feats <- data.frame(feature_id=rownames(exprs))
#' samps <- data.frame(observation_id=colnames(exprs))
#' state <- list(expression=exprs, features=feats, samples=samps)
#' config <- list()
#' h0testr::feature_median_expression(state, config)

feature_median_expression <- function(state, config) {

  if(!is.matrix(state$expression)) {
    f.err("feature_median_expression: !is.matrix(state$expression)", 
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

#' Number of expressed features per sample
#' @description
#'   Calculates the number of expressed features in each sample. 
#' @details 
#'   Features are considered to be expressed if \code{state$expression > 0}; 
#'     \code{NA}s count as no expression.
#' @param state A list with elements like that returned by `read_data()`:
#'   \tabular{ll}{
#'     \code{expression} \cr \tab Numeric matrix with non-negative expression values. \cr
#'     \code{features}   \cr \tab A data.frame with feature meta-data for rows of expression. \cr
#'     \code{samples}    \cr \tab A data.frame with observation meta-data for columns of expression. \cr
#'   }
#' @param config List with configuration values. Does not use any params, 
#'   so can pass empty list.
#' @return A numeric vector of length \code{ncol(state$expression)} with 
#'   number of features expressed in each sample.
#' @examples
#' set.seed(101)
#' exprs <- h0testr::sim1(n_obs=6, n_feats=12)$mat
#' feats <- data.frame(feature_id=rownames(exprs))
#' samps <- data.frame(observation_id=colnames(exprs))
#' state <- list(expression=exprs, features=feats, samples=samps)
#' config <- list()
#' h0testr::features_per_sample(state, config)

features_per_sample <- function(state, config) {

  if(!is.matrix(state$expression)) {
    f.err("features_per_sample: !is.matrix(state$expression)", config=config)
  }
  
  f <- function(v) {
    n <- sum(v > 0, na.rm=T)
    n[is.na(n)] <- 0      ## only if all(is.na(v))
    return(n)
  }
  n <- apply(state$expression, 2, f)
  
  return(n)
}

#' Prefilter data
#' @description
#'   Adds filtering-related statistics to \code{state$features}, 
#'     and \code{state$samples}.
#' @details 
#'   Wrapper for \code{samples_per_feature()}, \code{feature_median_expression()}, 
#'     \code{features_per_sample()}. Also reports quantiles of distributions. 
#'   See documentation for \code{h0testr::new_config()} 
#'     for more detailed description of configuration parameters. 
#' @param state List with elements formatted like the list returned by \code{read_data()}:
#'   \tabular{ll}{
#'     \code{expression} \cr \tab Numeric matrix with non-negative expression values. \cr
#'     \code{features}   \cr \tab A data.frame with feature meta-data for rows of expression. \cr
#'     \code{samples}    \cr \tab A data.frame with observation meta-data for columns of expression. \cr
#'   } 
#' @param config List with configuration values. Requires the following keys:
#'   \tabular{ll}{
#'     \code{feat_id_col}  \cr \tab Name of column (character) in \code{feature_file_in} that corresponds to rows of \code{data_file_in}. \cr
#'     \code{obs_id_col}   \cr \tab Name of column (character) in \code{sample_file_in} that corresponds to columns of \code{expression}. \cr
#'   }
#' @param n_samples_min minimum number of samples per feature; numeric >= 2.
#' @param n_features_min minimum number of features per sample; numeric >= 2.
#' @return A list (the filtered state) with the following elements:
#'   \tabular{ll}{
#'     \code{expression} \cr \tab Numeric matrix with non-negative expression values. \cr
#'     \code{features}   \cr \tab A data.frame with feature meta-data for rows of expression. \cr
#'     \code{samples}    \cr \tab A data.frame with observation meta-data for columns of expression. \cr
#'   } 
#' @examples
#' set.seed(101)
#' exprs <- h0testr::sim1(n_obs=6, n_feats=12, mcar_p=0.75)$mat
#' feats <- data.frame(feature_id=rownames(exprs))
#' samps <- data.frame(observation_id=colnames(exprs))
#' state <- list(expression=exprs, features=feats, samples=samps)
#' config <- list(feat_col="feature_id", obs_col="observation_id")
#' state2 <- h0testr::prefilter(state, config)
#' print(state2)

prefilter <- function(state, config, n_samples_min=2, n_features_min=2) {

  check_config(config)
  
  f.log_block("prefilter features and samples", config=config)
  f.msg("before filtering features:", config=config)
  f.report_state(state, config)
  
  state <- filter_features(state, config, n_samples_min=n_samples_min)
  f.msg("after filtering features", config=config)
  f.report_state(state, config)
  
  state <- filter_observations(state, config, n_features_min=n_features_min)
  f.msg("after filtering samples", config=config)
  f.report_state(state, config)
  
  f.check_state(state, config)
  
  return(state)
}

#' Add filter statistics
#' @description
#'   Adds filtering-related statistics to \code{state$features}, 
#'     and \code{state$samples}.
#' @details 
#'   Wrapper for \code{samples_per_feature()}, \code{feature_median_expression()}, 
#'     \code{features_per_sample()}. Also reports quantiles of distributions. 
#'   See documentation for \code{h0testr::new_config()} 
#'     for more detailed description of configuration parameters. 
#' @param state List with elements formatted like the list returned by \code{read_data()}:
#'   \tabular{ll}{
#'     \code{expression} \cr \tab Numeric matrix with non-negative expression values. \cr
#'     \code{features}   \cr \tab A data.frame with feature meta-data for rows of expression. \cr
#'     \code{samples}    \cr \tab A data.frame with observation meta-data for columns of expression. \cr
#'   } 
#' @param config List with configuration values. Requires the following keys:
#'   \tabular{ll}{
#'     \code{feat_col}  \cr \tab Name of column (character) in \code{feature_file_in} that corresponds to rows of \code{data_file_in}. \cr
#'     \code{obs_col}   \cr \tab Name of column (character) in \code{sample_file_in} that corresponds to columns of \code{expression}. \cr
#'   }
#' @return A list (the processed state) with the following elements:
#'   \tabular{ll}{
#'     \code{expression} \cr \tab Numeric matrix with non-negative expression values. \cr
#'     \code{features}   \cr \tab A data.frame with feature meta-data for rows of expression. \cr
#'     \code{samples}    \cr \tab A data.frame with observation meta-data for columns of expression. \cr
#'   } 
#' @examples
#' exprs <- h0testr::sim1(n_obs=6, n_feats=8)$mat
#' feats <- data.frame(feature_id=rownames(exprs))
#' samps <- data.frame(observation_id=colnames(exprs))
#' state <- list(expression=exprs, features=feats, samples=samps)
#' config <- h0testr::new_config()      ## defaults
#' config$save_state <- FALSE             ## default is TRUE
#' config$feat_col <- config$feat_id_col
#' config$obs_col <- config$obs_id_col
#' state <- h0testr::add_filter_stats(state, config)
#' print(state)

add_filter_stats <- function(state, config) {
  
  check_config(config)
  
  n <- samples_per_feature(state, config)
  if(!all(names(n) == state$features[[config$feat_col]])) {
    f.err("add_filter_stats: !all(names(n) == state$features[[config$feat_col]])", 
      config=config)
  }
  state$features[, config$n_samples_expr_col] <- n
  
  m <- feature_median_expression(state, config)
  if(!all(names(m) == state$features[[config$feat_col]])) {
    f.err("add_filter_stats: !all(names(m) == state$features[[config$feat_col]])", 
      config=config)
  }
  state$features[, config$median_raw_col] <- m
  
  n <- features_per_sample(state, config)
  if(!all(names(n) == state$samples[[config$obs_col]])) {
    f.err("add_filter_stats: !all(names(n) == state$samples[, config$obs_col])", 
      config=config)
  } 
  state$samples[, config$n_features_expr_col] <- n
  
  n <- apply(state$expression, 1, function(v) sum(v > 0, na.rm=T))
  f.msg("samples per feature:", config=config)
  f.quantile(n, config, digits=0)
  
  n <- apply(state$expression, 2, function(v) sum(v > 0, na.rm=T))
  f.msg("features per sample", config=config)
  f.quantile(n, config, digits=0)
  
  return(state)
}

#' Filter features and samples
#' @description
#'   Filter features and samples based on expression. 
#' @details 
#'   Filters out features with too few expressing samples, and filters out 
#'     samples with too few expressed features. Features are considered to be 
#'     expressed if \code{state$expression > 0}; \code{NA}s count as no expression.
#'   Features and/or samples considered constant if 
#'     \code{length(unique(expression_values)) \%in\% 1}.
#'   See documentation for \code{h0testr::new_config()} 
#'     for more detailed description of configuration parameters. 
#' @param state A list with elements like that returned by `read_data()`:
#'   \tabular{ll}{
#'     \code{expression} \cr \tab Numeric matrix with non-negative expression values. \cr
#'     \code{features}   \cr \tab A data.frame with feature meta-data for rows of expression. \cr
#'     \code{samples}    \cr \tab A data.frame with observation meta-data for columns of expression. \cr
#'   }
#' @param config List with configuration values. Uses the following keys:
#'   \tabular{ll}{
#'     \code{feat_col}            \cr \tab Name of column in \code{state$features} matching \code{rownames(state$expression)}.
#'     \code{obs_col}             \cr \tab Name of column in \code{state$samples} matching \code{colnames(state$expression)}.
#'     \code{n_features_min}      \cr \tab Minimum number (non-negative numeric) of features expressed in observation to keep observation. \cr
#'     \code{n_samples_min}       \cr \tab Minimum number (non-negative numeric) of samples expressing feature to keep feature. \cr
#'     \code{median_raw_col}      \cr \tab Name (character) of new column in feature metadata to hold median expression in expressing samples. \cr
#'     \code{n_samples_expr_col}  \cr \tab Name (character) of new column in feature metadata to hold number of expressing samples. \cr
#'     \code{n_features_expr_col} \cr \tab Name (character) of new column in sample metadata to hold number of expressed features. \cr
#'   }
#' @param remove_constant Logical scalar: if constant rows and columns of \code{state$expression} should be removed.
#' @return A list with elements like that returned by \code{read_data()}:
#'   \tabular{ll}{
#'     \code{expression} \cr \tab Numeric matrix with non-negative expression values. \cr
#'     \code{features}   \cr \tab A data.frame with feature meta-data for rows of expression. \cr
#'     \code{samples}    \cr \tab A data.frame with observation meta-data for columns of expression. \cr
#'   }
#' @examples
#' set.seed(101)
#' exprs <- h0testr::sim1(n_obs=6, n_feats=12, mcar_p=0.25)$mat
#' feats <- data.frame(feature_id=rownames(exprs))
#' samps <- data.frame(observation_id=colnames(exprs))
#' state <- list(expression=exprs, features=feats, samples=samps)
#'
#' ## assume default median_raw_col, n_samples_expr_col, and n_features_expr_col are ok:
#' config <- h0testr::new_config()      ## defaults
#' config$save_state <- FALSE             ## default is TRUE
#' config$feat_col <- config$feat_id_col
#' config$obs_col <- config$obs_id_col
#' config$n_features_min <- 6
#' config$n_samples_min <- 2
#' out <- h0testr::filter(state, config)
#' print(out$state)
#' str(out$config)

filter <- function(state, config, remove_constant=TRUE) {
  
  check_config(config)
  
  state <- filter_features(state, config, 
    remove_constant=remove_constant)  
  state <- filter_observations(state, config, 
    remove_constant=remove_constant)
  
  n <- samples_per_feature(state, config)
  state$features[[config$n_samples_expr_col]] <- n
  
  m <- feature_median_expression(state, config)
  state$features[[config$median_raw_col]] <- m
  
  n <- features_per_sample(state, config)
  state$samples[[config$n_features_expr_col]] <- n
  
  f.check_state(state, config)
  f.report_state(state, config)
  
  if(!is.null(config$run_order)) {
    i <- config$run_order %in% "filter"
    if(any(i)) {
      prfx <- paste0(which(i)[1] + 2, ".filtered")
    } else {
      prfx <- "filtered"
    }
  }
  f.save_state(state, config, prefix=prfx)
  
  return(list(state=state, config=config))
}
