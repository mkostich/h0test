## Helper to check that imputed values are usable:

f.pos_mat <- function(mat, config, is_log_transformed=NULL, fn_name="f.pos_mat") {

  is_log_transformed <- f.is_log_transformed(is_log_transformed, config,
    fn_name)

  if(!is_log_transformed) {

    i <- mat <= 0
    i[is.na(i)] <- F

    if(any(c(i))) {
      f.err(fn_name, ":", sum(i), "of", length(mat), "imputed values are at or",
        "below zero on the raw scale, where a value at or below zero cannot be",
        "a measurement;", "\n",
        "  min:", min(c(mat), na.rm=T), "; imputation cannot stand in a",
        "measurement here, so the values are not usable;", "\n",
        "  try a different config$impute_method, or normalize() first so that",
        "imputation happens on the log scale where these values are ordinary",
        config=config)
    }
  }

  ## infinite values unusable on either scale; reported rather than replaced:
  i <- is.infinite(mat)
  i[is.na(i)] <- F

  if(any(c(i))) {

    if(is.matrix(mat)) {
      pos <- which(i, arr.ind=T)
      rows <- rownames(mat)[pos[, 1]]
      if(is.null(rownames(mat))) rows <- as.character(pos[, 1])
      cols <- colnames(mat)[pos[, 2]]
      if(is.null(colnames(mat))) cols <- as.character(pos[, 2])
      where <- paste0("  features: ", length(unique(rows)), " affected, first: ",
        paste(utils::head(unique(rows), 3), collapse=", "),
        "; observations: ", length(unique(cols)), " affected, first: ",
        paste(utils::head(unique(cols), 3), collapse=", "))
    } else {
      where <- paste("  positions:",
        paste(utils::head(which(i), 5), collapse=", "))
    }

    f.err(fn_name, ":", sum(i), "of", length(mat), "imputed values are",
      "infinite, so the imputer had nothing to work from where they came;",
      "\n",
      "  an infinite value is not a measurement and is not substituted for:",
      "a missing value is a below-detection value, so a usable substitute",
      "would sit at the bottom of the scale, which is what the imputer could",
      "not supply;", "\n",
      where, "\n",
      "  the usual cause is a feature or observation with no measured value at",
      "all; drop those with h0testr::filter_state(), or try a different",
      "config$impute_method", config=config)
  }

  return(mat)
}

#' Impute missing values below a global LOD
#' @description
#'   Impute missing values by randomly drawing from uniform distribution
#'     below an estimated global LOD.
#' @details Imputed values are random draws from a uniform distribution over the
#'   interval \code{[floor, LOD]}, where \code{LOD} is quantile
#'   \code{impute_quantile} of the per-feature minima, an estimate of the global
#'   limit of detection. The \code{floor} depends on the scale of the data. For
#'   raw input (\code{config$is_log_transformed} \code{FALSE}) it is \code{0},
#'   zero abundance. For log input it is
#'   \code{min(state$expression, na.rm=TRUE) + config$impute_floor_offset}, since
#'   zero on a log scale is a single count rather than the bottom of the scale,
#'   and typically sits near the top of the range after normalization; the offset
#'   is also what gives the interval width when \code{impute_quantile} is
#'   \code{0}. An interval with no width is an error rather than a silent
#'   \code{NaN}. Only \code{NA} values are considered as missing, so if you want
#'   \code{0} to be considered missing, and have \code{0} in the data, do
#'   something like \code{state$expression[state$expression \%in\% 0] <- NA}
#'   prior to imputing.
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
#'     \code{impute_quantile}     \cr \tab Quantile of signal distribution to use as estimated limit of detection (LOD). \cr
#'     \code{is_log_transformed}  \cr \tab Whether the data are on a log-like scale; sets the lower bound of the draw. \cr
#'     \code{impute_floor_offset} \cr \tab Offset from the smallest measured value giving that lower bound, for log data. \cr
#'   }
#' @param impute_quantile Numeric between 0 and 1 specifying
#'   minimum non-zero/\code{NA} expression value for each feature to use as global
#'   LOD (limit of detection). If \code{NULL}, \code{config$impute_quantile} used.
#'   Default: \code{0}.
#' @return An updated `state` list with the following elements:
#'   \tabular{ll}{
#'     \code{expression} \cr \tab Numeric matrix with expression values. \cr
#'     \code{features}   \cr \tab A data.frame with feature meta-data for rows of expression. \cr
#'     \code{samples}    \cr \tab A data.frame with observation meta-data for columns of expression. \cr
#'   }
#' Returned \code{state$expression} matrix contains strictly positive values for
#'   raw input, and may contain values at or below zero for log input, where such
#'   a value is an ordinary dim measurement.
#' @examples
#' set.seed(101)
#' exprs <- h0testr::sim1(n_obs=6, n_feats=12)$mat
#' exprs <- log2(exprs + 1)
#' feats <- data.frame(feature_id=rownames(exprs))
#' samps <- data.frame(observation_id=colnames(exprs))
#' state <- list(expression=exprs, features=feats, samples=samps)
#' ## the data were log transformed above, so say so; imputed values are then
#' ##   allowed to be negative, as an ordinary small log intensity is:
#' config <- list(impute_quantile=0.05, is_log_transformed=TRUE)
#' state2 <- h0testr::impute_unif_global_lod(state, config)
#' summary(c(state$expression))   ## note number of NAs
#' summary(c(state2$expression))  ## note number of NAs
#' head(state$expression)
#' round(head(state2$expression))

impute_unif_global_lod <- function(state, config, impute_quantile=NULL) {

  check_config(config)

  if(!is.matrix(state$expression)) {
     f.err("impute_unif_global_lod: !is.matrix(state$expression)", "\n",
      "class(state$expression):", class(state$expression), config=config)
  }
  if(is.null(impute_quantile)) impute_quantile <- config$impute_quantile
  if(is.null(impute_quantile)) impute_quantile <- 0
  
  ## smallest value measured for each feature:
  v <- suppressWarnings(apply(state$expression, 1, min, na.rm=T))
  v <- v[is.finite(v)]

  if(!length(v)) {
    f.err("impute_unif_global_lod: no feature has a finite measured value, so",
      "there is no limit of detection to estimate;", "\n",
      "  features:", nrow(state$expression), "; observations:",
      ncol(state$expression), config=config)
  }

  max_val <- stats::quantile(v, probs=impute_quantile[1], na.rm=T)    ## quantile of min vals
  names(max_val) <- NULL

  ## bottom of the scale: 
  floor_val <- f.impute_floor(state$expression, config, "impute_unif_global_lod")

  if(max_val <= floor_val) {
    f.err("impute_unif_global_lod: the imputation interval has no width;",
      "lower bound:", floor_val, "; upper bound:", max_val, ";", "\n",
      "  the upper bound is quantile", impute_quantile[1], "of the per-feature",
      "minima, and the lower bound is zero for raw data, or the smallest",
      "measured value plus config$impute_floor_offset for log data;", "\n",
      "  to fix, either raise config$impute_quantile, or make",
      "config$impute_floor_offset more negative (log data only), or check that",
      "config$is_log_transformed describes the data: raw data with values at or",
      "below zero lands here", config=config)
  }

  f.msg("impute_unif_global_lod: impute_quantile:", impute_quantile[1],
    "; drawing from [", floor_val, ",", max_val, "]; width:",
    max_val - floor_val, config=config)

  mat <- state$expression
  i <- is.na(mat)
  if(any(i)) mat[i] <- stats::runif(sum(i), floor_val, max_val)
  mat <- f.pos_mat(mat, config, fn_name="impute_unif_global_lod")
  state$expression <- mat
  
  return(state)
}

#' Impute missing values below a sample LOD
#' @description
#'   Impute missing values by randomly drawing from uniform distribution below
#'     an estimated observation-specific limit of detection (LOD).
#' @details Imputed values are random draws from a uniform distribution over the
#'   interval \code{[floor, LOD]}, where \code{LOD} is quantile
#'   \code{impute_quantile} of the values in that observation, an estimate of the
#'   observation's limit of detection. The \code{floor} is shared by all
#'   observations, since the bottom of the scale is a property of the data as a
#'   whole: \code{0} for raw input, and
#'   \code{min(state$expression, na.rm=TRUE) + config$impute_floor_offset} for
#'   log input, where zero is a single count rather than the bottom of the scale.
#'   An interval with no width is an error rather than a silent \code{NaN}. Only
#'   \code{NA} values are considered as missing, so if you
#'   want \code{0} to be considered missing, and have \code{0} in the data, do
#'   something like \code{exprs[exprs \%in\% 0] <- NA} prior to imputing.
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
#'     \code{impute_quantile}     \cr \tab Quantile of signal distribution to use as estimated limit of detection (LOD). \cr
#'     \code{is_log_transformed}  \cr \tab Whether the data are on a log-like scale; sets the lower bound of the draw. \cr
#'     \code{impute_floor_offset} \cr \tab Offset from the smallest measured value giving that lower bound, for log data. \cr
#'   }
#' @param impute_quantile Numeric between 0 and 1 specifying
#'   minimum non-zero/\code{NA} expression value in each observation to use as the
#'   observation-specific LOD (limit of detection). Default: \code{0}.
#' @return An updated \code{state} list with the following elements:
#'   \tabular{ll}{
#'     \code{expression} \cr \tab Numeric matrix with expression values. \cr
#'     \code{features}   \cr \tab A data.frame with feature meta-data for rows of expression. \cr
#'     \code{samples}    \cr \tab A data.frame with observation meta-data for columns of expression. \cr
#'   } 
#' @examples
#' set.seed(101)
#' exprs <- h0testr::sim1(n_obs=6, n_feats=12)$mat
#' exprs <- log2(exprs + 1)
#' feats <- data.frame(feature_id=rownames(exprs))
#' samps <- data.frame(observation_id=colnames(exprs))
#' state <- list(expression=exprs, features=feats, samples=samps)
#' ## the data were log transformed above, so say so; imputed values are then
#' ##   allowed to be negative, as an ordinary small log intensity is:
#' config <- list(impute_quantile=0.05, is_log_transformed=TRUE)
#' state2 <- h0testr::impute_unif_sample_lod(state, config)
#' summary(c(state$expression))   ## note number of NAs
#' summary(c(state2$expression))  ## note number of NAs
#' head(state$expression)
#' round(head(state2$expression))

impute_unif_sample_lod <- function(state, config, impute_quantile=NULL) {

  check_config(config)

  if(!is.matrix(state$expression)) {
    f.err("impute_unif_sample_lod: !is.matrix(state$expression)", "\n",
      "class(state$expression):", class(state$expression), config=config)
  }
  if(is.null(impute_quantile)) impute_quantile <- config$impute_quantile
  if(is.null(impute_quantile)) impute_quantile <- 0
  
  i_all <- apply(state$expression, 2, function(v) all(is.na(v)))

  if(any(i_all)) {
    nom <- colnames(state$expression)[i_all]
    if(is.null(nom)) nom <- as.character(which(i_all))
    f.err("impute_unif_sample_lod:", sum(i_all), "of", length(i_all),
      "observations have no measured value, so they have no limit of detection",
      "to impute against;", "\n",
      "  observations:", paste(utils::head(nom, 5), collapse=", "),
      config=config)
  }

  ## one upper bound per observation, but single lower bound:
  max_vals <- apply(state$expression, 2, stats::quantile,
    probs=impute_quantile[1], na.rm=T)
  names(max_vals) <- colnames(state$expression)

  floor_val <- f.impute_floor(state$expression, config,
    "impute_unif_sample_lod")

  i_flat <- max_vals <= floor_val

  if(any(i_flat)) {
    nom <- names(max_vals)[i_flat]
    if(is.null(nom)) nom <- as.character(which(i_flat))
    f.err("impute_unif_sample_lod: the imputation interval has no width for",
      sum(i_flat), "of", length(max_vals), "observations; lower bound:",
      floor_val, "; smallest upper bound:", min(max_vals), ";", "\n",
      "  each upper bound is quantile", impute_quantile[1], "of its own",
      "observation, and the lower bound is zero for raw data, or the smallest",
      "measured value plus config$impute_floor_offset for log data;", "\n",
      "  to fix, either raise config$impute_quantile, or make",
      "config$impute_floor_offset more negative (log data only), or check that",
      "config$is_log_transformed describes the data: raw data with values at or",
      "below zero lands here;", "\n",
      "  observations:", paste(utils::head(nom, 5), collapse=", "),
      config=config)
  }

  f.msg("impute_unif_sample_lod: impute_quantile:", impute_quantile[1],
    "; drawing from [", floor_val, ", upper bound ]; upper bounds range over [",
    min(max_vals), ",", max(max_vals), "]; widths range over [",
    min(max_vals) - floor_val, ",", max(max_vals) - floor_val, "]",
    config=config)

  mat <- state$expression

  for(k in seq_len(ncol(mat))) {
    i <- is.na(mat[, k])
    if(any(i)) mat[i, k] <- stats::runif(sum(i), floor_val, max_vals[k])
  }

  mat <- f.pos_mat(mat, config, fn_name="impute_unif_sample_lod")
  state$expression <- mat
  
  return(state)
}

#' Impute missing values as the sample LOD
#' @description
#'   Impute missing values as the minimum observed value in the corresponding 
#'     observation.
#' @details Missing values are set to an estimate observation limit of
#'   detection (LOD). Only \code{NA} values are considered as missing, so if you
#'   want \code{0} to be considered missing, and have \code{0} in the data, do
#'   something like \code{exprs[exprs \%in\% 0] <- NA} prior to imputing.
#'   An observation with no measured value has no minimum to stand in, so it
#'     cannot be imputed and is an error rather than a column of \code{Inf}.
#'     \code{h0testr::filter_state()} removes such observations, so this normally only
#'     arises when the imputer is called on its own.
#' @param state A list with elements like that returned by \code{read_data()}:
#'   \tabular{ll}{
#'     \code{expression} \cr \tab Numeric matrix with non-negative expression values. \cr
#'     \code{features}   \cr \tab A data.frame with feature meta-data for rows of expression. \cr
#'     \code{samples}    \cr \tab A data.frame with observation meta-data for columns of expression. \cr
#'   }
#' @param config List with configuration values. Uses the following key:
#'   \tabular{ll}{
#'     \code{is_log_transformed} \cr \tab Whether the data are on a log-like scale; decides whether an imputed value at or below zero is an error. \cr
#'   }
#'   This function has no \code{is_log_transformed} argument of its own, the
#'     minimum of an observation being its minimum on either scale, so the key
#'     is the only place the scale can come from and must be set.
#'     \code{h0testr::init_state()} sets it and \code{h0testr::normalize()}
#'     updates it.
#' @return An updated \code{state} list with the following elements:
#'   \tabular{ll}{
#'     \code{expression} \cr \tab Numeric matrix with non-negative expression values. \cr
#'     \code{features}   \cr \tab A data.frame with feature meta-data for rows of expression. \cr
#'     \code{samples}    \cr \tab A data.frame with observation meta-data for columns of expression. \cr
#'   } 
#' @examples
#' set.seed(101)
#' exprs <- h0testr::sim1(n_obs=6, n_feats=12)$mat
#' exprs <- log2(exprs + 1)
#' feats <- data.frame(feature_id=rownames(exprs))
#' samps <- data.frame(observation_id=colnames(exprs))
#' state <- list(expression=exprs, features=feats, samples=samps)
#' ## the data were log transformed above, so say so; there is no
#' ##   is_log_transformed argument here, so the config is what says it:
#' config <- list(is_log_transformed=TRUE)
#' state2 <- h0testr::impute_sample_lod(state, config)
#' summary(c(state$expression))   ## note number of NAs
#' summary(c(state2$expression))  ## note number of NAs
#' head(state$expression)
#' round(head(state2$expression))

impute_sample_lod <- function(state, config) {

  check_config(config)

  if(!is.matrix(state$expression)) {
    f.err("impute_sample_lod: !is.matrix(state$expression)", "\n",
      "class(state$expression):", class(state$expression), config=config)
  }

  ## written a column at a time rather than through apply(., 2, ...);
  ##   because latter returns plain vector when matrix has single feature:

  mat <- state$expression

  for(k in seq_len(ncol(mat))) {
    i <- is.na(mat[, k])
    if(any(i)) mat[i, k] <- min(mat[!i, k])
  }

  mat <- f.pos_mat(mat, config, fn_name="impute_sample_lod")
  state$expression <- mat

  return(state)
}

#' Impute missing values near the feature mean
#' @description
#'   Impute missing values normally distributed around the feature mean.
#' @details Missing values are draws from 
#'   \code{normal(mean=mean(exprs[feature, ]), sd=(scale. * sd(exprs[feature, ])))}. 
#'   Only non-\code{NA} values are used in calculation of \code{mean} and 
#'   \code{sd}. If you want \code{0} to be considered missing, and have 
#'   \code{0} in the data, do something like \code{exprs[exprs \%in\% 0] <- NA}
#'   prior to imputing.
#'   The dispersion has a lower bound, for a feature with fewer than two observed
#'     values: \code{sqrt(mean)} on the raw scale, a tenth of the whole matrix
#'     standard deviation on the log scale. Raw scale draws are resampled until
#'     non-negative; log scale draws are kept, a value below zero being ordinary there.
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
#'     \code{impute_scale}  \cr \tab Factor (numeric) rescaling variance of normal distribution from which draws are made. See details.
#'     \code{is_log_transformed} \cr \tab Whether the data are on a log-like scale; sets the dispersion floor and whether a value below zero is allowed. See details.
#'   }
#' @param is_log_transformed Logical scalar: whether \code{state$expression} has
#'   been log transformed. Defaults to \code{config$is_log_transformed}, which
#'   \code{h0testr::init_state()} and \code{h0testr::normalize()} maintain;
#'   passing both is an error unless they agree.
#' @param scale. Numeric greater than zero, linearly scaling the
#'   dispersion around the feature mean. Default: \code{1.0}.
#' @return An updated \code{state} list with the following elements:
#'   \tabular{ll}{
#'     \code{expression} \cr \tab Numeric matrix with non-negative expression values. \cr
#'     \code{features}   \cr \tab A data.frame with feature meta-data for rows of expression. \cr
#'     \code{samples}    \cr \tab A data.frame with observation meta-data for columns of expression. \cr
#'   } 
#' @examples
#' set.seed(101)
#' exprs <- h0testr::sim1(n_obs=6, n_feats=12)$mat
#' exprs <- log2(exprs + 1)
#' feats <- data.frame(feature_id=rownames(exprs))
#' samps <- data.frame(observation_id=colnames(exprs))
#' state <- list(expression=exprs, features=feats, samples=samps)
#' 
#' summary(c(state$expression))     ## note number of NAs
#' head(state$expression)
#'
#' ## impute using default scale. parameter:
#' config <- list(is_log_transformed=TRUE)   ## the data were log transformed above
#' state2 <- h0testr::impute_rnorm_feature(state, config)
#' summary(c(state2$expression))    ## note number of NAs
#' round(head(state2$expression))
#'
#' ## impute using passed scale. parameter:
#' state2 <- h0testr::impute_rnorm_feature(state, config, scale.=2)
#' summary(c(state2$expression))    ## note number of NAs
#' round(head(state2$expression))
#'
#' ## impute using impute_scale parameter from config:
#' config <- list(impute_scale=0.5, is_log_transformed=TRUE)
#' state2 <- h0testr::impute_rnorm_feature(state, config)
#' summary(c(state2$expression))    ## note number of NAs
#' round(head(state2$expression))
#'
#' ## raw scale input; here a value below zero would be an error:
#' state$expression <- 2^state$expression - 1
#' config$is_log_transformed <- FALSE   ## argument and config must agree
#' state2 <- h0testr::impute_rnorm_feature(state, config, is_log_transformed=FALSE)
#' summary(c(state2$expression))    ## note number of NAs

impute_rnorm_feature <- function(state, config, is_log_transformed=NULL,
    scale.=NULL) {

  check_config(config)

  if(!is.matrix(state$expression)) {
    f.err("impute_rnorm_feature: !is.matrix(state$expression)", "\n",
      "class(state$expression):", class(state$expression), config=config)
  }

  is_log_transformed <- f.is_log_transformed(is_log_transformed, config,
    "impute_rnorm_feature")

  if(!is_log_transformed) {
    i <- c(state$expression) < 0
    i[is.na(i)] <- F
    if(any(i)) {
      f.err("impute_rnorm_feature: state$expression contains negative values on the",
        "raw scale, where a value below zero cannot be a measurement;", "\n",
        "  min:", min(c(state$expression), na.rm=T), "; if these are log scale",
        "values, set config$is_log_transformed to TRUE", config=config)
    }
  }

  if(is.null(scale.)) scale. <- config$impute_scale
  if(is.null(scale.)) scale. <- 1.0

  s_all <- stats::sd(c(state$expression), na.rm=T)
  if(!is.finite(s_all) || s_all <= 0) s_all <- 1

  f <- function(v) {
    i <- is.na(v)
    if(all(i)) {
      f.err("impute_rnorm_feature: all(is.na(state$expression[row, ]))",
        config=config)
    }
    if(any(i)) {
      m <- mean(v[!i])
      ## sqrt(mean) is the raw scale dispersion and means nothing on the log scale:
      if(is_log_transformed) {
        s_min <- 0.1 * s_all
      } else {
        s_min <- sqrt(m)
      }
      if(sum(!i) >= 2) {
        s <- stats::sd(v[!i]) * scale.
      } else {
        s <- s_min
      }
      if(s < s_min) s <- s_min

      v[i] <- stats::rnorm(sum(i), mean=m, sd=s)

      ## resampling on the log scale would bias the imputed values upward:
      if(!is_log_transformed) {
        i <- v < 0
        i[is.na(i)] <- F
        while(any(i)) {
          v[i] <- stats::rnorm(sum(i), mean=m, sd=s)
          i <- v < 0
          i[is.na(i)] <- F
        }
      }
    }
    return(v)
  }
  state$expression <- t(apply(state$expression, 1, f))
  state$expression <- f.pos_mat(state$expression, config, is_log_transformed,
    "impute_rnorm_feature")

  return(state)
}

#' Impute by drawing from \code{p(missing|intensity)} from binomial glm
#' @description
#'   Impute missing values by randomly drawing from an estimated density 
#'     of \code{p(missing|intensity)}.
#' @details 
#'   Impute random draws by drawing from the estimated density of 
#'     \code{p(missing|intensity)}. Estimates \code{p(missing|intensity)} using 
#'     a binomial glm with formula \code{cbind(n.missing, n.found) ~ intensity},
#'     where intensity is assumed to be previously log transformed. If you  
#'     want \code{0} to be considered missing, and have \code{0} in the data, do 
#'     something like \code{exprs[exprs \%in\% 0] <- NA} prior to imputing. 
#'   See documentation for \code{h0testr::new_config()} 
#'     for more detailed description of configuration parameters. 
#' @param state A list with elements like that returned by \code{read_data()}:
#'   \tabular{ll}{
#'     \code{expression} \cr \tab Numeric matrix with non-negative expression values. \cr
#'     \code{features}   \cr \tab A data.frame with feature meta-data for rows of expression. \cr
#'     \code{samples}    \cr \tab A data.frame with observation meta-data for columns of expression. \cr
#'   } 
#' It is assumed that state$expression has been previously \code{log2(x+1)} transformed.
#' @param config List with configuration values. Uses the following keys:
#'   \tabular{ll}{
#'     \code{impute_n_pts}  \cr \tab Numeric greater than one. Determines granularity of imputation. Larger values lead to finer grain.
#'   }
#' @param is_log_transformed Logical scalar: whether \code{state$expression} has
#'   been log transformed. Defaults to \code{config$is_log_transformed}, which
#'   \code{h0testr::init_state()} and \code{h0testr::normalize()} maintain;
#'   passing both is an error unless they agree.
#' @param n_pts Numeric greater than one. Granularity of prediction 
#'   grid. Larger values lead to less chance of duplicate imputed values.
#'   Larger values require more compute time and memory. Default: \code{1e7}.
#' @param off Numeric offset for calculating 
#'   \code{p.missing = (n.missing + off) / (n.total + off)}.
#' @param f_mid Function to use for calculating central tendency of
#'   feature expression across samples.
#' @param min_fit_pts Numeric scalar (at least two): the number of distinct
#'   feature intensities that must remain after features measured in no sample
#'   have been dropped, before the missingness fit is attempted. Stops with an 
#'   error when too few distinct intensities remain. Default: \code{10}.
#' @return An updated `state` list with the following elements:
#'   \tabular{ll}{
#'     \code{expression} \cr \tab Numeric matrix with non-negative expression values. \cr
#'     \code{features}   \cr \tab A data.frame with feature meta-data for rows of expression. \cr
#'     \code{samples}    \cr \tab A data.frame with observation meta-data for columns of expression. \cr
#'   }
#' @examples
#' set.seed(101)
#' exprs <- h0testr::sim1(n_obs=6, n_feats=12)$mat
#' feats <- data.frame(feature_id=rownames(exprs))
#' samps <- data.frame(observation_id=colnames(exprs))
#' state <- list(expression=exprs, features=feats, samples=samps)
#' config <- list(impute_n_pts=1e7)
#' 
#' ## untransformed example:
#' state2 <- h0testr::impute_glm_binom(state, config, is_log_transformed=FALSE)
#' summary(c(state$expression))    ## Note number of NAs
#' summary(c(state2$expression))   ## Note number of NAs
#' head(state$expression)
#' round(head(state2$expression))
#' 
#' ## log-transformed example:
#' state$expression <- log2(state$expression + 1)
#' state2 <- h0testr::impute_glm_binom(state, config, is_log_transformed=TRUE)
#' summary(c(state$expression))    ## Note number of NAs
#' summary(c(state2$expression))   ## Note number of NAs
#' head(state$expression)
#' round(head(state2$expression))

impute_glm_binom <- function(state, config, is_log_transformed=NULL,
    n_pts=NULL, off=1, f_mid=stats::median, min_fit_pts=10) {
  
  check_config(config)
  
  if(!is.matrix(state$expression)) {
    f.err("impute_glm_binom: !is.matrix(state$expression)", "\n",
      "class(state$expression):", class(state$expression), config=config)
  }
  
  is_log_transformed <- f.is_log_transformed(is_log_transformed, config,
    "impute_glm_binom")
  
  mat <- state$expression
  if(!is_log_transformed) mat <- log2(mat + 1) 
  
  if(is.null(n_pts)) n_pts <- config$impute_n_pts
  if(is.null(n_pts)) n_pts <- 1e7
  if(n_pts <= 0) {
    f.err("impute_glm_binom: n_pts <= 0; n_pts:", n_pts, config=config)
  }
  
  ## NA is only indicator of missing value:

  m <- apply(mat, 1, f_mid, na.rm=T)

  n0 <- apply(mat, 1, function(v) sum(is.na(v)))             ## n.missing
  n1 <- ncol(mat) - n0
  p <- (n0 + off) / (ncol(mat) + off)
  dat <- data.frame(n0=n0, n1=n1, p=p, m=m)
  dat <- f.drop_unfittable(dat, config, "impute_glm_binom", min_fit_pts)

  fit <- stats::glm(cbind(n0, n1) ~ m, data=dat, family="binomial")
  
  m_new <- seq(from=max(c(mat), na.rm=T) / n_pts, 
    to=max(c(mat), na.rm=T), length.out=n_pts)
  
  p_hat <- stats::predict(fit, newdata=data.frame(m=m_new), type="response")
  p_hat[is.na(p_hat)] <- min(p_hat, na.rm=T)
  i_na <- is.na(c(mat))
  mat[i_na] <- sample(m_new, sum(i_na), replace=T, prob=p_hat)
  
  if(!is_log_transformed) mat <- 2^mat
  mat <- f.pos_mat(mat, config, is_log_transformed, "impute_glm_binom")
  state$expression <- mat
  
  return(state)
}

#' Impute by drawing from \code{p(missing|intensity)} estimated with loess
#' @description
#'   Impute missing values by randomly drawing from an estimated density 
#'     of \code{p(missing|intensity)}.
#' @details 
#'   Impute random draws by drawing from the estimated density of 
#'     \code{p(missing|intensity)}. Estimates \code{p(missing|intensity)} using 
#'     a loess fit with form \code{log(p_missing/p_found) ~ intensity}, where 
#'     intensity is assumed to be previously log transformed. If you want 
#'     \code{0} to be considered missing, and have \code{0} in the data, do 
#'     something like \code{exprs[exprs \%in\% 0] <- NA} prior to imputing. 
#'   See documentation for \code{h0testr::new_config()} 
#'     for more detailed description of configuration parameters. 
#' @param state A list with elements like that returned by \code{read_data()}:
#'   \tabular{ll}{
#'     \code{expression} \cr \tab Numeric matrix with non-negative expression values. \cr
#'     \code{features}   \cr \tab A data.frame with feature meta-data for rows of expression. \cr
#'     \code{samples}    \cr \tab A data.frame with observation meta-data for columns of expression. \cr
#'   } 
#' It is assumed that state$expression has been previously \code{log2(x+1)} transformed.
#' @param config List with configuration values. Uses the following keys:
#'   \tabular{ll}{
#'     \code{impute_n_pts}  \cr \tab Numeric greater than one. Determines granularity of imputation. Larger values lead to finer grain.
#'     \code{impute_span}   \cr \tab Span (numeric between zero and 1) for \code{loess} fit.
#'   }
#' @param span Span for loess fit. Scalar with \code{0 < span < 1}. Default: \code{0.5}.
#' @param n_pts Numeric greater than one. Granularity of prediction 
#'   grid. Larger values lead to less chance of duplicate imputed values.
#'   Larger values require more compute time and memory. Default: \code{1e7}.
#' @param off Numeric offset for calculating 
#'   \code{p.missing = (n.missing + off) / (n.total + off)}.
#' @param f_mid Function to use for calculating central tendency of 
#'   feature expression across samples. 
#' @param degree Numeric in set \code{c(1, 2)}. Degree for loess fit.
#' @param fam Character in set \code{c("symmetric", "gaussian")}. Family
#'   for \code{loess} fit.
#' @param min_fit_pts Numeric scalar (at least two): the number of distinct
#'   feature intensities that must remain after features measured in no sample
#'   have been dropped, before the missingness fit is attempted. Those features
#'   have no intensity to fit against, so they are excluded from the fit and
#'   imputed from the curve the rest determine. Stops with an error when too few
#'   distinct intensities remain, since the fit is extrapolated across the whole
#'   intensity range and imputed values are drawn from it. Default: \code{10}.
#' @return An updated \code{state} list with the following elements:
#'   \tabular{ll}{
#'     \code{expression} \cr \tab Numeric matrix with non-negative expression values. \cr
#'     \code{features}   \cr \tab A data.frame with feature meta-data for rows of expression. \cr
#'     \code{samples}    \cr \tab A data.frame with observation meta-data for columns of expression. \cr
#'   } 
#' @examples
#' set.seed(101)
#' exprs <- h0testr::sim1(n_obs=6, n_feats=25)$mat
#' exprs <- log2(exprs + 1)
#' feats <- data.frame(feature_id=rownames(exprs))
#' samps <- data.frame(observation_id=colnames(exprs))
#' state <- list(expression=exprs, features=feats, samples=samps)
#' config <- list(impute_n_pts=1e7, impute_span=0.5)
#' state2 <- h0testr::impute_loess_logit(state, config)
#' summary(c(state$expression))     ## note number of NAs
#' summary(c(state2$expression))    ## note number of NAs
#' head(state$expression)
#' round(head(state2$expression))

impute_loess_logit <- function(state, config, span=NULL, n_pts=NULL,
    off=0.1, f_mid=stats::median, degree=1, fam="symmetric", min_fit_pts=10) {
  
  check_config(config)
  
  if(!is.matrix(state$expression)) {
    f.err("impute_loess_logit: !is.matrix(state$expression)", "\n",
      "class(state$expression):", class(state$expression), config=config)
  }
  if(is.null(span)) span <- config$impute_span
  if(is.null(span)) span <- 0.5
  if(is.null(n_pts)) n_pts <- config$impute_n_pts
  if(is.null(n_pts)) n_pts <- 1e7
  
  ## here NA means feature measured nowhere:

  m <- apply(state$expression, 1, f_mid, na.rm=T)
  if(!any(is.finite(m))) {
    f.err("impute_loess_logit: !any(is.finite(m))", "\n",
      "sort(m):", sort(m), config=config)
  }
  m[is.infinite(m)] <- max(m[is.finite(m)])
  
  ## NA is only indicator of missing value:

  f <- function(v) sum(is.na(v))
  n0 <- apply(state$expression, 1, f)                        ## n.missing
  if(all(n0 %in% 0)) return(state)                           ## nothing missing to model
  
  n0[is.na(n0)] <- ncol(state$expression)
  n1 <- ncol(state$expression) - n0                          ## n.found
  p <- (n0 + off) / (ncol(state$expression) + 2 * off)       ## p.missing
  logitp <- log(p / (1 - p))
  if(length(unique(logitp)) %in% 1) {
    f.err("impute_loess_logit: length(unique(logitp)) %in% 1", "\n", 
      "unique(logitp):", unique(logitp), config=config)
  }
  
  if(!all(is.finite(logitp))) {
    f.err("impute_loess_logit: !all(is.finite(logitp))", "\n",
      "sort(logitp):", sort(logitp), config=config)
  }
  dat <- data.frame(n0=n0, n1=n1, p=p, logitp=logitp, m=m)
  dat <- f.drop_unfittable(dat, config, "impute_loess_logit", min_fit_pts)

  fit <- stats::loess(logitp ~ m, data=dat, span=span,
    degree=degree, family=fam)
  
  m_new <- seq(from=max(c(state$expression), na.rm=T) / n_pts, 
    to=max(c(state$expression), na.rm=T), length.out=n_pts)
  
  p_hat <- stats::predict(fit, newdata=data.frame(m=m_new))  ## on logit scale
  p_hat = exp(p_hat) / (1 + exp(p_hat))                      ## inverse logit
  
  i_p <- !is.na(p_hat)
  i_p[is.na(i_p)] <- F
  if(any(i_p)) {
    p_hat[is.na(p_hat)] <- min(p_hat[i_p], na.rm=T) / 2      ## is.na -> low p
  } else {
    ## p_hat <- rep(1 / n_pts, length(p_hat))                ## maybe uniform?
    f.err("impute_loess_logit: all probabilities NA", 
      config=config)                                         ## for now
  }
  
  i_na <- is.na(c(state$expression))                         ## values to be imputed
  state$expression[i_na] <- sample(m_new, sum(i_na), replace=T, prob=p_hat)
  
  return(state)
}

## helper for augmenting training data with linearly transformed
##   (multiplied by mult) versions of itself;
##   mult: numeric multiplier > 0;
##   steps: number of augmentatons, each one x' = x * mult^step:

f.augment_affine <- function(exprs, mult=1, add=0, steps=1) {

  xs <- list()
  xs[[length(xs) + 1]] <- exprs

  while(steps > 0) {
    exprs <- (exprs * mult) + add
    xs[[length(xs) + 1]] <- exprs
    steps = steps - 1
  }
  exprs <- do.call(cbind, xs)
  
  return(exprs)
}

#' Impute missing values using \code{randomForest}
#' @description
#'   Impute missing values using the \code{randomForest} package.
#' @details Imputes missing values using a \code{randomForest::randomForest} 
#'   model trained using observations in which the feature was expressed. 
#'   Iterates through features beginning with those having fewest missing 
#'   values. If \code{aug_steps > 0}, augments observations with affine 
#'   transormed versions. This is meant to enable extrapolation outside of 
#'   observed intensity range. If \code{aug_steps > 0}, assumes expression 
#'   data have been previously \code{log(x+1)} transformed. If you want 
#'   \code{0} to be considered missing, and have \code{0} in the data, do 
#'   something like \code{exprs[exprs \%in\% 0] <- NA} prior to imputing. 
#'   Augmentation uses \code{f.augment_affine()}.
#' @param state A list with elements like that returned by \code{read_data()}:
#'   \tabular{ll}{
#'     \code{expression} \cr \tab Numeric matrix with non-negative expression values. \cr
#'     \code{features}   \cr \tab A data.frame with feature meta-data for rows of expression. \cr
#'     \code{samples}    \cr \tab A data.frame with observation meta-data for columns of expression. \cr
#'   } 
#' It is assumed that state$expression has been previously \code{log2(x+1)} transformed.
#' @param config List with configuration values. Does not use any keys, so can pass empty list.
#' @param is_log_transformed Logical scalar: whether \code{state$expression} has
#'   been log transformed. Defaults to \code{config$is_log_transformed}, which
#'   \code{h0testr::init_state()} and \code{h0testr::normalize()} maintain;
#'   passing both is an error unless they agree.
#' @param f_imp Function to use for initial rough imputation.
#' @param ntree Numeric (greater than 0) number of trees in random forest.
#'   \code{randomForest} \code{ntree} parameter.
#' @param mtry Numeric (greater than 0) number of features to sample at each 
#'   node in each tree. \code{randomForest} \code{mtry} parameter.
#' @param aug_mult Numeric affine shift for training data augmentation. 
#' @param aug_add Numeric affine shift for training data augmentation.
#' @param aug_steps Numeric (non-negative) number of augmentation steps. Set to
#'   \code{0} to skip augmentation. Default: \code{3}.
#' @param verbose Logical if TRUE, emits progress messages.
#' @return A list with the following elements:
#'   \tabular{llll}{
#'     \code{state} \cr 
#'       \tab A list with the following elements: \cr
#'       \tab \code{expression} \cr 
#'       \tab \tab \tab Numeric matrix with non-negative expression values. \cr
#'       \tab \code{features}   \cr 
#'       \tab \tab \tab A data.frame with feature meta-data for rows of expression. \cr
#'       \tab \code{samples}    \cr 
#'       \tab \tab \tab A data.frame with observation meta-data for columns of expression. \cr
#'     \code{log} \cr 
#'       \tab A data.frame logging statistics for each fit. \cr
#'   } 
#' @examples
#' set.seed(101)
#' exprs <- h0testr::sim1(n_obs=20, n_feats=30)$mat
#' exprs <- log2(exprs + 1)
#' feats <- data.frame(feature_id=rownames(exprs))
#' samps <- data.frame(observation_id=colnames(exprs))
#' state <- list(expression=exprs, features=feats, samples=samps)
#' config <- list()
#'
#' ## untransformed example:
#' out <- h0testr::impute_rf(state, config, is_log_transformed=FALSE, verbose=FALSE)
#' state2 <- out$state
#' summary(c(state$expression))    ## Note number of NAs
#' summary(c(state2$expression))   ## Note number of NAs
#' head(state$expression)
#' round(head(state2$expression))
#' print(out$log)
#' 
#' ## log-transformed example:
#' state$expression <- log2(state$expression + 1)
#' out <- h0testr::impute_rf(state, config, is_log_transformed=TRUE, verbose=FALSE)
#' state2 <- out$state
#' summary(c(state$expression))    ## Note number of NAs
#' summary(c(state2$expression))   ## Note number of NAs
#' head(state$expression)
#' round(head(state2$expression))
#' print(out$log)

impute_rf <- function(state, config, is_log_transformed=NULL, 
    f_imp=impute_sample_lod, ntree=100, mtry=NULL, aug_mult=0.33, 
    aug_add=0, aug_steps=NULL, verbose=T) {
  
  check_config(config)

  f.need_pkgs("randomForest", "impute_rf", config)
  
  is_log_transformed <- f.is_log_transformed(is_log_transformed, config,
    "impute_rf")

  config$is_log_transformed <- is_log_transformed
  
  if(is.null(aug_steps)) aug_steps <- config$impute_aug_steps
  if(is.null(aug_steps)) aug_steps <- 3
  if(aug_steps < 0) {
    f.err("impute_rf: aug_steps < 0; aug_steps:", aug_steps, config=config)
  }
  
  if(!is.matrix(state$expression)) {
    f.err("impute_rf: !is.matrix(state$expression)", "\n",
      "class(state$expression):", class(state$expression), config=config)
  }
  
  ## NA is only indicator of missing value:

  n_miss <- apply(state$expression, 1, function(v) sum(is.na(v)))
  o <- order(n_miss, decreasing=F)
  tbl <- NULL
  
  for(idx_feat in o) {
    
    if(n_miss[idx_feat] < 1) next
    f.msg("processing", rownames(state$expression)[idx_feat], config=config)
    
    x_train <- f_imp(state, config)$expression
    if(any(is.na(c(x_train)))) {
      f.err("impute_rf: f_imp returned NAs", config=config)
    }
    if(is_log_transformed) x_train <- 2^x_train - 1
    x_train <- f.augment_affine(x_train, mult=aug_mult, add=aug_add, steps=aug_steps)
    if(is_log_transformed) x_train <- log2(x_train + 1)
    
    y <- state$expression[idx_feat, , drop=T]
    i_miss <- is.na(y)
    col_names_miss <- colnames(state$expression)[i_miss]
    i_miss_train <- colnames(x_train) %in% col_names_miss
    x_train_i <- x_train[-idx_feat, !i_miss_train, drop=F]
    y_train_i <- x_train[idx_feat, !i_miss_train, drop=T]
    if(is.null(mtry)) { mtry0 <- round(sqrt(nrow(x_train_i))) } else { mtry0 <- mtry }
    
    fit <- randomForest::randomForest(x=t(x_train_i), y=y_train_i, mtry=mtry0, 
      ntree=ntree, importance=F, replace=T)
    
    x_miss <- x_train[-idx_feat, which(i_miss), drop=F]
    y_miss <- stats::predict(fit, newdata=t(x_miss), type="response")
    
    if(any(y_miss < 0)) {
      f.msg("impute_rf: y_miss < 0; deferred to f_imp; y_miss:", y_miss, config=config)
      y_miss[y_miss < 0] <- NA
    }
    state$expression[idx_feat, i_miss] <- y_miss
    
    tm_stmp <- format(Sys.time(), format='%Y%m%d%H%M%S')
    tbl_i <- data.frame(
      time=tm_stmp, feat=rownames(state$expression)[idx_feat], 
      mtry=mtry0, n_miss1=sum(i_miss), n_miss2=sum(is.na(y_miss)), 
      mean0=mean(y[!i_miss]), mean1=mean(y_miss), 
      pvar0=fit$rsq[1], pvar1=fit$rsq[ntree], 
      mse0=fit$mse[1], mse1=fit$mse[ntree]
    )
    
    if(verbose) f.log_obj(tbl_i, config=config)
    tbl <- rbind(tbl, tbl_i)
  }
  ## fall-back:
  if(any(is.na(state$expression))) {
    f.msg("impute_rf: fall-back imputation for", 
      sum(is.na(state$expression)), "features", config=config)
    state <- f_imp(state, config)
  }
  
  return(list(state=state, log=tbl))
}

#' Impute missing values using \code{glmnet}
#' @description
#'   Impute missing values using \code{glmnet} package.
#' @details 
#'   Imputes missing values using a \code{glmnet::cv.glmnet} model trained 
#'     using observations in which the feature was expressed. Iterates through 
#'     features beginning with those having fewest missing values. If 
#'     \code{aug_steps > 0}, augments observations with affine transormed versions.
#'     This is meant to enable extrapolation outside of observed intensity 
#'     range. If \code{aug_steps > 0}, assumes expression data have been previously
#'     \code{log(x+1)} transformed. If you want \code{0} to be considered missing, 
#'     and have \code{0} in the data, do something like 
#'     \code{exprs[exprs \%in\% 0] <- NA} prior to imputing. Augmentation uses 
#'     \code{f.aug_mult()}. 
#' @param state A list with elements like that returned by \code{read_data()}:
#'   \tabular{ll}{
#'     \code{expression} \cr \tab Numeric matrix with non-negative expression values. \cr
#'     \code{features}   \cr \tab A data.frame with feature meta-data for rows of expression. \cr
#'     \code{samples}    \cr \tab A data.frame with observation meta-data for columns of expression. \cr
#'   }
#' It is assumed that state$expression has been previously \code{log2(x+1)} transformed.
#' @param config List with configuration values. Does not use any keys so can pass empty list.
#' @param is_log_transformed Logical scalar: whether \code{state$expression} has
#'   been log transformed. Defaults to \code{config$is_log_transformed}, which
#'   \code{h0testr::init_state()} and \code{h0testr::normalize()} maintain;
#'   passing both is an error unless they agree.
#' @param f_imp Function to use for initial rough imputation.
#' @param alpha Numeric (between 0 and 1) number of trees in random forest.
#'  Default: \code{1.0}.
#' @param nfolds Numeric (greater than or equal to 2) number of folds for 
#'   cross-validation tuning of \code{lambda} value for \code{glmnet}.
#' @param measure Character in set \code{c("mae", "mse")}. Loss function used 
#'   during cross-validation tuning of \code{lambda} value for \code{glmnet}.
#' @param aug_mult Numeric affine shift for training data augmentation. 
#' @param aug_add Numeric affine shift for training data augmentation.
#' @param aug_steps Numeric (non-negative) number of augmentation steps. Set to
#'   \code{0} to skip augmentation. Default: \code{3}.
#' @param verbose Logical if \code{TRUE}, emits progress messages.
#' @return A list with the following elements:
#'   \tabular{llll}{
#'     \code{state} \cr 
#'       \tab A list with the following elements: \cr
#'       \tab \code{expression} \cr 
#'       \tab \tab \tab Numeric matrix with non-negative expression values. \cr
#'       \tab \code{features}   \cr 
#'       \tab \tab \tab A data.frame with feature meta-data for rows of expression. \cr
#'       \tab \code{samples}    \cr 
#'       \tab \tab \tab A data.frame with observation meta-data for columns of expression. \cr
#'     \code{log} \cr 
#'       \tab A data.frame logging statistics for each fit. \cr
#'   } 
#' @examples
#' set.seed(101)
#' exprs <- h0testr::sim1(n_obs=20, n_feats=30)$mat
#' feats <- data.frame(feature_id=rownames(exprs))
#' samps <- data.frame(observation_id=colnames(exprs))
#' state <- list(expression=exprs, features=feats, samples=samps)
#' config <- list()
#'
#' ## example with untransformed data:
#' out <- h0testr::impute_glmnet(state, config, is_log_transformed=FALSE, verbose=FALSE)
#' state2 <- out$state
#' summary(c(state$expression))    ## Note number of NAs
#' summary(c(state2$expression))   ## Note number of NAs
#' head(state$expression)
#' round(head(state2$expression))
#' print(out$log)
#' 
#' ## example with log-transformed data:
#' state$expression <- log2(state$expression + 1)
#' out <- h0testr::impute_glmnet(state, config, is_log_transformed=TRUE, verbose=FALSE)
#' state2 <- out$state
#' summary(c(state$expression))    ## Note number of NAs
#' summary(c(state2$expression))   ## Note number of NAs
#' head(state$expression)
#' round(head(state2$expression))
#' print(out$log)

impute_glmnet <- function(state, config, is_log_transformed=NULL,
    f_imp=impute_unif_sample_lod, nfolds=5, alpha=NULL, measure="mae", 
    aug_mult=0.33, aug_add=0, aug_steps=NULL, verbose=T) {
  
  check_config(config)

  f.need_pkgs("glmnet", "impute_glmnet", config)
  
  is_log_transformed <- f.is_log_transformed(is_log_transformed, config,
    "impute_glmnet")

  ## f_imp is another imputer:

  config$is_log_transformed <- is_log_transformed
  
  if(is.null(alpha)) alpha <- config$impute_alpha
  if(is.null(alpha)) alpha <- 1.0
  if(alpha < 0 || alpha > 1) {
    f.err("impute_glmnet: alpha < 0 || alpha > 1; alpha:", 
      alpha, config=config)
  }
  
  if(is.null(aug_steps)) aug_steps <- config$impute_aug_steps
  if(is.null(aug_steps)) aug_steps <- 3
  if(aug_steps < 0) {
    f.err("impute_glmnet: aug_steps < 0; aug_steps:", aug_steps, config=config)
  }
  
  if(!is.matrix(state$expression)) {
    f.err("impute_glmnet: !is.matrix(state$expression)", "\n",
      "class(state$expression):", class(state$expression), config=config)
  }
  
  ## NA is only indicator of a missing value:

  n_miss <- apply(state$expression, 1, function(v) sum(is.na(v)))
  o <- order(n_miss, decreasing=F)
  tbl <- NULL
  
  for(idx_feat in o) {
    
    if(n_miss[idx_feat] < 1) next
    f.msg("processing", rownames(state$expression)[idx_feat], config=config)
    
    x_train <- f_imp(state, config)$expression
    if(is_log_transformed) x_train <- 2^x_train - 1
    x_train <- f.augment_affine(x_train, mult=aug_mult, add=aug_add, steps=aug_steps)
    if(is_log_transformed) x_train <- log2(x_train + 1)   
    
    y <- state$expression[idx_feat, , drop=T]
    i_miss <- is.na(y)
    col_names_miss <- colnames(state$expression)[i_miss]
    i_miss_train <- colnames(x_train) %in% col_names_miss
    x_train_i <- x_train[-idx_feat, !i_miss_train, drop=F]
    y_train_i <- x_train[idx_feat, !i_miss_train, drop=T]
    
    fit <- glmnet::cv.glmnet(x=t(x_train_i), y=y_train_i, family="gaussian", 
      alpha=alpha, type.measure=measure, nfolds=nfolds, parallel=F)
    
    x_miss <- x_train[-idx_feat, which(i_miss), drop=F]
    y_miss <- stats::predict(fit, newx=t(x_miss), s=fit$lambda.1se, type="response")
    if(any(y_miss < 0)) {
      f.msg("impute_glmnet: y_miss < 0; deferred to f_imp; y_miss:", y_miss, config=config)
      y_miss[y_miss < 0] <- NA
    }
    state$expression[idx_feat, i_miss] <- y_miss
    
    idx_lambda <- which(fit$lambda %in% fit$lambda.1se)[1] 
    tm_stmp <- format(Sys.time(), format='%Y%m%d%H%M%S')
    
    tbl_i <- data.frame(
      time=tm_stmp, feat=rownames(state$expression)[idx_feat], 
      alpha=alpha, nfolds=nfolds, 
      n_miss1=sum(i_miss), n_miss2=sum(is.na(y_miss)), 
      mean0=mean(y[!i_miss]), mean1=mean(y_miss, na.rm=T), 
      cvm0=max(fit$cvm), cvm1=fit$cvm[idx_lambda], 
      cvup0=max(fit$cvup), cvup1=fit$cvup[idx_lambda]
    )
    if(verbose) f.log_obj(tbl_i, config=config)
    tbl <- rbind(tbl, tbl_i)
  }
  ## fall-back:
  if(any(is.na(state$expression))) {
    f.msg("impute_glmnet: fall-back imputation for", 
      sum(is.na(state$expression)), "features", config=config)
    state <- f_imp(state, config)
  }
  
  return(list(state=state, log=tbl))
}

#' Impute missing values with \code{impute::impute.knn()}
#' @description
#'   Impute missing values with \code{impute::impute.knn()}.
#' @details
#'   Blocks larger than \code{maxp} are recursively divided into smaller 
#'     sub-blocks prior to imputation.
#' @param state A list with elements like that returned by \code{read_data()}:
#'   \tabular{ll}{
#'     \code{expression} \cr \tab Numeric matrix with non-negative expression values. \cr
#'     \code{features}   \cr \tab A data.frame with feature meta-data for rows of expression. \cr
#'     \code{samples}    \cr \tab A data.frame with observation meta-data for columns of expression. \cr
#'   } 
#' @param config List with configuration values. Does not use any keys so can pass empty list.
#' @param k Number of nearest neighbor features to average, \code{k >= 1}.
#'   Defaults to \code{config$impute_k}, and to \code{round(sqrt(n_features))}
#'   when that is unset. \code{round(sqrt(n_features))} is also a cap: a larger
#'   value, from either source, is reduced to it and reported, rather than being
#'   refused or passed through.
#' @param rowmax Maximum proportion missing per row, else use row means; \code{0 < rowmax <= 1}.
#' @param colmax Maximum proportion missing per column, else use column means; \code{0 < colmax <= 1}.
#' @param maxp Scalar max number of genes per imputation block; \code{0 < maxp <= n_features}.
#' @return 
#'   An updated \code{state} list with the following elements:
#'   \tabular{ll}{
#'     \code{expression} \cr \tab Numeric matrix with non-negative imputed 
#'       expression values. \cr
#'     \code{features}   \cr \tab Feature meta-data \code{data.frame} 
#'       corresponding to rows of \code{expression}. \cr
#'     \code{samples}    \cr \tab Observation meta-data \code{data.frame} 
#'       corresponding to columns of \code{expression}. \cr
#'   } 
#' @examples
#' ## setup state and config, including prefiltering. No effects are planted, these
#' ##   examples being about missing values rather than about testing:
#' set.seed(101)
#' samps <- h0testr::sim_samples(factors=list(grp=c("ctl", "trt"), sex=c("F", "M")),
#'   n_per_cell=3)
#' sim <- h0testr::sim_design(samps, frm=~grp + sex, test_term="grp", n_genes=100)
#' state <- sim$state
#' config <- sim$config    ## frm, the id columns, and is_log_transformed=FALSE (raw data)
#' rm(samps, sim)
#' state <- h0testr::filter_features(state, config, n_samples_min=3)
#' state <- h0testr::filter_observations(state, config, n_features_min=30)
#'
#' ## impute:
#' state2 <- h0testr::impute_knn(state, config)
#' summary(c(state$expression))    ## Note number of NAs
#' summary(c(state2$expression))   ## Note number of NAs
#' head(state$expression)
#' round(head(state2$expression))

impute_knn <- function(state, config, k=NULL, rowmax=0.5, colmax=0.8, maxp=1500) {
  
  check_config(config)

  f.need_pkgs("impute", "impute_knn", config)
  
  if(!is.matrix(state$expression)) {
    f.err("impute_knn: !is.matrix(state$expression)", "\n",
      "class(state$expression):", class(state$expression), config=config)
  }
  
  k_max <- round(sqrt(nrow(state$expression)))
  if(is.null(k)) k <- config$impute_k
  if(is.null(k)) k <- k_max
  
  if(k > k_max) {
    f.msg("impute_knn: k > k_max; k:", k, "\n", 
      "setting k to k_max:", k_max, config=config)
    k <- k_max
  }
  
  out <- impute::impute.knn(state$expression, k=k, 
    rowmax=rowmax, colmax=colmax, maxp=maxp)
    
  mat <- out$data
  mat <- f.pos_mat(mat, config, fn_name="impute_knn")
  state$expression <- mat

  return(state)
}

#' Impute missing values with \code{imputeLCMD::impute.MinDet()}
#' @description
#'   Impute missing values with \code{imputeLCMD::impute.MinDet()}.
#' @details
#'   Only \code{NAs} considered missing. If you want \code{0} to be considered 
#'     missing, do something like \code{exprs[exprs \%in\% 0] <- NA} prior to 
#'     imputing. 
#' @param state A list with elements like that returned by \code{read_data()}:
#'   \tabular{ll}{
#'     \code{expression} \cr \tab Numeric matrix with non-negative expression values. \cr
#'     \code{features}   \cr \tab A data.frame with feature meta-data for rows of expression. \cr
#'     \code{samples}    \cr \tab A data.frame with observation meta-data for columns of expression. \cr
#'   } 
#' @param config List with configuration values. Does not require any keys so 
#'   can pass empty list.
#' @param impute_quantile Quantile for imputation; scalar with 
#'   \code{0 <= q < 1.0}. Default: \code{0.01}.
#' @return 
#'   An updated \code{state} list with the following elements:
#'   \tabular{ll}{
#'     \code{expression} \cr \tab Numeric matrix with non-negative imputed 
#'       expression values. \cr
#'     \code{features}   \cr \tab Feature meta-data \code{data.frame} 
#'       corresponding to rows of \code{expression}. \cr
#'     \code{samples}    \cr \tab Observation meta-data \code{data.frame} 
#'       corresponding to columns of \code{expression}. \cr
#'   } 
#' @examples
#' ## setup state and config, including prefiltering. No effects are planted, these
#' ##   examples being about missing values rather than about testing:
#' set.seed(101)
#' samps <- h0testr::sim_samples(factors=list(grp=c("ctl", "trt"), sex=c("F", "M")),
#'   n_per_cell=3)
#' sim <- h0testr::sim_design(samps, frm=~grp + sex, test_term="grp", n_genes=100)
#' state <- sim$state
#' config <- sim$config    ## frm, the id columns, and is_log_transformed=FALSE (raw data)
#' rm(samps, sim)
#' state <- h0testr::filter_features(state, config, n_samples_min=3)
#' state <- h0testr::filter_observations(state, config, n_features_min=30)
#'
#' ## impute:
#' state2 <- h0testr::impute_min_det(state, config)
#' summary(c(state$expression))    ## Note number of NAs
#' summary(c(state2$expression))   ## Note number of NAs
#' head(state$expression)
#' round(head(state2$expression))

impute_min_det <- function(state, config, impute_quantile=NULL) {
  
  check_config(config)

  f.need_pkgs("imputeLCMD", "impute_min_det", config)
  
  if(!is.matrix(state$expression)) {
    f.err("impute_min_det: !is.matrix(state$expression)", "\n",
      "class(state$expression):", class(state$expression), config=config)
  }
  
  if(is.null(impute_quantile)) impute_quantile <- config$impute_quantile
  if(is.null(impute_quantile)) impute_quantile <- 0.01
  
  mat <- imputeLCMD::impute.MinDet(state$expression, q=impute_quantile)
  mat <- f.pos_mat(mat, config, fn_name="impute_min_det")
  state$expression <- mat
  
  return(state)
}

#' Impute missing values with \code{imputeLCMD::impute.MinProb()}
#' @description
#'   Impute missing values with \code{imputeLCMD::impute.MinProb()}.
#' @details
#'   Only \code{NAs} considered missing. If you want \code{0} to be considered 
#'     missing, do something like \code{exprs[exprs \%in\% 0] <- NA} prior to 
#'     imputing. 
#' @param state A list with elements like that returned by \code{read_data()}:
#'   \tabular{ll}{
#'     \code{expression} \cr \tab Numeric matrix with non-negative expression values. \cr
#'     \code{features}   \cr \tab A data.frame with feature meta-data for rows of expression. \cr
#'     \code{samples}    \cr \tab A data.frame with observation meta-data for columns of expression. \cr
#'   } 
#' @param config List with configuration values. Does not use any keys so 
#'   can pass empty list.
#' @param is_log_transformed Logical scalar: whether \code{state$expression} has
#'   been log transformed. Defaults to \code{config$is_log_transformed}, which
#'   \code{h0testr::init_state()} and \code{h0testr::normalize()} maintain;
#'   passing both is an error unless they agree.
#' @param impute_quantile Quantile to use for imputation; scalar with 
#'   \code{0 <= impute_quantile < 1.0}. Default: \code{0.01}.
#' @param scale. Scale parameter for normal distribution; scalar 
#'   with \code{0 < scale.}. Default: \code{1.0}.
#' @return 
#'   An updated \code{state} list with the following elements:
#'   \tabular{ll}{
#'     \code{expression} \cr \tab Numeric matrix with non-negative imputed 
#'       expression values. \cr
#'     \code{features}   \cr \tab Feature meta-data \code{data.frame} 
#'       corresponding to rows of \code{expression}. \cr
#'     \code{samples}    \cr \tab Observation meta-data \code{data.frame} 
#'       corresponding to columns of \code{expression}. \cr
#'   } 
#' @examples
#' ## setup state and config, including prefiltering. No effects are planted, these
#' ##   examples being about missing values rather than about testing:
#' set.seed(101)
#' samps <- h0testr::sim_samples(factors=list(grp=c("ctl", "trt"), sex=c("F", "M")),
#'   n_per_cell=3)
#' sim <- h0testr::sim_design(samps, frm=~grp + sex, test_term="grp", n_genes=100)
#' state <- sim$state
#' config <- sim$config    ## frm, the id columns, and is_log_transformed=FALSE (raw data)
#' rm(samps, sim)
#' state <- h0testr::filter_features(state, config, n_samples_min=3)
#' state <- h0testr::filter_observations(state, config, n_features_min=30)
#'
#' ## untransformed example:
#' state2 <- h0testr::impute_min_prob(state, config, is_log_transformed=FALSE)
#' summary(c(state$expression))    ## Note number of NAs
#' summary(c(state2$expression))   ## Note number of NAs
#' head(state$expression)
#' round(head(state2$expression))
#' 
#' ## log-transformed example:
#' state$expression <- log2(state$expression + 1)
#' config$is_log_transformed <- TRUE   ## normalize() would; argument and config must agree
#' state2 <- h0testr::impute_min_prob(state, config, is_log_transformed=TRUE)
#' summary(c(state$expression))    ## Note number of NAs
#' summary(c(state2$expression))   ## Note number of NAs
#' head(state$expression)
#' round(head(state2$expression))

impute_min_prob <- function(state, config, is_log_transformed=NULL, 
    impute_quantile=NULL, scale.=NULL) {
    
  check_config(config)

  f.need_pkgs("imputeLCMD", "impute_min_prob", config)
  
  is_log_transformed <- f.is_log_transformed(is_log_transformed, config,
    "impute_min_prob")
  
  if(!is.matrix(state$expression)) {
    f.err("impute_min_prob: !is.matrix(state$expression)", "\n",
      "class(state$expression):", class(state$expression), config=config)
  }
  
  if(is.null(impute_quantile)) impute_quantile <- config$impute_quantile
  if(is.null(impute_quantile)) impute_quantile <- 0.01
  if(is.null(scale.)) scale. <- config$impute_scale
  if(is.null(scale.)) scale. <- 1.0
  
  mat <- state$expression
  if(!is_log_transformed) mat <- log2(mat + 1)  ## otherwise can get negative
  
  mat <- imputeLCMD::impute.MinProb(mat, q=impute_quantile, tune.sigma=scale.)

  if(!is_log_transformed) mat <- (2^mat) - 1
  mat <- f.pos_mat(mat, config, is_log_transformed, "impute_min_prob")
  
  state$expression <- mat
  
  return(state)
}

#' Impute missing values with \code{imputeLCMD::impute.QRILC()}
#' @description
#'   Impute missing values with \code{imputeLCMD::impute.QRILC()}.
#' @details
#'   Only \code{NAs} considered missing. If you want \code{0} to be considered 
#'     missing, do something like \code{exprs[exprs \%in\% 0] <- NA} prior to 
#'     imputing. 
#' @param state A list with elements like that returned by \code{read_data()}:
#'   \tabular{ll}{
#'     \code{expression} \cr \tab Numeric matrix with non-negative expression values. \cr
#'     \code{features}   \cr \tab A data.frame with feature meta-data for rows of expression. \cr
#'     \code{samples}    \cr \tab A data.frame with observation meta-data for columns of expression. \cr
#'   } 
#' @param config List with configuration values. Does not use any keys so can 
#'   pass empty list.
#' @param is_log_transformed Logical scalar: whether \code{state$expression} has
#'   been log transformed. Defaults to \code{config$is_log_transformed}, which
#'   \code{h0testr::init_state()} and \code{h0testr::normalize()} maintain;
#'   passing both is an error unless they agree.
#' @param scale. Scaling parameter for normal distribution; numeric scalar 
#'   with \code{0 < scale.}. Default: \code{1.0}.
#' @return 
#'   An updated \code{state} list with the following elements:
#'   \tabular{ll}{
#'     \code{expression} \cr \tab Numeric matrix with non-negative imputed 
#'       expression values. \cr
#'     \code{features}   \cr \tab Feature meta-data \code{data.frame} 
#'       corresponding to rows of \code{expression}. \cr
#'     \code{samples}    \cr \tab Observation meta-data \code{data.frame} 
#'       corresponding to columns of \code{expression}. \cr
#'   } 
#' @examples
#' ## setup state and config, including prefiltering. No effects are planted, these
#' ##   examples being about missing values rather than about testing:
#' set.seed(101)
#' samps <- h0testr::sim_samples(factors=list(grp=c("ctl", "trt"), sex=c("F", "M")),
#'   n_per_cell=3)
#' sim <- h0testr::sim_design(samps, frm=~grp + sex, test_term="grp", n_genes=100)
#' state <- sim$state
#' config <- sim$config    ## frm, the id columns, and is_log_transformed=FALSE (raw data)
#' rm(samps, sim)
#' state <- h0testr::filter_features(state, config, n_samples_min=3)
#' state <- h0testr::filter_observations(state, config, n_features_min=30)
#' 
#' ## untransformed example. QRILC draws from a truncated distribution fitted in
#' ##   log space, and a draw below zero there back-transforms to a raw value
#' ##   below zero, which cannot be a measurement, so wrap in try(). 
#' state2 <- try(h0testr::impute_qrilc(state, config, is_log_transformed=FALSE))
#' summary(c(state$expression))    ## Note number of NAs
#' if(!inherits(state2, "try-error")) summary(c(state2$expression))
#' 
#' ## log-transformed example:
#' state$expression <- log2(state$expression + 1)
#' config$is_log_transformed <- TRUE   ## normalize() would; argument and config must agree
#' state2 <- h0testr::impute_qrilc(state, config, is_log_transformed=TRUE)
#' summary(c(state$expression))    ## Note number of NAs
#' summary(c(state2$expression))   ## Note number of NAs
#' head(state$expression)
#' round(head(state2$expression))

impute_qrilc <- function(state, config, is_log_transformed=NULL, scale.=NULL) {
  
  check_config(config)

  f.need_pkgs("imputeLCMD", "impute_qrilc", config)
  
  is_log_transformed <- f.is_log_transformed(is_log_transformed, config,
    "impute_qrilc")
  
  if(!is.matrix(state$expression)) {
    f.err("impute_qrilc: !is.matrix(state$expression)", "\n",
      "class(state$expression):", class(state$expression), config=config)
  }
  
  if(is.null(scale.)) scale. <- config$impute_scale
  if(is.null(scale.)) scale. <- 1.0
  
  mat <- state$expression
  if(!is_log_transformed) mat <- log2(mat + 1)  ## otherwise can get negative
  
  obj <- imputeLCMD::impute.QRILC(mat, tune.sigma=scale.)
  
  mat <- obj[[1]]
  if(!is_log_transformed) mat <- (2^mat) - 1
  mat <- f.pos_mat(mat, config, is_log_transformed, "impute_qrilc")
  state$expression <- mat
  
  return(state)
}

#' Impute missing values with \code{pcaMethods::pca()}
#' @description
#'   Impute missing values with \code{pcaMethods::pca()}.
#' @details
#'   Only \code{NAs} considered missing. If you want \code{0} to be considered 
#'     missing, do something like \code{exprs[exprs \%in\% 0] <- NA} prior to 
#'     imputing. 
#'   Options for parameter \code{method} are described more fully in 
#'     the \code{pcaMethods} package documentation. Acceptable values include:
#'   \tabular{ll}{
#'     \code{bpca}      \cr \tab Bayesian PCA (https://doi.org/10.1093/bioinformatics/btg287). \cr
#'     \code{ppca}      \cr \tab Probabilistic PCA (https://dl.acm.org/doi/10.5555/3008904.3008993). \cr
#'     \code{svdImpute} \cr \tab Use svdImpute (https://doi.org/10.1093/bioinformatics/17.6.520). \cr
#'   }
#' @param state A list with elements like that returned by \code{read_data()}:
#'   \tabular{ll}{
#'     \code{expression} \cr \tab Numeric matrix with non-negative expression values. \cr
#'     \code{features}   \cr \tab A data.frame with feature meta-data for rows of expression. \cr
#'     \code{samples}    \cr \tab A data.frame with observation meta-data for columns of expression. \cr
#'   } 
#' @param config List with configuration values. Does not use any keys so can pass empty list.
#' @param is_log_transformed Logical scalar: whether \code{state$expression} has
#'   been log transformed. Defaults to \code{config$is_log_transformed}, which
#'   \code{h0testr::init_state()} and \code{h0testr::normalize()} maintain;
#'   passing both is an error unless they agree.
#' @param n_pcs Number (scalar numeric >= 1) of principle components to compute. Default: \code{5}.
#' @param method Method to use. Scalar character in \code{c("bpca", "ppca", "svdImpute")}.
#' @return 
#'   An updated \code{state} list with the following elements:
#'   \tabular{ll}{
#'     \code{expression} \cr \tab Numeric matrix with non-negative imputed 
#'       expression values. \cr
#'     \code{features}   \cr \tab Feature meta-data \code{data.frame} 
#'       corresponding to rows of \code{expression}. \cr
#'     \code{samples}    \cr \tab Observation meta-data \code{data.frame} 
#'       corresponding to columns of \code{expression}. \cr
#'   } 
#' @examples
#' ## setup state and config, including prefiltering:
#' set.seed(101)
#' samps <- h0testr::sim_samples(factors=list(grp=c("ctl", "trt"), sex=c("F", "M")),
#'   n_per_cell=3)
#' sim <- h0testr::sim_design(samps, frm=~grp + sex, test_term="grp", n_genes=100)
#' state <- sim$state
#' config <- sim$config    ## frm, the id columns, and is_log_transformed=FALSE (raw data)
#' rm(samps, sim)
#' state <- h0testr::filter_features(state, config, n_samples_min=3)
#' state <- h0testr::filter_observations(state, config, n_features_min=30)
#'
#' summary(c(state$expression))    ## Note number of NAs
#' head(state$expression)
#'
#' ## impute using bayesian pca:
#' state2 <- h0testr::impute_pca(state, config, method="bpca", 
#'   is_log_transformed=FALSE)
#' summary(c(state2$expression))   ## Note number of NAs
#' round(head(state2$expression))
#'
#' ## impute using probabilistic pca:
#' state2 <- h0testr::impute_pca(state, config, method="ppca", 
#'   is_log_transformed=FALSE)
#' summary(c(state2$expression))   ## Note number of NAs
#' round(head(state2$expression))
#'
#' ## raw example impute as linear combo of n_pcs eigengenes,
#' ##   which can yield below 0 measurements so wrap in try():
#' try(h0testr::impute_pca(state, config, method="svdImpute",
#'   is_log_transformed=FALSE))
#'
#' ## the same method on the log scale, where the reconstruction has no floor to
#' ##   fall below and every imputed value is an ordinary log intensity:
#' state_log <- state
#' state_log$expression <- log2(state_log$expression + 1)
#' config$is_log_transformed <- TRUE   ## normalize() would; argument and config must agree
#' state2 <- h0testr::impute_pca(state_log, config, method="svdImpute",
#'   is_log_transformed=TRUE)
#' summary(c(state2$expression))   ## Note number of NAs
#' round(head(state2$expression), 2)

impute_pca <- function(state, config, is_log_transformed=NULL,
    n_pcs=NULL, method="bpca") {
  
  check_config(config)

  f.need_pkgs("pcaMethods", "impute_pca", config)
  
  is_log_transformed <- f.is_log_transformed(is_log_transformed, config,
    "impute_pca")
  
  if(!is.matrix(state$expression)) {
    f.err("impute_pca: !is.matrix(state$expression)", "\n",
      "class(state$expression):", class(state$expression), config=config)
  }
  
  allowed <- c("bpca", "ppca", "svdImpute")
  if(!(method %in% allowed)) {
    f.err("impute_pca: !(method %in% allowed); method:", method, config=config)
  }
  
  if(is.null(n_pcs)) n_pcs <- config$impute_npcs
  n_pcs_max <- round(sqrt(nrow(state$expression)))
  if(is.null(n_pcs)) n_pcs <- n_pcs_max
  if(n_pcs > n_pcs_max) {
    f.msg(
      "impute_pca: n_pcs > n_pcs_max; n_pcs:", n_pcs, "\n",
      "setting n_pcs to n_pcs_max:", n_pcs_max, config=config
    )
    n_pcs <- n_pcs_max
  }
  
  mat <- state$expression
  if(!is_log_transformed) mat <- log2(mat + 1)  ## otherwise can get negative
  
  ## wants sample rows and 'variables' columns:
  obj <- pcaMethods::pca(t(mat), nPcs=n_pcs, method=method)
  mat <- pcaMethods::completeObs(obj)

  if(!is_log_transformed) mat <- (2^mat) - 1
  mat <- t(mat)
  mat <- f.pos_mat(mat, config, is_log_transformed, "impute_pca")
  state$expression <- mat
  
  return(state)
}

#' Impute missing values with \code{pcaMethods::llsImpute()}
#' @description
#'   Impute missing values with \code{pcaMethods::llsImpute()}.
#' @details
#'   Impute using local least squares.
#'   Only \code{NAs} considered missing. If you want \code{0} to be considered 
#'     missing, do something like \code{exprs[exprs \%in\% 0] <- NA} prior to 
#'     imputing. 
#' @param state A list with elements like that returned by \code{read_data()}:
#'   \tabular{ll}{
#'     \code{expression} \cr \tab Numeric matrix with non-negative expression values. \cr
#'     \code{features}   \cr \tab A data.frame with feature meta-data for rows of expression. \cr
#'     \code{samples}    \cr \tab A data.frame with observation meta-data for columns of expression. \cr
#'   } 
#' @param config List with configuration values. Does not use any keys so can 
#'   pass empty list.
#' @param is_log_transformed Logical scalar: whether \code{state$expression} has
#'   been log transformed. Defaults to \code{config$is_log_transformed}, which
#'   \code{h0testr::init_state()} and \code{h0testr::normalize()} maintain;
#'   passing both is an error unless they agree.
#' @param method Name of correlation method; scalar character in 
#'   \code{c("pearson", "kendall", "spearman")}
#' @param k Number (scalar numeric with \code{k >= 2}) of features per local cluster. Default: \code{5}.
#' @param maxit Maximum number of iterations. Scalar numeric with \code{maxit >= 2}.
#' @return 
#'   An updated \code{state} list with the following elements:
#'   \tabular{ll}{
#'     \code{expression} \cr \tab Numeric matrix with non-negative imputed 
#'       expression values. \cr
#'     \code{features}   \cr \tab Feature meta-data \code{data.frame} 
#'       corresponding to rows of \code{expression}. \cr
#'     \code{samples}    \cr \tab Observation meta-data \code{data.frame} 
#'       corresponding to columns of \code{expression}. \cr
#'   } 
#' @examples
#' ## setup state and config, including prefiltering. No effects are planted, these
#' ##   examples being about missing values rather than about testing:
#' set.seed(101)
#' samps <- h0testr::sim_samples(factors=list(grp=c("ctl", "trt"), sex=c("F", "M")),
#'   n_per_cell=3)
#' sim <- h0testr::sim_design(samps, frm=~grp + sex, test_term="grp", n_genes=100)
#' state <- sim$state
#' config <- sim$config    ## frm, the id columns, and is_log_transformed=FALSE (raw data)
#' rm(samps, sim)
#' state <- h0testr::filter_features(state, config, n_samples_min=3)
#' state <- h0testr::filter_observations(state, config, n_features_min=30)
#'
#' ## example with untransformed data. LLS regresses a feature on its correlated
#' ##   neighbors, and a reconstruction can land below zero, so this is wrapped
#' ##   in try().
#' summary(c(state$expression))    ## Note number of NAs
#' state2 <- try(h0testr::impute_lls(state, config, is_log_transformed=FALSE))
#' if(!inherits(state2, "try-error")) summary(c(state2$expression))
#'
#' ## example with log-transformed data, where the reconstruction has no floor to
#' ##   fall below and every imputed value is an ordinary log intensity:
#' state$expression <- log2(state$expression + 1)
#' config$is_log_transformed <- TRUE   ## normalize() would; argument and config must agree
#' state2 <- h0testr::impute_lls(state, config, is_log_transformed=TRUE)
#' summary(c(state$expression))    ## Note number of NAs
#' summary(c(state2$expression))   ## Note number of NAs
#' head(state$expression)
#' round(head(state2$expression))

impute_lls <- function(state, config, is_log_transformed=NULL, 
    k=NULL, method="pearson", maxit=100) {
  
  check_config(config)

  f.need_pkgs("pcaMethods", "impute_lls", config)
  
  is_log_transformed <- f.is_log_transformed(is_log_transformed, config,
    "impute_lls")
  
  if(is.null(k)) k <- config$impute_k
  if(is.null(k)) k <- 5
  
  if(!is.matrix(state$expression)) {
    f.err("impute_lls: !is.matrix(state$expression)", "\n",
      "class(state$expression):", class(state$expression), config=config)
  }
  
  mat <- state$expression
  ## NA is the only indicator of missing value:

  if(sum(is.na(c(mat))) %in% 0) {
    f.msg("no missing values found; returning input state", config=config)
    return(state)
  }
  if(!is_log_transformed) mat <- log2(mat + 1)  ## otherwise can get negative
  
  ## wants sample rows and 'variables' as columns:
  obj <- pcaMethods::llsImpute(t(mat), k=k, center=F, completeObs=T,
    correlation=method, allVariables=F, maxSteps=maxit, xval=NULL, verbose=F)

  mat <- pcaMethods::completeObs(obj)
  if(!is_log_transformed) mat <- (2^mat) - 1
  mat <- t(mat)
  mat <- f.pos_mat(mat, config, is_log_transformed, "impute_lls")
  state$expression <- mat
  
  return(state)
}

#' Impute missing values with \code{missForest::missForest()}
#' @description
#'   Impute missing values with \code{missForest::missForest()}.
#' @details
#'   Only \code{NAs} considered missing. If you want \code{0} to be considered 
#'     missing, do something like \code{exprs[exprs \%in\% 0] <- NA} prior to 
#'     imputing. 
#' @param state A list with elements like that returned by \code{read_data()}:
#'   \tabular{ll}{
#'     \code{expression} \cr \tab Numeric matrix with non-negative expression values. \cr
#'     \code{features}   \cr \tab A data.frame with feature meta-data for rows of expression. \cr
#'     \code{samples}    \cr \tab A data.frame with observation meta-data for columns of expression. \cr
#'   } 
#' @param config List with configuration values. Does not use any keys so can pass empty list.
#' @param maxit Maximum number of iterations used during fitting.
#' @param ntree Number of trees to grow in forest.
#' @return 
#'   An updated \code{state} list with the following elements:
#'   \tabular{ll}{
#'     \code{expression} \cr \tab Numeric matrix with non-negative imputed 
#'       expression values. \cr
#'     \code{features}   \cr \tab Feature meta-data \code{data.frame} 
#'       corresponding to rows of \code{expression}. \cr
#'     \code{samples}    \cr \tab Observation meta-data \code{data.frame} 
#'       corresponding to columns of \code{expression}. \cr
#'   } 
#' @examples
#' ## setup state and config, including prefiltering. No effects are planted, these
#' ##   examples being about missing values rather than about testing:
#' set.seed(101)
#' samps <- h0testr::sim_samples(factors=list(grp=c("ctl", "trt"), sex=c("F", "M")),
#'   n_per_cell=3)
#' sim <- h0testr::sim_design(samps, frm=~grp + sex, test_term="grp", n_genes=100)
#' state <- sim$state
#' config <- sim$config    ## frm, the id columns, and is_log_transformed=FALSE (raw data)
#' rm(samps, sim)
#' state <- h0testr::filter_features(state, config, n_samples_min=3)
#' state <- h0testr::filter_observations(state, config, n_features_min=30)
#'
#' ## impute:
#' state2 <- h0testr::impute_missforest(state, config)
#' summary(state$expression)     ## Note number of NAs
#' summary(state2$expression)    ## Note number of NAs
#' head(state$expression)
#' round(head(state2$expression))

impute_missforest <- function(state, config, maxit=10, ntree=100) {

  f.need_pkgs("missForest", "impute_missforest", config)
    
  if(!is.matrix(state$expression)) {
    f.err("impute_missforest: !is.matrix(state$expression)", "\n",
      "class(state$expression):", class(state$expression), config=config)
  }
  
  obj <- missForest::missForest(t(state$expression), maxiter=maxit, ntree=ntree)
  
  mat <- t(obj$ximp)
  mat <- f.pos_mat(mat, config, fn_name="impute_missforest")
  state$expression <- mat

  return(state)
}

#' Get choices for \code{method} parameter of \code{h0testr::impute()}
#' @description 
#'   Returns valid values for \code{method} parameter of \code{h0testr::impute()}
#' @details
#'   Corresponding methods have prefix \code{impute_}.
#' @return
#'   Character vector with options for \code{method} parameter of \code{h0testr::impute()}.
#' @examples
#' impute_methods <- h0testr::impute_methods()
#' cat("available methods:\n")
#' for(method in impute_methods) {
#'   cat("method:", method, "\n")
#' }

impute_methods <- function() {
  return(
    c("sample_lod", "unif_global_lod", "unif_sample_lod", "qrilc", "bpca", 
      "ppca", "svdImpute", "missforest", "knn", "lls", "min_det", "min_prob", 
      "glm_binom", "loess_logit", "rf", "glmnet", "rnorm_feature", "none"
    )
  )
}

#' Impute missing values
#' @description
#'   Impute missing values in \code{state$expression} according to settings
#'     in \code{config}.
#' @details
#'   Imputes missing values using a method and parameters specified in
#'     \code{config}. The scale of the data is not assumed: it comes from the
#'     \code{is_log_transformed} argument or from
#'     \code{config$is_log_transformed}, which \code{h0testr::init_state()} sets
#'     and \code{h0testr::normalize()} updates, and whichever answers is recorded
#'     in the returned \code{config} and passed to the method, so that methods
#'     without an \code{is_log_transformed} argument of their own still see it.
#'   Only \code{NA} is missing. If you want \code{0} to be considered missing,
#'     and have \code{0} in the data, do something like
#'     \code{exprs[exprs \%in\% 0] <- NA} prior to imputing;
#'     \code{h0testr::init_state()} does this for raw input.
#'   An imputed value that cannot be a measurement is an error.
#'   See documentation for \code{h0testr::new_config()} 
#'     for more detailed description of configuration parameters. 
#' @param state A list with elements like that returned by \code{read_data()}:
#'   \tabular{ll}{
#'     \code{expression} \cr \tab Numeric matrix with non-negative expression values. \cr
#'     \code{features}   \cr \tab A data.frame with feature meta-data for rows of expression. \cr
#'     \code{samples}    \cr \tab A data.frame with observation meta-data for columns of expression. \cr
#'   }
#' @param config List with configuration values. Required keys are 
#'   \code{feat_col} and \code{obs_col}. Recognizes the following keys:
#'   \tabular{ll}{
#'     \code{feat_col}          \cr \tab Name of column in \code{state$features} matching \code{rownames(state$expression)}. \cr
#'     \code{obs_col}           \cr \tab Name of column in \code{state$samples} matching \code{colnames(state$expression)}. \cr
#'     \code{impute_method}     \cr \tab Name of a method in list returned by \code{h0testr::impute_methods()}. \cr
#'     \code{impute_quantile}   \cr \tab Used if parameter \code{impute_quantile} is unset. \cr
#'     \code{impute_scale}      \cr \tab Used if parameter \code{scale.} is unset. \cr
#'     \code{impute_span}       \cr \tab Used if parameter \code{span} is unset. \cr
#'     \code{impute_k}          \cr \tab Used if parameter \code{k} is unset. \cr
#'     \code{impute_npcs}       \cr \tab Used if parameter \code{n_pcs} is unset. \cr
#'     \code{impute_aug_steps}  \cr \tab Used if parameter \code{aug_steps} is unset. \cr
#'     \code{impute_n_pts}      \cr \tab Used if parameter \code{n_pts} is unset. \cr
#'   }
#' @param method Method to use (required). A character scalar from the list 
#'   returned by \code{h0testr::impute_methods()}.
#' @param is_log_transformed Logical scalar: whether \code{state$expression} has
#'   been log transformed. Defaults to \code{config$is_log_transformed}, which
#'   \code{h0testr::init_state()} and \code{h0testr::normalize()} maintain;
#'   passing both is an error unless they agree. It is an error for both to be
#'   unset.
#' @param k Number of nearest neighbors passed to methods for 
#'   \code{c("knn", "lls")}.
#' @param span Span passed to method for \code{"loess_logit"}.
#' @param n_pcs Number of PCs passed to methods for 
#'   \code{c("bpca", "ppca", "svdImpute")}.
#' @param impute_quantile Passed to methods for 
#'   \code{c("unif_global_lod", "unif_sample_lod", "min_det", "min_prob")}.
#' @param scale. Scale parameter passed to methods for 
#'   \code{c("rnorm_feature", "qrilc")}.
#' @param aug_steps Number of augmentations passed to methods for 
#'   \code{c("glmnet", "rf")}.
#' @param alpha Mixing parameter passed to method for \code{"glmnet"}.
#' @param n_pts Points in prediction grid; passed to methods for 
#'   \code{c("glm_binom", "loess_logit")}.
#' @param verbose Logical scalar passed to methods for 
#'   \code{c("glmnet", "rf")}.
#' @return Updated \code{state} list with the following elements:
#'   \tabular{ll}{
#'     \code{expression} \cr \tab Numeric matrix with non-negative expression values. \cr
#'     \code{features}   \cr \tab A data.frame with feature meta-data for rows of expression. \cr
#'     \code{samples}    \cr \tab A data.frame with observation meta-data for columns of expression. \cr
#'   } 
#' @examples
#' set.seed(101)
#' exprs <- h0testr::sim1(n_obs=8, n_feats=12)$mat
#' exprs <- log2(exprs + 1)
#' feats <- data.frame(feature_id=rownames(exprs))
#' samps <- data.frame(observation_id=colnames(exprs))
#' state <- list(expression=exprs, features=feats, samples=samps)
#' config <- list(feat_col="feature_id", obs_col="observation_id")
#'
#' cat("available methods:\n")
#' print(h0testr::impute_methods())
#'
#' ## impute using method passed as parameter:
#' out <- h0testr::impute(state, config, method="unif_sample_lod", 
#'   is_log_transformed=FALSE)
#' summary(c(state$expression))        ## note number of NAs
#' summary(c(out$state$expression))    ## note number of NAs
#' head(state$expression)
#' round(head(out$state$expression))
#'
#' ## impute using method passed in configuration:
#' config$impute_method <- "unif_sample_lod"
#' out <- h0testr::impute(state, config, is_log_transformed=FALSE)
#' summary(c(state$expression))        ## note number of NAs
#' summary(c(out$state$expression))    ## note number of NAs
#' head(state$expression)
#' round(head(out$state$expression))

impute <- function(state, config, method=NULL, is_log_transformed=NULL, 
    k=NULL, span=NULL, n_pcs=NULL, impute_quantile=NULL, scale.=NULL, 
    aug_steps=NULL, alpha=NULL, n_pts=NULL, verbose=NULL) {

  check_config(config)
  f.check_state(state, config)
  
  if(is.null(method)) method <- config$impute_method
  if(is.null(method)) {
    f.err("impute: both method and config$impute_method are unset", 
      config=config)
  }
  if(!(method %in% impute_methods())) {
    f.err("impute: !(method %in% impute_methods()); method:", 
      method, config=config)
  }
  
  is_log_transformed <- f.is_log_transformed(is_log_transformed, config, "impute")

  config$is_log_transformed <- is_log_transformed

  ## corresponding methods should have reasonable defaults for NULLs in config:
  if(is.null(k)) k <- config$impute_k
  if(is.null(span)) span <- config$impute_span
  if(is.null(n_pcs)) n_pcs <- config$impute_npcs
  if(is.null(impute_quantile)) impute_quantile <- config$impute_quantile
  if(is.null(scale.)) scale. <- config$impute_scale
  if(is.null(aug_steps)) aug_steps <- config$impute_aug_steps
  if(is.null(alpha)) alpha <- config$impute_alpha
  if(is.null(n_pts)) n_pts <- config$impute_n_pts
  if(is.null(verbose)) verbose <- config$verbose
  if(is.null(verbose)) verbose <- TRUE
  
  f.msg("impute: method:", method, 
    "; is_log_transformed:", is_log_transformed,
    "; k:", k, "; span:", span, "; n_pcs:", n_pcs, 
    "impute_quantile:", impute_quantile, "; scale.:", scale.,
    "; aug_steps:", aug_steps, "; alpha:", alpha, "; n_pts:", n_pts,
    config=config)

  if(method %in% "unif_global_lod") {
    state <- impute_unif_global_lod(state, config, 
      impute_quantile=impute_quantile)
  } else if(method %in% "sample_lod") {
    state <- impute_sample_lod(state, config)
  } else if(method %in% "unif_sample_lod") {
    state <- impute_unif_sample_lod(state, config, 
      impute_quantile=impute_quantile)
  } else if(method %in% "rnorm_feature") {
    state <- impute_rnorm_feature(state, config,
      is_log_transformed=is_log_transformed, scale.=scale.)
  } else if(method %in% "glm_binom") {
    state <- impute_glm_binom(state, config, 
      is_log_transformed=is_log_transformed, n_pts=n_pts)
  } else if(method %in% "loess_logit") {
    state <- impute_loess_logit(state, config, span=span, n_pts=n_pts)
  } else if(method %in% "glmnet") {
    out <- impute_glmnet(state, config, 
      is_log_transformed=is_log_transformed, 
      alpha=alpha, aug_steps=aug_steps, verbose=verbose)
    state <- out$state
  } else if(method %in% "rf") {
    out <- impute_rf(state, config, is_log_transformed=is_log_transformed, 
      aug_steps=aug_steps, verbose=verbose)
    state <- out$state
  } else if(method %in% "missforest") {
    state <- impute_missforest(state, config)
  } else if(method %in% "knn") {
    state <- impute_knn(state, config, k=k)
  } else if(method %in% "lls") {
    state <- impute_lls(state, config, 
      is_log_transformed=is_log_transformed, k=k)
  } else if(method %in% "bpca") {
    state <- impute_pca(state, config, is_log_transformed=is_log_transformed, 
      n_pcs=n_pcs, method="bpca")
  } else if(method %in% "ppca") {
    state <- impute_pca(state, config, is_log_transformed=is_log_transformed, 
      n_pcs=n_pcs, method="ppca")
  } else if(method %in% "svdImpute") {
    state <- impute_pca(state, config, is_log_transformed=is_log_transformed, 
      n_pcs=n_pcs, method="svdImpute")
  } else if(method %in% "min_det") {
    state <- impute_min_det(state, config, impute_quantile=impute_quantile)
  } else if(method %in% "min_prob") {
    state <- impute_min_prob(state, config, 
      is_log_transformed=is_log_transformed, 
      impute_quantile=impute_quantile, scale.=scale.)
  } else if(method %in% "qrilc") {
    state <- impute_qrilc(state, config, 
      is_log_transformed=is_log_transformed, scale.=scale.)
  } else if(method %in% "none") {
    f.msg("skipping imputation: config$impute_method %in% 'none'", 
      config=config)
  } else {
    f.err("impute: unexpected method:", method, config=config)
  }
    
  f.check_state(state, config)
  f.report_state(state, config)
  
  prfx <- "imputed"
  if(!is.null(config$run_order)) {
    i <- config$run_order %in% "impute"
    if(any(i)) {
      prfx <- paste0(which(i)[1] + 2, ".imputed")
    }
  }
  f.save_state(state, config, prefix=prfx)
  
  return(list(state=state, config=config))
}
