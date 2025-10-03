## helper to ensure positive and finite imputed values:

f.pos_mat <- function(mat, config, maxit=10) {
  
  n <- 0
  i <- mat <= 0
  i[is.na(i)] <- F
  
  ## heuristic to ensure all strictly positive; 
  ##   iteratively nudge higher till positive:
  while(any(c(i))) {
    n <- n + 1
    if(n > maxit) break
    mat[i] <- log2(2^mat[i] + 1)
    i <- mat <= 0
    i[is.na(i)] <- F
  }
  
  if(any(c(mat) <= 0)) {
    f.err("f.pos_mat: min(c(mat)):", min(c(mat)), "after", n-1, 
      "iterations.", config=config)
  }
  
  ## heuristic to ensure all finite; randomly assign between 
  ##   max(finite_vals) and 2*max(finite_vals):
  i <- is.infinite(mat)
  i[is.na(i)] <- F
  maxval <- max(mat[!i], na.rm=T)
  mat[i] <- maxval * (1 + stats::runif(sum(i)))
  
  return(mat)
}

#' Impute missing values between 0 and global LOD
#' @description
#'   Impute missing values by randomly drawing from uniform distribution 
#'     below an estimated global LOD. 
#' @details Imputed values are random draws from uniform distribution over the 
#'   interval \code{[0, LOD]}, where \code{LOD} is an estimate global limit of 
#'   detection. Only \code{NA} values are considered as missing, so if you want 
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
#'     \code{impute_quantile}  \cr \tab Quantile of signal distribution to use as estimated limit of detection (LOD).
#'   }
#' @param impute_quantile Numeric between 0 and 1 specifying 
#'   minimum non-zero/\code{NA} expression value for each feature to use as global
#'   LOD (limit of detection). If \code{NULL}, \code{config$impute_quantile} used.
#'   Default: \code{0}.
#' @return An updated `state` list with the following elements:
#'   \tabular{ll}{
#'     \code{expression} \cr \tab Numeric matrix with non-negative expression values. \cr
#'     \code{features}   \cr \tab A data.frame with feature meta-data for rows of expression. \cr
#'     \code{samples}    \cr \tab A data.frame with observation meta-data for columns of expression. \cr
#'   } 
#' Returned \code{state$xpression} matrix contains strictly positive numeric values.
#' @examples
#' set.seed(101)
#' exprs <- h0testr::sim1(n_obs=6, n_feats=12)$mat
#' exprs <- log2(exprs + 1)
#' feats <- data.frame(feature_id=rownames(exprs))
#' samps <- data.frame(observation_id=colnames(exprs))
#' state <- list(expression=exprs, features=feats, samples=samps)
#' config <- list(impute_quantile=0.05)
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
  
  v <- apply(state$expression, 1, min, na.rm=T)
  v <- v[!is.na(v)]
  v <- v[v > 0]        ## minimum non-zero values for each gene (with non-0 vals)
  
  max_val <- stats::quantile(v, probs=impute_quantile[1], na.rm=T)    ## quantile of min vals
  if(is.na(max_val)) {
    f.err("impute_unif_global_lod: is.na(max_val)", config=config)
  }
  names(max_val) <- NULL
  f.msg("impute_unif_global_lod: inpute_quantile: ", impute_quantile, 
    "; max_val:", max_val, config=config)
  
  f <- function(v, max_val) {
    i <- is.na(v)
    if(any(i)) {
      v[i] <- stats::runif(sum(i), 0, max_val)
    }
    v
  }
  mat <- t(apply(state$expression, 1, f, max_val))
  mat <- f.pos_mat(mat, config)
  state$expression <- mat
  
  return(state)
}

#' Impute missing values between 0 and sample LOD
#' @description
#'   Impute missing values by randomly drawing from uniform distribution below 
#'     an estimated observation-specific limit of detection (LOD). 
#' @details Imputed values are random draws from uniform distribution over the 
#'   interval \code{[0, LOD]}, where LOD is an estimate observation limit of 
#'   detection. Only \code{NA} values are considered as missing, so if you 
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
#'     \code{impute_quantile}  \cr \tab Quantile of signal distribution to use as estimated limit of detection (LOD).
#'   }
#' @param impute_quantile Numeric between 0 and 1 specifying 
#'   minimum non-zero/\code{NA} expression value in each observation to use as the
#'   observation-specific LOD (limit of detection). Default: \code{0}.
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
#' config <- list(impute_quantile=0.05)
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
  
  f <- function(v) {
    i <- is.na(v)
    if(all(i)) {
      f.err("impute_unif_sample_lod: all(is.na(state$expression[, column]))", 
        config=config)
    }
    if(any(i)) {
      max_val <- stats::quantile(v, probs=impute_quantile, na.rm=T)
      v[i] <- stats::runif(sum(i), 0, max_val)
    }
    v
  }
  mat <- apply(state$expression, 2, f)
  mat <- f.pos_mat(mat, config)
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
#' @param state A list with elements like that returned by \code{read_data()}:
#'   \tabular{ll}{
#'     \code{expression} \cr \tab Numeric matrix with non-negative expression values. \cr
#'     \code{features}   \cr \tab A data.frame with feature meta-data for rows of expression. \cr
#'     \code{samples}    \cr \tab A data.frame with observation meta-data for columns of expression. \cr
#'   } 
#' @param config List with configuration values. No keys used, so can send empty list.
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
#' config <- list()
#' state2 <- h0testr::impute_sample_lod(state, config)
#' summary(c(state$expression))   ## note number of NAs
#' summary(c(state2$expression))  ## note number of NAs
#' head(state$expression)
#' round(head(state2$expression))

impute_sample_lod <- function(state, config) {

  if(!is.matrix(state$expression)) {
    f.err("impute_sample_lod: !is.matrix(state$expression)", "\n",
      "class(state$expression):", class(state$expression), config=config)
  }
  
  f <- function(v) {
    v[is.na(v)] <- min(v, na.rm=T)
    v
  }
  
  state$expression <- apply(state$expression, 2, f)
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
#'   prior to imputing. All values are gauranteed non-negative.
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
#'   }
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
#' config <- list()    
#' state2 <- h0testr::impute_rnorm_feature(state, config)
#' summary(c(state2$expression))    ## note number of NAs
#' round(head(state2$expression))
#'
#' ## impute using passed scale. parameter:
#' config <- list()    
#' state2 <- h0testr::impute_rnorm_feature(state, config, scale.=2)
#' summary(c(state2$expression))    ## note number of NAs
#' round(head(state2$expression))
#'
#' ## impute using impute_scale parameter from config:
#' config <- list(impute_scale=0.5)    
#' state2 <- h0testr::impute_rnorm_feature(state, config)
#' summary(c(state2$expression))    ## note number of NAs
#' round(head(state2$expression))

impute_rnorm_feature <- function(state, config, scale.=NULL) {

  check_config(config)

  if(!is.matrix(state$expression)) {
    f.err("impute_rnorm_feature: !is.matrix(state$expression)", "\n",
      "class(state$expression):", class(state$expression), config=config)
  } 
  i <- c(state$expression) < 0
  i[is.na(i)] <- F
  if(any(i)) {
    f.err("impute_rnorm_feature: state$expression contains negative values", 
      config=config)
  }
  if(is.null(scale.)) scale. <- config$impute_scale
  if(is.null(scale.)) scale. <- 1.0
  
  f <- function(v) {
    i <- is.na(v)
    if(all(i)) {
      f.err("impute_rnorm_feature: all(is.na(state$expression[row, ]))",
        config=config)
    }
    if(any(i)) {
      m <- mean(v[!i])
      s_min <- sqrt(m)
      if(sum(!i) >= 2) {
        s <- stats::sd(v[!i]) * scale.
      } else {
        s <- s_min
      }
      if(s < s_min) s <- s_min

      v[i] <- stats::rnorm(sum(i), mean=m, sd=s)
      i <- v < 0
      i[is.na(i)] <- F
      while(any(i)) {
        v[i] <- stats::rnorm(sum(i), mean=m, sd=s)
        i <- v < 0
        i[is.na(i)] <- F
      }
    }
    return(v)
  }
  state$expression <- t(apply(state$expression, 1, f))
  
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
#' @param is_log_transformed Logical scalar: if \code{state$expression} has 
#'   been log transformed.
#' @param n_pts Numeric greater than one. Granularity of prediction 
#'   grid. Larger values lead to less chance of duplicate imputed values.
#'   Larger values require more compute time and memory. Default: \code{1e7}.
#' @param off Numeric offset for calculating 
#'   \code{p.missing = (n.missing + off) / (n.total + off)}.
#' @param f_mid Function to use for calculating central tendency of 
#'   feature expression across samples. 
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
    n_pts=NULL, off=1, f_mid=stats::median) {
  
  check_config(config)
  
  if(!is.matrix(state$expression)) {
    f.err("impute_glm_binom: !is.matrix(state$expression)", "\n",
      "class(state$expression):", class(state$expression), config=config)
  }
  
  if(!is.logical(is_log_transformed)) {
    f.err("impute_glm_binom: !is.logical(is_log_transformed)", "\n",
      "is_log_transformed:", is_log_transformed, "\n",
      "typeof(is_log_transformed):", typeof(is_log_transformed), 
      config=config
    )
  }
  
  mat <- state$expression
  if(!is_log_transformed) mat <- log2(mat + 1) 
  
  if(is.null(n_pts)) n_pts <- config$impute_n_pts
  if(is.null(n_pts)) n_pts <- 1e7
  if(n_pts <= 0) {
    f.err("impute_glm_binom: n_pts <= 0; n_pts:", n_pts, config=config)
  }
  
  m <- apply(mat, 1, f_mid, na.rm=T)
  m[is.na(m)] <- 0
  n0 <- apply(mat, 1, function(v) sum(is.na(v) | v %in% 0))  ## n.missing
  n1 <- ncol(mat) - n0                                     
  p <- (n0 + off) / (ncol(mat) + off)
  dat <- data.frame(n0=n0, n1=n1, p=p, m=m)
  
  fit <- stats::glm(cbind(n0, n1) ~ m, data=dat, family="binomial")
  
  m_new <- seq(from=max(c(mat), na.rm=T) / n_pts, 
    to=max(c(mat), na.rm=T), length.out=n_pts)
  
  p_hat <- stats::predict(fit, newdata=data.frame(m=m_new), type="response")
  p_hat[is.na(p_hat)] <- min(p_hat, na.rm=T)
  i_na <- is.na(c(mat))
  mat[i_na] <- sample(m_new, sum(i_na), replace=T, prob=p_hat)
  
  if(!is_log_transformed) mat <- 2^mat
  mat <- f.pos_mat(mat, config)
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
    off=0.1, f_mid=stats::median, degree=1, fam="symmetric") {
  
  check_config(config)
  
  if(!is.matrix(state$expression)) {
    f.err("impute_loess_logit: !is.matrix(state$expression)", "\n",
      "class(state$expression):", class(state$expression), config=config)
  }
  if(is.null(span)) span <- config$impute_span
  if(is.null(span)) span <- 0.5
  if(is.null(n_pts)) n_pts <- config$impute_n_pts
  if(is.null(n_pts)) n_pts <- 1e7
  
  m <- apply(state$expression, 1, f_mid, na.rm=T)
  m[is.na(m)] <- 0
  if(!any(is.finite(m))) {
    f.err("impute_loess_logit: !any(is.finite(m))", "\n",
      "sort(m):", sort(m), config=config)
  }
  m[is.infinite(m)] <- max(m[is.finite(m)])
  
  f <- function(v) sum(is.na(v) | v %in% 0)
  n0 <- apply(state$expression, 1, f)                        ## n.missing
  if(all(n0 %in% 0)) return(state)                           ## cannot model all zeros
  
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
#' @param is_log_transformed Logical scalar: if \code{state$expression} has 
#'   been log transformed.
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
  
  if(!is.logical(is_log_transformed)) {
    f.err("impute_rf: !is.logical(is_log_transformed)", "\n",
      "is_log_transformed:", is_log_transformed, "\n",
      "typeof(is_log_transformed):", typeof(is_log_transformed), config=config)
  }
  
  if(is.null(aug_steps)) aug_steps <- config$impute_aug_steps
  if(is.null(aug_steps)) aug_steps <- 3
  if(aug_steps < 0) {
    f.err("impute_rf: aug_steps < 0; aug_steps:", aug_steps, config=config)
  }
  
  if(!is.matrix(state$expression)) {
    f.err("impute_rf: !is.matrix(state$expression)", "\n",
      "class(state$expression):", class(state$expression), config=config)
  }
  
  n_miss <- apply(state$expression, 1, function(v) sum(is.na(v) | v %in% 0))
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
    i_miss <- is.na(y) | y %in% 0
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
#' @param is_log_transformed Logical scalar: if \code{state$expression} has 
#'   been log transformed.
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
  
  if(!is.logical(is_log_transformed)) {
    f.err(
      "impute_glmnet: !is.logical(is_log_transformed)", "\n",
      "is_log_transformed:", is_log_transformed, "\n",
      "typeof(is_log_transformed):", typeof(is_log_transformed), 
      config=config
    )
  }
  
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
  
  n_miss <- apply(state$expression, 1, function(v) sum(is.na(v) | v %in% 0))
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
    i_miss <- is.na(y) | y %in% 0
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
#' @param k Number of nearest neighbors (features) to use, with \code{1 <= k < n_features}. Default: \code{10}.
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
#' ## setup state and config, including prefiltering:
#' set.seed(101)
#' nsamps <- 6
#' sim <- h0testr::sim2(n_samps1=nsamps, n_samps2=nsamps, n_genes=100, n_genes_signif=20, 
#'   fold_change=1, peps_per_gene=1, reps_per_sample=1)
#' exprs <- sim$mat
#' feats <- data.frame(gene=rownames(exprs))
#' samps <- data.frame(
#'   obs=colnames(exprs), 
#'   grp=factor(c(rep("ctl", nsamps), rep("trt", nsamps))),
#'   sex=factor(rep(c("M", "F"), nsamps))
#' )
#' state <- list(expression=exprs, features=feats, samples=samps)
#' rm(nsamps, sim, exprs, feats, samps)
#' config <- list(frm=~grp+sex)    ## need frm for filter_features with filter_by_formula=TRUE
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
  mat <- f.pos_mat(mat, config)
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
#' ## setup state and config, including prefiltering:
#' set.seed(101)
#' nsamps <- 6
#' sim <- h0testr::sim2(n_samps1=nsamps, n_samps2=nsamps, n_genes=100, n_genes_signif=20, 
#'   fold_change=1, peps_per_gene=1, reps_per_sample=1)
#' exprs <- sim$mat
#' feats <- data.frame(gene=rownames(exprs))
#' samps <- data.frame(
#'   obs=colnames(exprs), 
#'   grp=factor(c(rep("ctl", nsamps), rep("trt", nsamps))),
#'   sex=factor(rep(c("M", "F"), nsamps))
#' )
#' state <- list(expression=exprs, features=feats, samples=samps)
#' rm(nsamps, sim, exprs, feats, samps)
#' config <- list(frm=~grp+sex)    ## need frm for filter_features with filter_by_formula=TRUE
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
  
  if(!is.matrix(state$expression)) {
    f.err("impute_min_det: !is.matrix(state$expression)", "\n",
      "class(state$expression):", class(state$expression), config=config)
  }
  
  if(is.null(impute_quantile)) impute_quantile <- config$impute_quantile
  if(is.null(impute_quantile)) impute_quantile <- 0.01
  
  mat <- imputeLCMD::impute.MinDet(state$expression, q=impute_quantile)
  mat <- f.pos_mat(mat, config)
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
#' @param is_log_transformed Logical scalar: if \code{state$expression} has 
#'   been log transformed.
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
#' ## setup state and config, including prefiltering:
#' set.seed(101)
#' nsamps <- 6
#' sim <- h0testr::sim2(n_samps1=nsamps, n_samps2=nsamps, n_genes=100, 
#'   n_genes_signif=20, fold_change=1, peps_per_gene=1, reps_per_sample=1)
#' exprs <- sim$mat
#' feats <- data.frame(gene=rownames(exprs))
#' samps <- data.frame(
#'   obs=colnames(exprs), 
#'   grp=factor(c(rep("ctl", nsamps), rep("trt", nsamps))),
#'   sex=factor(rep(c("M", "F"), nsamps))
#' )
#' state <- list(expression=exprs, features=feats, samples=samps)
#' rm(nsamps, sim, exprs, feats, samps)
#' config <- list(frm=~grp+sex)    ## need frm for filter_features with filter_by_formula=TRUE
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
#' state2 <- h0testr::impute_min_prob(state, config, is_log_transformed=TRUE)
#' summary(c(state$expression))    ## Note number of NAs
#' summary(c(state2$expression))   ## Note number of NAs
#' head(state$expression)
#' round(head(state2$expression))

impute_min_prob <- function(state, config, is_log_transformed=NULL, 
    impute_quantile=NULL, scale.=NULL) {
    
  check_config(config)
  
  if(!is.logical(is_log_transformed)) {
    f.err("impute_min_prob: !is.logical(is_log_transformed)", config=config)
  }
  
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
  mat <- f.pos_mat(mat, config)
  
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
#' @param is_log_transformed Logical scalar: if \code{state$expression} has 
#'   been log transformed.
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
#' ## setup state and config, including prefiltering:
#' set.seed(101)
#' nsamps <- 6
#' sim <- h0testr::sim2(n_samps1=nsamps, n_samps2=nsamps, n_genes=100, 
#'   n_genes_signif=20, fold_change=1, peps_per_gene=1, reps_per_sample=1)
#' exprs <- sim$mat
#' feats <- data.frame(gene=rownames(exprs))
#' samps <- data.frame(
#'   obs=colnames(exprs), 
#'   grp=factor(c(rep("ctl", nsamps), rep("trt", nsamps))),
#'   sex=factor(rep(c("M", "F"), nsamps))
#' )
#' state <- list(expression=exprs, features=feats, samples=samps)
#' rm(nsamps, sim, exprs, feats, samps)
#' config <- list(frm=~grp+sex)    ## need frm for filter_features with filter_by_formula=TRUE
#' state <- h0testr::filter_features(state, config, n_samples_min=3)
#' state <- h0testr::filter_observations(state, config, n_features_min=30)
#' 
#' ## untransformed example:
#' state2 <- h0testr::impute_qrilc(state, config, is_log_transformed=FALSE)
#' summary(c(state$expression))    ## Note number of NAs
#' summary(c(state2$expression))   ## Note number of NAs
#' head(state$expression)
#' round(head(state2$expression))
#' 
#' ## log-transformed example:
#' state$expression <- log2(state$expression + 1)
#' state2 <- h0testr::impute_qrilc(state, config, is_log_transformed=TRUE)
#' summary(c(state$expression))    ## Note number of NAs
#' summary(c(state2$expression))   ## Note number of NAs
#' head(state$expression)
#' round(head(state2$expression))

impute_qrilc <- function(state, config, is_log_transformed=NULL, scale.=NULL) {
  
  check_config(config)
  
  if(!is.logical(is_log_transformed)) {
    f.err("impute_qrilc: !is.logical(is_log_transformed)", config=config)
  }
  
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
  mat <- f.pos_mat(mat, config)
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
#' @param is_log_transformed Logical scalar: if \code{state$expression} has 
#'   been log transformed.
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
#' nsamps <- 6
#' sim <- h0testr::sim2(n_samps1=nsamps, n_samps2=nsamps, n_genes=100, n_genes_signif=20, 
#'   fold_change=1, peps_per_gene=1, reps_per_sample=1)
#' exprs <- sim$mat
#' feats <- data.frame(gene=rownames(exprs))
#' samps <- data.frame(
#'   obs=colnames(exprs), 
#'   grp=factor(c(rep("ctl", nsamps), rep("trt", nsamps))),
#'   sex=factor(rep(c("M", "F"), nsamps))
#' )
#' state <- list(expression=exprs, features=feats, samples=samps)
#' rm(nsamps, sim, exprs, feats, samps)
#' config <- list(frm=~grp+sex)    ## need frm for filter_features with filter_by_formula=TRUE
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
#' ## impute as linear combo of \code{n_pcs} eigengenes:
#' state2 <- h0testr::impute_pca(state, config, method="svdImpute", 
#'   is_log_transformed=FALSE)
#' summary(c(state2$expression))   ## Note number of NAs
#' round(head(state2$expression))

impute_pca <- function(state, config, is_log_transformed=NULL,
    n_pcs=NULL, method="bpca") {
  
  check_config(config)
  
  if(!is.logical(is_log_transformed)) {
    f.err("impute_pca: !is.logical(is_log_transformed)", config=config)
  }
  
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
  
  ## wants sample rows and 'variables' as columns:
  obj <- pcaMethods::pca(t(mat), nPcs=n_pcs, method=method)
  mat <- pcaMethods::completeObs(obj) 
  
  if(!is_log_transformed) mat <- (2^mat)
  mat <- t(mat)
  mat <- f.pos_mat(mat, config)
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
#' @param is_log_transformed Logical scalar: if \code{state$expression} has 
#'   been log transformed.
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
#' ## setup state and config, including prefiltering:
#' set.seed(101)
#' nsamps <- 6
#' sim <- h0testr::sim2(n_samps1=nsamps, n_samps2=nsamps, n_genes=100, n_genes_signif=20, 
#'   fold_change=1, peps_per_gene=1, reps_per_sample=1)
#' exprs <- sim$mat
#' feats <- data.frame(gene=rownames(exprs))
#' samps <- data.frame(
#'   obs=colnames(exprs), 
#'   grp=factor(c(rep("ctl", nsamps), rep("trt", nsamps))),
#'   sex=factor(rep(c("M", "F"), nsamps))
#' )
#' state <- list(expression=exprs, features=feats, samples=samps)
#' rm(nsamps, sim, exprs, feats, samps)
#' config <- list(frm=~grp+sex)    ## need frm for filter_features with filter_by_formula=TRUE
#' state <- h0testr::filter_features(state, config, n_samples_min=3)
#' state <- h0testr::filter_observations(state, config, n_features_min=30)
#'
#' ## example with untransformed data:
#' state2 <- h0testr::impute_lls(state, config, is_log_transformed=FALSE)
#' summary(c(state$expression))    ## Note number of NAs
#' summary(c(state2$expression))   ## Note number of NAs
#' head(state$expression)
#' round(head(state2$expression))
#' 
#' ## example with log-transformed data:
#' state$expression <- log2(state$expression + 1)
#' state2 <- h0testr::impute_lls(state, config, is_log_transformed=TRUE)
#' summary(c(state$expression))    ## Note number of NAs
#' summary(c(state2$expression))   ## Note number of NAs
#' head(state$expression)
#' round(head(state2$expression))

impute_lls <- function(state, config, is_log_transformed=NULL, 
    k=NULL, method="pearson", maxit=100) {
  
  check_config(config)
  
  if(!is.logical(is_log_transformed)) {
    f.err("impute_lls: !is.logical(is_log_transformed)", config=config)
  }
  
  if(is.null(k)) k <- config$impute_k
  if(is.null(k)) k <- 5
  
  if(!is.matrix(state$expression)) {
    f.err("impute_lls: !is.matrix(state$expression)", "\n",
      "class(state$expression):", class(state$expression), config=config)
  }
  
  mat <- state$expression
  if(sum(is.na(c(mat)) | c(mat) %in% 0) %in% 0) {
    f.msg("no missing values found; returning input state", config=config)
    return(state)
  }
  if(!is_log_transformed) mat <- log2(mat + 1)  ## otherwise can get negative
  
  ## wants sample rows and 'variables' as columns:
  obj <- pcaMethods::llsImpute(t(mat), k=k, center=F, completeObs=T,
    correlation=method, allVariables=F, maxSteps=maxit, xval=NULL, verbose=F)
  
  mat <- pcaMethods::completeObs(obj) 
  if(!is_log_transformed) mat <- (2^mat)
  mat <- t(mat)
  mat <- f.pos_mat(mat, config)
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
#' ## setup state and config, including prefiltering:
#' set.seed(101)
#' nsamps <- 6
#' sim <- h0testr::sim2(n_samps1=nsamps, n_samps2=nsamps, n_genes=100, n_genes_signif=20, 
#'   fold_change=1, peps_per_gene=1, reps_per_sample=1)
#' exprs <- sim$mat
#' feats <- data.frame(gene=rownames(exprs))
#' samps <- data.frame(
#'   obs=colnames(exprs), 
#'   grp=factor(c(rep("ctl", nsamps), rep("trt", nsamps))),
#'   sex=factor(rep(c("M", "F"), nsamps))
#' )
#' state <- list(expression=exprs, features=feats, samples=samps)
#' rm(nsamps, sim, exprs, feats, samps)
#' config <- list(frm=~grp+sex)  ## need frm for filter_features with filter_by_formula=TRUE
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
    
  if(!is.matrix(state$expression)) {
    f.err("impute_missforest: !is.matrix(state$expression)", "\n",
      "class(state$expression):", class(state$expression), config=config)
  }
  
  obj <- missForest::missForest(t(state$expression), maxiter=maxit, ntree=ntree)
  
  mat <- t(obj$ximp)
  mat <- f.pos_mat(mat, config)
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
#'     \code{config}. Assumes expression data have been previously
#'     \code{log(x+1)} transformed. If you want \code{0} to be considered missing, 
#'     and have \code{0} in the data, do something like 
#'     \code{exprs[exprs \%in\% 0] <- 0} prior to imputing. See invidual 
#'     \code{impute_*} methods for more details.
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
#'     \code{impute_n_pcs}      \cr \tab Used if parameter \code{n_pcs} is unset. \cr
#'     \code{impute_aug_steps}  \cr \tab Used if parameter \code{aug_steps} is unset. \cr
#'     \code{impute_n_pts}      \cr \tab Used if parameter \code{n_pts} is unset. \cr
#'   }
#' @param method Method to use (required). A character scalar from the list 
#'   returned by \code{h0testr::impute_methods()}.
#' @param is_log_transformed Logical scalar (required): if 
#'   \code{state$expression} has been log transformed.
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
  
  if(!is.logical(is_log_transformed)) {
    f.err("impute: !is.logical(is_log_transformed)", config=config)
  }
  
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
    state <- impute_rnorm_feature(state, config, scale.=scale.)
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
  
  if(!is.null(config$run_order)) {
    i <- config$run_order %in% "impute"
    if(any(i)) {
      prfx <- paste0(which(i)[1] + 2, ".imputed")
    } else {
      prfx <- "imputed"
    }
  }
  f.save_state(state, config, prefix=prfx)
  
  return(list(state=state, config=config))
}
