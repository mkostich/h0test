#' Impute missing values from between 0 and global LOD.
#' @description
#'   Impute missing values by randomly drawing from uniform distribution 
#'     below an estimated global LOD. 
#' @details Imputed values are random draws from uniform distribution over the 
#'   interval \code{[0, LOD]}, where \code{LOD} is an estimate global limit of 
#'   detection. Only \code{NA} values are considered as missing, so if you want 
#'   \code{0} to be considered missing, and have \code{0} in the data, do 
#'   something like \code{state$expression[state$expression \%in\% 0] <- NA} 
#'   prior to imputing.
#' @param state A list with elements like that returned by \code{f.read_data()}:
#'   \tabular{ll}{
#'     \code{expression} \cr \tab Numeric matrix with non-negative expression values. \cr
#'     \code{features}   \cr \tab A data.frame with feature meta-data for rows of expression. \cr
#'     \code{samples}    \cr \tab A data.frame with observation meta-data for columns of expression. \cr
#'   } 
#' @param config List with configuration values. Uses the following keys:
#'   \tabular{ll}{
#'     \code{log_file}         \cr \tab Path to log file (character); \code{log_file=""} outputs to console.
#'     \code{impute_quantile}  \cr \tab Quantile of signal distribution to use as estimated limit of detection (LOD).
#'   }
#' @param impute_quantile Numeric between 0 and 1 specifying 
#'   minimum non-zero/\code{NA} expression value for each feature to use as global
#'   LOD (limit of detection). If \code{NULL}, \code{config$impute_quantile} used.
#' @return An updated `state` list with the following elements:
#'   \tabular{ll}{
#'     \code{expression} \cr \tab Numeric matrix with non-negative expression values. \cr
#'     \code{features}   \cr \tab A data.frame with feature meta-data for rows of expression. \cr
#'     \code{samples}    \cr \tab A data.frame with observation meta-data for columns of expression. \cr
#'   } 
#' Returned \code{state$xpression} matrix contains strictly positive numeric values.
#' @examples
#' set.seed(101)
#' exprs <- h0testr::f.sim1(n_obs=6, n_feats=12)$mat
#' exprs <- log2(exprs + 1)
#' feats <- data.frame(feature_id=rownames(exprs))
#' samps <- data.frame(observation_id=colnames(exprs))
#' state <- list(expression=exprs, features=feats, samples=samps)
#' config <- list(log_file="", impute_quantile=0.05)
#' state2 <- h0testr::f.impute_unif_global_lod(state, config)
#' summary(c(state$expression))   ## note number of NAs
#' summary(c(state2$expression))  ## note number of NAs

f.impute_unif_global_lod <- function(state, config, impute_quantile=NULL) {

  if(!is.matrix(state$expression)) {
     f.err("f.impute_unif_global_lod: !is.matrix(state$expression)", config=config)
  }
  if(is.null(impute_quantile)) impute_quantile <- config$impute_quantile
  if(is.null(impute_quantile)) {
    f.err("f.impute_unif_global_lod: impute_quantile unset", config=config)
  }
  v <- apply(state$expression, 1, min, na.rm=T)
  v <- v[!is.na(v)]
  v <- v[v > 0]        ## minimum non-zero values for each gene (with non-0 vals)
  
  max_val <- stats::quantile(v, probs=impute_quantile[1], na.rm=T)    ## quantile of min vals
  if(is.na(max_val)) f.err("f.impute_unif_global_lod: is.na(max_val)", config=config)
  names(max_val) <- NULL
  f.msg("f.impute_unif_global_lod: inpute_quantile: ", impute_quantile, 
    "; max_val:", max_val, config=config)
  
  f <- function(v, max_val) {
    i <- is.na(v)
    if(any(i)) {
      v[i] <- stats::runif(sum(i), 0, max_val)
    }
    v
  }
  state$expression <- t(apply(state$expression, 1, f, max_val))
  
  return(state)
}

#' Impute missing values from between 0 and sample LOD.
#' @description
#'   Impute missing values by randomly drawing from uniform distribution below 
#'     an estimated observation-specific limit of detection (LOD). 
#' @details Imputed values are random draws from uniform distribution over the 
#'   interval \code{[0, LOD]}, where LOD is an estimate observation limit of 
#'   detection. Only \code{NA} values are considered as missing, so if you 
#'   want \code{0} to be considered missing, and have \code{0} in the data, do 
#'   something like \code{exprs[exprs \%in\% 0] <- NA} prior to imputing.
#' @param state A list with elements like that returned by \code{f.read_data()}:
#'   \tabular{ll}{
#'     \code{expression} \cr \tab Numeric matrix with non-negative expression values. \cr
#'     \code{features}   \cr \tab A data.frame with feature meta-data for rows of expression. \cr
#'     \code{samples}    \cr \tab A data.frame with observation meta-data for columns of expression. \cr
#'   } 
#' @param config List with configuration values. Uses the following keys:
#'   \tabular{ll}{
#'     \code{log_file}         \cr \tab Path to log file (character); \code{log_file=""} outputs to console.
#'     \code{impute_quantile}  \cr \tab Quantile of signal distribution to use as estimated limit of detection (LOD).
#'   }
#' @param impute_quantile Numeric between 0 and 1 specifying 
#'   minimum non-zero/\code{NA} expression value in each observation to use as the
#'   observation-specific LOD (limit of detection).
#' @return An updated \code{state} list with the following elements:
#'   \tabular{ll}{
#'     \code{expression} \cr \tab Numeric matrix with non-negative expression values. \cr
#'     \code{features}   \cr \tab A data.frame with feature meta-data for rows of expression. \cr
#'     \code{samples}    \cr \tab A data.frame with observation meta-data for columns of expression. \cr
#'   } 
#' @examples
#' set.seed(101)
#' exprs <- h0testr::f.sim1(n_obs=6, n_feats=12)$mat
#' exprs <- log2(exprs + 1)
#' feats <- data.frame(feature_id=rownames(exprs))
#' samps <- data.frame(observation_id=colnames(exprs))
#' state <- list(expression=exprs, features=feats, samples=samps)
#' config <- list(log_file="", impute_quantile=0.05)
#' state2 <- h0testr::f.impute_unif_sample_lod(state, config)
#' summary(c(state$expression))   ## note number of NAs
#' summary(c(state2$expression))  ## note number of NAs

f.impute_unif_sample_lod <- function(state, config, impute_quantile=NULL) {

  if(!is.matrix(state$expression)) {
    f.err("f.impute_unif_sample_lod: !is.matrix(state$expression)", 
      config=config)
  }
  if(is.null(impute_quantile)) impute_quantile <- config$impute_quantile
  if(is.null(impute_quantile)) {
    f.err("f.impute_unif_sample_lod: impute_quantile unset", config=config)
  }
  
  f <- function(v) {
    i <- is.na(v)
    if(all(i)) {
      f.err("f.impute_unif_sample_lod: all(is.na(state$expression[, column]))", 
        config=config)
    }
    if(any(i)) {
      max_val <- stats::quantile(v, probs=impute_quantile, na.rm=T)
      v[i] <- stats::runif(sum(i), 0, max_val)
    }
    v
  }
  state$expression <- apply(state$expression, 2, f)
  
  return(state)
}

#' Impute missing values as the sample LOD.
#' @description
#'   Impute missing values as the minimum observed value in the corresponding 
#'     observation.
#' @details Missing values are set to an estimate observation limit of 
#'   detection (LOD). Only \code{NA} values are considered as missing, so if you 
#'   want \code{0} to be considered missing, and have \code{0} in the data, do 
#'   something like \code{exprs[exprs \%in\% 0] <- NA} prior to imputing.
#' @param state A list with elements like that returned by \code{f.read_data()}:
#'   \tabular{ll}{
#'     \code{expression} \cr \tab Numeric matrix with non-negative expression values. \cr
#'     \code{features}   \cr \tab A data.frame with feature meta-data for rows of expression. \cr
#'     \code{samples}    \cr \tab A data.frame with observation meta-data for columns of expression. \cr
#'   } 
#' @param config List with configuration values. Uses the following keys:
#'   \tabular{ll}{
#'     \code{log_file}         \cr \tab Path to log file (character); \code{log_file=""} outputs to console.
#'   }
#' @return An updated \code{state} list with the following elements:
#'   \tabular{ll}{
#'     \code{expression} \cr \tab Numeric matrix with non-negative expression values. \cr
#'     \code{features}   \cr \tab A data.frame with feature meta-data for rows of expression. \cr
#'     \code{samples}    \cr \tab A data.frame with observation meta-data for columns of expression. \cr
#'   } 
#' @examples
#' set.seed(101)
#' exprs <- h0testr::f.sim1(n_obs=6, n_feats=12)$mat
#' exprs <- log2(exprs + 1)
#' feats <- data.frame(feature_id=rownames(exprs))
#' samps <- data.frame(observation_id=colnames(exprs))
#' state <- list(expression=exprs, features=feats, samples=samps)
#' config <- list(log_file="")
#' state2 <- h0testr::f.impute_sample_lod(state, config)
#' summary(c(state$expression))   ## note number of NAs
#' summary(c(state2$expression))  ## note number of NAs

f.impute_sample_lod <- function(state, config) {

  if(!is.matrix(state$expression)) {
    f.err("f.impute_sample_lod: !is.matrix(state$expression)", 
      config=config)
  }
  
  f <- function(v) {
    v[is.na(v)] <- min(v, na.rm=T)
    v
  }
  
  state$expression <- apply(state$expression, 2, f)
  return(state)
}

#' Impute missing values near the feature mean.
#' @description
#'   Impute missing values normally distributed around the feature mean.
#' @details Missing values are draws from 
#'   \code{normal(mean=mean(exprs[feature, ]), sd=scale.*sd(exprs[feature, ]))}. 
#'   Only non-\code{NA} values are used in calculation of \code{mean} and 
#'   \code{sd}. If you want \code{0} to be considered missing, and have 
#'   \code{0} in the data, do something like \code{exprs[exprs \%in\% 0] <- NA} 
#'   prior to imputing. All values are gauranteed non-negative.
#' @param state A list with elements like that returned by \code{f.read_data()}:
#'   \tabular{ll}{
#'     \code{expression} \cr \tab Numeric matrix with non-negative expression values. \cr
#'     \code{features}   \cr \tab A data.frame with feature meta-data for rows of expression. \cr
#'     \code{samples}    \cr \tab A data.frame with observation meta-data for columns of expression. \cr
#'   } 
#' @param config List with configuration values. Uses the following keys:
#'   \tabular{ll}{
#'     \code{log_file}      \cr \tab Path to log file (character); \code{log_file=""} outputs to console.
#'     \code{impute_scale}  \cr \tab Factor (numeric) rescaling variance of normal distribution from which draws are made. See details.
#'   }
#' @param scale. Numeric greater than zero, linearly scaling the
#'   dispersion around the feature mean.
#' @return An updated \code{state} list with the following elements:
#'   \tabular{ll}{
#'     \code{expression} \cr \tab Numeric matrix with non-negative expression values. \cr
#'     \code{features}   \cr \tab A data.frame with feature meta-data for rows of expression. \cr
#'     \code{samples}    \cr \tab A data.frame with observation meta-data for columns of expression. \cr
#'   } 
#' @examples
#' set.seed(101)
#' exprs <- h0testr::f.sim1(n_obs=6, n_feats=12)$mat
#' exprs <- log2(exprs + 1)
#' feats <- data.frame(feature_id=rownames(exprs))
#' samps <- data.frame(observation_id=colnames(exprs))
#' state <- list(expression=exprs, features=feats, samples=samps)
#' config <- list(log_file="", impute_scale.=1)
#' state2 <- h0testr::f.impute_rnorm_feature(state, config)
#' summary(c(state$expression))   ## note number of NAs
#' summary(c(state2$expression))  ## note number of NAs

f.impute_rnorm_feature <- function(state, config, scale.=NULL) {

  if(!is.matrix(state$expression)) {
    f.err("f.impute_rnorm_feature: !is.matrix(state$expression)", 
      config=config)
  } 
  i <- c(state$expression) < 0
  i[is.na(i)] <- F
  if(any(i)) {
    f.err("f.impute_rnorm_feature: state$expression contains negative values", 
      config=config)
  }
  if(is.null(scale.)) scale. <- config$impute_scale
  if(is.null(scale.)) {
    f.err("f.impute_rnorm_feature: scale. and config$impute_scale unset", 
      config=config)
  }
  
  f <- function(v) {
    i <- is.na(v)
    if(all(i)) {
      f.err("f.impute_rnorm_feature: all(is.na(state$expression[row,]))")
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

#' Impute by drawing from `p(missing|intensity)` from binomial glm.
#' @description
#'   Impute missing values by randomly drawing from an estimated density 
#'     of \code{p(missing|intensity)}.
#' @details Impute random draws by drawing from the estimated density of 
#'   \code{p(missing|intensity)}. Estimates \code{p(missing|intensity)} using 
#'   a binomial glm with formula \code{cbind(n.missing, n.found) ~ intensity},
#'   where intensity is assumed to be previously log transformed. If you  
#'   want \code{0} to be considered missing, and have \code{0} in the data, do 
#'   something like \code{exprs[exprs \%in\% 0] <- NA} prior to imputing. 
#' @param state A list with elements like that returned by \code{f.read_data()}:
#'   \tabular{ll}{
#'     \code{expression} \cr \tab Numeric matrix with non-negative expression values. \cr
#'     \code{features}   \cr \tab A data.frame with feature meta-data for rows of expression. \cr
#'     \code{samples}    \cr \tab A data.frame with observation meta-data for columns of expression. \cr
#'   } 
#' It is assumed that state$expression has been previously \code{log2(x+1)} transformed.
#' @param config List with configuration values. Uses the following keys:
#'   \tabular{ll}{
#'     \code{log_file}      \cr \tab Path to log file (character); \code{log_file=""} outputs to console.
#'     \code{impute_n_pts}  \cr \tab Numeric greater than one. Determines granularity of imputation. Larger values lead to finer grain.
#'   }
#' @param n_pts Numeric greater than one. Granularity of prediction 
#'   grid. Larger values lead to less chance of duplicate imputed values.
#'   Larger values require more compute time and memory.
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
#' exprs <- h0testr::f.sim1(n_obs=6, n_feats=12)$mat
#' exprs <- log2(exprs + 1)
#' feats <- data.frame(feature_id=rownames(exprs))
#' samps <- data.frame(observation_id=colnames(exprs))
#' state <- list(expression=exprs, features=feats, samples=samps)
#' config <- list(log_file="", impute_n_pts=1e7)
#' state2 <- h0testr::f.impute_glm_binom(state, config)
#' summary(c(state$expression))   ## note number of NAs
#' summary(c(state2$expression))  ## note number of NAs

f.impute_glm_binom <- function(state, config, n_pts=NULL, off=1, 
    f_mid=stats::median) {

  if(!is.matrix(state$expression)) {
    f.err("f.impute_glm_binom: !is.matrix(state$expression)", config=config)
  }
  if(is.null(n_pts)) n_pts <- config$impute_n_pts
  if(is.null(n_pts)) f.err("f.impute_glm_binom: n_pts and config$impute_n_pts unset", 
    config=config)
  if(n_pts <= 0) f.err("f.impute_glm_binom: n_pts <= 0", config=config)

  m <- apply(state$expression, 1, f_mid, na.rm=T)
  m[is.na(m)] <- 0
  n0 <- apply(state$expression, 1, function(v) sum(is.na(v) | v %in% 0))  ## n.missing
  n1 <- ncol(state$expression) - n0                                     
  p <- (n0 + off) / (ncol(state$expression) + off)
  dat <- data.frame(n0=n0, n1=n1, p=p, m=m)

  fit <- stats::glm(cbind(n0, n1) ~ m, data=dat, family="binomial")
  m <- seq(from=max(c(state$expression), na.rm=T) / n_pts, 
    to=max(c(state$expression), na.rm=T), length.out=n_pts)
  p <- stats::predict(fit, newdata=data.frame(m=m), type="response")
  p[is.na(p)] <- min(p, na.rm=T)
  i_na <- is.na(c(state$expression))
  state$expression[i_na] <- sample(m, sum(i_na), replace=T, prob=p)

  return(state)
}

#' Impute by drawing from \code{p(missing|intensity)} estimated with loess.
#' @description
#'   Impute missing values by randomly drawing from an estimated density 
#'     of \code{p(missing|intensity)}.
#' @details Impute random draws by drawing from the estimated density of 
#'   \code{p(missing|intensity)}. Estimates \code{p(missing|intensity)} using 
#'   a loess fit with form \code{log(p_missing/p_found) ~ intensity}, where 
#'   intensity is assumed to be previously log transformed. If you want 
#'   \code{0} to be considered missing, and have \code{0} in the data, do 
#'   something like \code{exprs[exprs \%in\% 0] <- NA} prior to imputing. 
#' @param state A list with elements like that returned by \code{f.read_data()}:
#'   \tabular{ll}{
#'     \code{expression} \cr \tab Numeric matrix with non-negative expression values. \cr
#'     \code{features}   \cr \tab A data.frame with feature meta-data for rows of expression. \cr
#'     \code{samples}    \cr \tab A data.frame with observation meta-data for columns of expression. \cr
#'   } 
#' It is assumed that state$expression has been previously \code{log2(x+1)} transformed.
#' @param config List with configuration values. Uses the following keys:
#'   \tabular{ll}{
#'     \code{log_file}      \cr \tab Path to log file (character); \code{log_file=""} outputs to console.
#'     \code{impute_n_pts}  \cr \tab Numeric greater than one. Determines granularity of imputation. Larger values lead to finer grain.
#'     \code{impute_span}   \cr \tab Span (numeric between zero and 1) for \code{loess} fit.
#'   }
#' @param span. Span for loess fit. Numeric in the open interval \code{(0, 1)}.
#' @param n_pts Numeric greater than one. Granularity of prediction 
#'   grid. Larger values lead to less chance of duplicate imputed values.
#'   Larger values require more compute time and memory.
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
#' exprs <- h0testr::f.sim1(n_obs=6, n_feats=12)$mat
#' exprs <- log2(exprs + 1)
#' feats <- data.frame(feature_id=rownames(exprs))
#' samps <- data.frame(observation_id=colnames(exprs))
#' state <- list(expression=exprs, features=feats, samples=samps)
#' config <- list(log_file="", impute_n_pts=1e7, impute_span=0.4)
#' state2 <- h0testr::f.impute_loess_logit(state, config)
#' summary(c(state$expression))   ## note number of NAs
#' summary(c(state2$expression))  ## note number of NAs

f.impute_loess_logit <- function(state, config, span.=NULL, n_pts=NULL, 
    off=1, f_mid=stats::median, degree=1, fam="symmetric") {

  if(!is.matrix(state$expression)) {
    f.err("f.impute_loess_logit: !is.matrix(state$expression)", config=config)
  }
  if(is.null(span.)) span. <- config$impute_span
  if(is.null(span.)) {
    f.err("f.impute_loess_logit: span. and config$impute_span unset", 
      config=config)
  }
  if(is.null(n_pts)) n_pts <- config$impute_n_pts
  if(is.null(n_pts)) {
    f.err("f.impute_loess_logit: n_pts and config$impute_n_pts unset", 
      config=config)
  }

  m <- apply(state$expression, 1, f_mid, na.rm=T)
  m[is.na(m)] <- 0
  f <- function(v) sum(is.na(v) | v %in% 0)
  n0 <- apply(state$expression, 1, f)                        ## n.missing
  n1 <- ncol(state$expression) - n0                          ## n.found
  p <- (n0 + off) / (ncol(state$expression) + off)           ## p.missing
  dat <- data.frame(n0=n0, n1=n1, p=p, m=m)
  
  fit <- stats::loess(log(p/(1-p)) ~ m, data=dat, span=span., 
    degree=degree, family=fam)
    
  m_new <- seq(from=max(c(state$expression), na.rm=T) / n_pts, 
    to=max(c(state$expression), na.rm=T), length.out=n_pts)
      
  p_hat <- stats::predict(fit, newdata=data.frame(m=m_new))  ## on logit scale
  p_hat[is.na(p_hat)] <- min(p_hat[p_hat > 0], na.rm=T)      ## is.na -> low p
  p_hat = exp(p_hat) / (1 + exp(p_hat))                      ## inverse logit
  i_na <- is.na(c(state$expression))                         ## values to be imputed
  state$expression[i_na] <- sample(m, sum(i_na), replace=T, prob=p)

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

#' Impute missing values using \code{randomForest}.
#' @description
#'   Impute missing values using \code{randomForest}.
#' @details Imputes missing values using a \code{randomForest} model trained 
#'   using observations in which the feature was expressed. Iterates through 
#'   features beginning with those having fewest missing values. If 
#'   \code{aug_steps > 0}, augments observations with affine transormed 
#'   versions. This is meant to enable extrapolation outside of observed 
#'   intensity range. If \code{aug_steps > 0}, assumes expression data have been 
#'   previously \code{log(x+1)} transformed. If you want \code{0} to be 
#'   considered missing, and have \code{0} in the data, do something like 
#'   \code{exprs[exprs \%in\% 0] <- NA} prior to imputing. Augmentation uses 
#'   \code{f.aug_mult()}.
#' @param state A list with elements like that returned by \code{f.read_data()}:
#'   \tabular{ll}{
#'     \code{expression} \cr \tab Numeric matrix with non-negative expression values. \cr
#'     \code{features}   \cr \tab A data.frame with feature meta-data for rows of expression. \cr
#'     \code{samples}    \cr \tab A data.frame with observation meta-data for columns of expression. \cr
#'   } 
#' It is assumed that state$expression has been previously \code{log2(x+1)} transformed.
#' @param config List with configuration values. Uses the following keys:
#'   \tabular{ll}{
#'     \code{log_file}   \cr \tab Path to log file (character); \code{log_file=""} outputs to console.
#'   }
#' @param f_imp Function to use for initial rough imputation.
#' @param ntree Numeric (greater than 0) number of trees in random forest.
#'   \code{randomForest} \code{ntree} parameter.
#' @param mtry Numeric (greater than 0) number of features to sample at each 
#'   node in each tree. \code{randomForest} \code{mtry} parameter.
#' @param aug_mult Numeric affine shift for training data augmentation. 
#' @param aug_add Numeric affine shift for training data augmentation.
#' @param aug_steps Numeric (non-negative) number of augmentation steps. Set to
#'   \code{0} to skip augmentation.
#' @param unlog2 Logical if \code{TRUE}, treats \code{state$expression} as if 
#'   \code{log2} transformed (the normal case). Set to \code{FALSE} for raw 
#'   data.
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
#' exprs <- h0testr::f.sim1(n_obs=20, n_feats=30)$mat
#' exprs <- log2(exprs + 1)
#' feats <- data.frame(feature_id=rownames(exprs))
#' samps <- data.frame(observation_id=colnames(exprs))
#' state <- list(expression=exprs, features=feats, samples=samps)
#' config <- list(log_file="", impute_quantile=0)
#' out <- h0testr::f.impute_rf(state, config, verbose=FALSE)
#' summary(c(state$expression))   ## note number of NAs
#' summary(c(out$state$expression))  ## note number of NAs
#' print(out$log)

f.impute_rf <- function(state, config, f_imp=f.impute_unif_sample_lod, ntree=100, 
    mtry=NULL, aug_mult=0.33, aug_add=0, aug_steps=3, unlog2=T, verbose=T) {
  
  if(!is.matrix(state$expression)) {
    f.err("f.impute_rf: !is.matrix(exprs)", config=config)
  }
  
  n_miss <- apply(state$expression, 1, function(v) sum(is.na(v) | v %in% 0))
  o <- order(n_miss, decreasing=F)
  tbl <- NULL
  
  for(idx_feat in o) {
    
    if(n_miss[idx_feat] < 1) next
    f.msg("processing", rownames(state$expression)[idx_feat], config=config)
    
    x_train <- f_imp(state, config)$expression
    if(any(is.na(c(x_train)))) {
      f.err("f.impute_rf: f_imp returned NAs", config=config)
    }
    if(unlog2) x_train <- 2^x_train - 1
    x_train <- f.augment_affine(x_train, mult=aug_mult, add=aug_add, steps=aug_steps)
    if(unlog2) x_train <- log2(x_train + 1)
    
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
      f.msg("f.impute_rf: y_miss < 0; deferred to f_imp; y_miss:", y_miss, config=config)
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
    f.msg("f.impute_rf: fall-back imputation for", 
      sum(is.na(state$expression)), "features", config=config)
    state <- f_imp(state, config)
  }
  
  return(list(state=state, log=tbl))
}

#' Impute missing values using \code{glmnet}.
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
#' @param state A list with elements like that returned by \code{f.read_data()}:
#'   \tabular{ll}{
#'     \code{expression} \cr \tab Numeric matrix with non-negative expression values. \cr
#'     \code{features}   \cr \tab A data.frame with feature meta-data for rows of expression. \cr
#'     \code{samples}    \cr \tab A data.frame with observation meta-data for columns of expression. \cr
#'   } 
#' It is assumed that state$expression has been previously \code{log2(x+1)} transformed.
#' @param config List with configuration values. Uses the following keys:
#'   \tabular{ll}{
#'     \code{log_file}            \cr \tab Path to log file (character); \code{log_file=""} outputs to console.
#'   }
#' @param f_imp Function to use for initial rough imputation.
#' @param alpha Numeric (between 0 and 1) number of trees in random forest.
#' @param nfolds Numeric (greater than or equal to 2) number of folds for 
#'   cross-validation tuning of \code{lambda} value for \code{glmnet}.
#' @param measure Character in set \code{c("mae", "mse")}. Loss function used 
#'   during cross-validation tuning of \code{lambda} value for \code{glmnet}.
#' @param aug_mult Numeric affine shift for training data augmentation. 
#' @param aug_add Numeric affine shift for training data augmentation.
#' @param aug_steps Numeric (non-negative) number of augmentation steps. Set to
#'   \code{0} to skip augmentation.
#' @param unlog2 Logical if \code{TRUE}, treats \code{state$expression} as if 
#'   \code{log2} transformed (the usual case). Set to \code{FALSE} for
#'   raw data.
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
#' exprs <- h0testr::f.sim1(n_obs=20, n_feats=30)$mat
#' exprs <- log2(exprs + 1)
#' feats <- data.frame(feature_id=rownames(exprs))
#' samps <- data.frame(observation_id=colnames(exprs))
#' state <- list(expression=exprs, features=feats, samples=samps)
#' config <- list(log_file="", impute_quantile=0)
#' out <- h0testr::f.impute_glmnet(state, config, verbose=FALSE)
#' summary(c(state$expression))   ## note number of NAs
#' summary(c(out$state$expression))  ## note number of NAs
#' print(out$log)

f.impute_glmnet <- function(state, config, f_imp=f.impute_unif_sample_lod, 
    nfolds=5, alpha=1, measure="mae", aug_mult=0.33, aug_add=0, 
    aug_steps=3, unlog2=T, verbose=T) {
  
  if(!is.matrix(state$expression)) {
    f.err("f.impute_glmnet: !is.matrix(state$expression)", config=config)
  }
  
  n_miss <- apply(state$expression, 1, function(v) sum(is.na(v) | v %in% 0))
  o <- order(n_miss, decreasing=F)
  tbl <- NULL
  
  for(idx_feat in o) {
    
    if(n_miss[idx_feat] < 1) next
    f.msg("processing", rownames(state$expression)[idx_feat], config=config)
    
    x_train <- f_imp(state, config)$expression
    if(unlog2) x_train <- 2^x_train - 1
    x_train <- f.augment_affine(x_train, mult=aug_mult, add=aug_add, steps=aug_steps)
    if(unlog2) x_train <- log2(x_train + 1)   
    
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
      f.msg("f.impute_glmnet: y_miss < 0; deferred to f_imp; y_miss:", y_miss, config=config)
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
    f.msg("f.impute_glmnet: fall-back imputation for", 
      sum(is.na(state$expression)), "features", config=config)
    state <- f_imp(state, config)
  }
  
  return(list(state=state, log=tbl))
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
#'     \code{exprs[exprs \%in\% 0] <- 0} prior to imputing. 
#' @param state A list with elements like that returned by \code{f.read_data()}:
#'   \tabular{ll}{
#'     \code{expression} \cr \tab Numeric matrix with non-negative expression values. \cr
#'     \code{features}   \cr \tab A data.frame with feature meta-data for rows of expression. \cr
#'     \code{samples}    \cr \tab A data.frame with observation meta-data for columns of expression. \cr
#'   }
#' It is assumed that \code{state$expression} has been previously \code{log2(x+1)} transformed.
#' @param config List with configuration values. Uses the following keys:
#'   \tabular{ll}{
#'     \code{log_file}         \cr \tab Path to log file (character); \code{log_file=""} outputs to console. \cr
#'     \code{impute_method}    \cr \tab In \code{c("unif_global_lod","unif_sample_lod","sample_lod","rnorm_feature","glm_binom","loess_logit","glmnet","rf","none")} \cr
#'     \code{impute_quantile}  \cr \tab Quantile of signal distribution to use as estimated limit of detection (LOD). \cr
#'     \code{impute_scale}     \cr \tab Factor (numeric) rescaling variance of normal distribution from which draws are made for \code{impute_method \%in\% "rnorm"}. \cr
#'     \code{impute_n_pts}     \cr \tab Numeric greater than one. Determines granularity of imputation. Larger values lead to finer grain. \cr
#'     \code{impute_span}      \cr \tab Span (numeric between zero and 1) for \code{loess} fit. \cr
#'   }
#' @return Updated \code{state} list with the following elements:
#'   \tabular{ll}{
#'     \code{expression} \cr \tab Numeric matrix with non-negative expression values. \cr
#'     \code{features}   \cr \tab A data.frame with feature meta-data for rows of expression. \cr
#'     \code{samples}    \cr \tab A data.frame with observation meta-data for columns of expression. \cr
#'   } 
#' @examples
#' set.seed(101)
#' exprs <- h0testr::f.sim1(n_obs=8, n_feats=12)$mat
#' exprs <- log2(exprs + 1)
#' feats <- data.frame(feature_id=rownames(exprs))
#' samps <- data.frame(observation_id=colnames(exprs))
#' state <- list(expression=exprs, features=feats, samples=samps)
#' config <- list(log_file="", impute_method="unif_sample_lod", impute_quantile=0, save_state=FALSE)
#' out <- h0testr::f.impute(state, config)
#' summary(c(state$expression))      ## note number of NAs
#' summary(c(out$state$expression))  ## note number of NAs

f.impute <- function(state, config) {

  if(config$impute_method %in% "unif_global_lod") {
    state <- f.impute_unif_global_lod(state, config)
  } else if(config$impute_method %in% "sample_lod") {
    state <- f.impute_sample_lod(state, config)
  } else if(config$impute_method %in% "unif_sample_lod") {
    state <- f.impute_unif_sample_lod(state, config)
  } else if(config$impute_method %in% "rnorm_feature") {
    state <- f.impute_rnorm_feature(state, config)
  } else if(config$impute_method %in% "glm_binom") {
    state <- f.impute_glm_binom(state, config)
  } else if(config$impute_method %in% "loess_logit") {
    state <- f.impute_loess_logit(state, config)
  } else if(config$impute_method %in% "glmnet") {
    out <- f.impute_glmnet(state, config)
    state <- out$state
  } else if(config$impute_method %in% "rf") {
    out <- f.impute_rf(state, config)
    state <- out$state
  } else if(config$impute_method %in% "none") {
    f.msg("skipping imputation: config$impute_method %in% 'none'", 
      config=config)
  } else {
    f.err("f.impute: unexpected config$impute_method:", 
      config$impute_method, config=config)
  }
  
  f.check_state(state, config)
  f.report_state(state, config)
  i <- config$run_order %in% "impute"
  prfx <- paste0(which(i)[1] + 2, ".imputed")
  f.save_state(state, config, prefix=prfx)
  
  return(list(state=state, config=config))
}

