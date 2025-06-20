## helper for f.sim1 and f.sim2:

f.sim0 <- function(n_obs, feat_means, feat_sds, mnar_c0, mnar_c1, mnar_off=0.0001) {

  mat <- NULL

  for(i_obs in 1:n_obs) {

    v <- feat_means + stats::rnorm(length(feat_means), mean=0, sd=feat_sds)
    v[v < 0] <- 0

    resp <- mnar_c0 + mnar_c1 * log(v + mnar_off)   ## logit(p_mnar) ~ c0 + c1 * log(intensity)
    p_mnar = exp(resp) / (1 + exp(resp))            ## inverse logit

    i_mnar <- as.logical(stats::rbinom(length(v), 1, p_mnar))
    v[i_mnar] <- NA

    mat <- cbind(mat, v)
  }

  return(mat)
}

#' Simulate a one-condition dataset
#' @description 
#'   Simulates non-negative data matrix with \code{n_feats} rows and 
#    \code{n_obs} columns. 
#' @details 
#'   Feature means: \code{log(feature_mean) ~ rnorm(mean=log_m_mean, sd=log_m_sd)}.
#'   Dispersion of feature CVs: \code{log(feature_cv) ~ rnorm(mean=log_cv_mean, sd=log_cv_sd)}.
#'   MNAR: \code{logit(p(mnar|log(m))) ~ mnar_c0 + mnar_c1 * log(m + mnar_off)}.
#'   For no MNAR, set \code{mnar_c0=-Inf, mnar_c1=0}.
#'   For no MCAR, set \code{mcar_p=0}.
#' @param n_obs Number of observations to simulate. Scalar numeric, with \code{n_obs >= 1}.
#' @param n_feats Number of features to simulate. Scalar numeric, with \code{n_feats >= 1}.
#' @param log_m_mean Mean of feature means. Scalar numeric, with \code{log_m_mean > 0}.
#' @param log_m_sd Standard deviation of feature means around \code{log_m_mean}. 
#'   Scalar numeric, with \code{log_m_sd > 0}.
#' @param log_cv_mean Mean coefficient of variation of features around respective 
#'   feature means. Scalar numeric, with \code{log_cv_mean > 0}.
#' @param log_cv_sd Standard deviation of coefficient of variation of features 
#'   around respective feature means. Scalar numeric, with \code{log_cv_sd > 0}.
#' @param mnar_c0 Intercept of logistic fit of \code{p(mnar|log(intensity))}. Scalar numeric.
#' @param mnar_c1 Slope of logistic fit of \code{p(mnar|log(intensity))}. Scalar numeric.
#' @param mnar_off Offset for taking logit of \code{p(mnar|log(intensity))}. 
#'   Scalar numeric, with \code{0 < mnar_off < 1}.
#' @param mcar_p Probability of missing completely at random. 
#'   Numeric, with \code{0 <= mcar_p <= 1};
#' @return List with elements:
#'   \tabular{ll}{
#'     \code{mat}       \cr \tab Non-negative numeric matrix with \code{n_feats} rows and \code{n_obs} columns. \cr
#'     \code{feat_mean} \cr \tab Numeric vector of non-negative parameter means for each feature in \code{mat}. \cr
#'     \code{feat_cv}   \cr \tab Numeric vector of non-negative parameter CVs for each feature in \code{mat}. \cr
#'   }
#' @examples
#' ## default missing value settings:
#' rslt <- h0testr::f.sim1(n_obs=6, n_feats=8)
#' print(rslt$mat)
#' print(rslt$feat_mean)
#' print(rslt$feat_cv)
#'
#' ## no missing values:
#' rslt <- h0testr::f.sim1(n_obs=6, n_feats=8, mnar_c0=-Inf, mnar_c1=0, mcar_p=0)
#' print(rslt$mat)
#' print(rslt$feat_mean)
#' print(rslt$feat_cv)

f.sim1 <- function(n_obs, n_feats, log_m_mean=11, log_m_sd=2.7,
    log_cv_mean=-0.75, log_cv_sd=0.5, mnar_c0=4.65, mnar_c1=-0.5,
    mnar_off=0.0001, mcar_p=0.002) {

  ## feature mean, cv, and sd:
  m <- exp(stats::rnorm(n=n_feats, mean=log_m_mean, sd=log_m_sd))
  cv <- exp(stats::rnorm(n=n_feats, mean=log_cv_mean, sd=log_cv_sd))
  s <- m * cv
  
  mat <- f.sim0(n_obs=n_obs, feat_means=m, feat_sds=s, 
    mnar_c0=mnar_c0, mnar_c1=mnar_c1, mnar_off=mnar_off)

  rownames(mat) <- names(m) <- names(cv) <- paste0("feat_", 1:n_feats)
  colnames(mat) <- paste0("obs_", 1:n_obs)

  i_mcar <- as.logical(stats::rbinom(length(mat), 1, mcar_p))
  i_mcar <- matrix(i_mcar, nrow=nrow(mat), ncol=ncol(mat))
  mat[i_mcar] <- NA

  return(list(mat=mat, feat_mean=m, feat_cv=cv))
}

#' Simulate a two-condition dataset
#' @description 
#'   Simulates non-negative data matrix with \code{n_feats} rows and 
#    \code{n_obs1 + n_obs2} columns. 
#' @details 
#'   For \code{n_sig} features out of \code{n_feats} will have a 50:50 chance 
#'     of an increase or a decrease. 
#'   For an increase, \code{mean2 <- mean1 * 2^(fold_change)}.
#'   For a decrease, \code{mean2 <- mean1 * 2^(-fold_change)}.
#'   Feature means: 
#'     \code{log(feature_mean) ~ rnorm(mean=log_m_mean, sd=log_m_sd)}.
#'   Dispersion of feature CVs: 
#'     \code{log(feature_cv) ~ rnorm(mean=log_cv_mean, sd=log_cv_sd)}.
#'   MNAR: \code{logit(p(mnar|log(m))) ~ mnar_c0 + mnar_c1 * log(m + mnar_off)}.
#'   For no MNAR, set \code{mnar_c0=-Inf, mnar_c1=0}.
#'   For no MCAR, set \code{mcar_p=0}.
#' @param n_obs1 Number of observations to simulate for condition1. 
#'   Scalar numeric, with \code{n_obs1 >= 1}.
#' @param n_obs2 Number of observations to simulate for condition2. 
#'   Scalar numeric, with \code{n_obs2 >= 1}.
#' @param n_feats Number of features to simulate. Scalar numeric, 
#'   with \code{n_feats >= 1}.
#' @param n_sig Number of significant (affected) features. Scalar numeric, 
#'   with \code{n_sig >= 0}.
#' @param fold_change Fold change for \code{n_sig} features. Scalar numeric.
#' @param log_m_mean Mean of feature means. Scalar numeric, 
#'   with \code{log_m_mean > 0}.
#' @param log_m_sd Standard deviation of feature means around 
#'   \code{log_m_mean}. Scalar numeric, with \code{log_m_sd > 0}.
#' @param log_cv_mean Mean coefficient of variation of features around respective 
#'   feature means. Scalar numeric, with \code{log_cv_mean > 0}.
#' @param log_cv_sd Standard deviation of coefficient of variation of features 
#'   around respective feature means. Scalar numeric, with \code{log_cv_sd > 0}.
#' @param mnar_c0 Intercept of logistic fit of \code{p(mnar|log(intensity))}. 
#'   Scalar numeric.
#' @param mnar_c1 Slope of logistic fit of \code{p(mnar|log(intensity))}. 
#'   Scalar numeric.
#' @param mnar_off Offset for taking logit of \code{p(mnar|log(intensity))}. 
#'   Scalar numeric, with \code{0 < mnar_off < 1}.
#' @param mcar_p Probability of missing completely at random. 
#'   Numeric, with \code{0 <= mcar_p <= 1};
#' @return List with elements:
#'   \tabular{ll}{
#'     \code{mat}       \cr \tab Non-negative numeric matrix with \code{n_feats} rows and \code{n_obs} columns. \cr
#'     \code{feat_mean} \cr \tab Numeric vector of non-negative parameter means for each feature in \code{mat}. \cr
#'     \code{feat_cv}   \cr \tab Numeric vector of non-negative parameter CVs for each feature in \code{mat}. \cr
#'   }
#' @examples
#' ## default missing value settings:
#' rslt <- h0testr::f.sim2(n_obs1=3, n_obs2=3, n_feats=8, n_sig=2)
#' print(rslt)
#'
#' ## no missing values:
#' rslt <- h0testr::f.sim2(n_obs1=3, n_obs2=3, n_feats=8, n_sig=2, mnar_c0=-Inf, mnar_c1=0, mcar_p=0)
#' print(rslt)

f.sim2 <- function(n_obs1, n_obs2, n_feats, n_sig, fold_change=0.5, 
    log_m_mean=11, log_m_sd=2.7, log_cv_mean=-0.75, log_cv_sd=0.5, 
    mnar_c0=4.65, mnar_c1=-0.5, mnar_off=0.0001, mcar_p=0.002) {

  ## feature mean, cv, and sds:
  cv <- exp(stats::rnorm(n=n_feats, mean=log_cv_mean, sd=log_cv_sd))
  m1 <- exp(stats::rnorm(n=n_feats, mean=log_m_mean, sd=log_m_sd))
  s1 <- m1 * cv

  ## group "a":
  mat1 <- f.sim0(n_obs=n_obs1, feat_means=m1, feat_sds=s1, 
    mnar_c0=mnar_c0, mnar_c1=mnar_c1, mnar_off=mnar_off)
  colnames(mat1) <- paste0("obs_a_", 1:n_obs1)

  ## effects:
  m2 <- m1
  idx_sig <- sample(1:length(m1), n_sig)
  i_neg <- as.logical(stats::rbinom(n_sig, 1, prob=0.5))
  m2[idx_sig][i_neg]  <- m2[idx_sig][i_neg] * 2^(-fold_change)
  m2[idx_sig][!i_neg] <- m2[idx_sig][!i_neg] * 2^(fold_change)
  s2 <- m2 * cv

  ## group "b":
  mat2 <- f.sim0(n_obs=n_obs2, feat_means=m2, feat_sds=s2, 
    mnar_c0=mnar_c0, mnar_c1=mnar_c1, mnar_off=mnar_off)
  colnames(mat2) <- paste0("obs_b_", 1:n_obs2)

  mat <- cbind(mat1, mat2)
  rownames(mat) <- paste0("feat_", 1:n_feats)

  i_mcar <- as.logical(stats::rbinom(length(mat), 1, mcar_p))
  i_mcar <- matrix(i_mcar, nrow=nrow(mat), ncol=ncol(mat))
  mat[i_mcar] <- NA

  return(mat=mat)
}
