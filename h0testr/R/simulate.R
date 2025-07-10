## helper for f.sim1 and f.sim2; ensures all values strictly positive (> 0):

f.sim0 <- function(n_obs, feat_means, feat_sds, mnar_c0, mnar_c1, mnar_off=0.0001) {

  mat <- NULL

  for(i_obs in 1:n_obs) {

    v <- stats::rnorm(length(feat_means), mean=feat_means, sd=feat_sds)
    i <- v <= 0
    while(any(i)) {
      v[i] <- stats::rnorm(sum(i), mean=feat_means[i], sd=feat_sds[i])
      i <- v <= 0
    }

    mat <- cbind(mat, v)
  }

  return(mat)
}

## helper for f.sim1 and f.sim2; estimate p(missing|log(intensity)) using 
##   logit model. Drops cells in mat randomly based on p(missing|log(intensity)):

f.mnar <- function(mat, mnar_c0, mnar_c1, mnar_off=0.0001) {

  f <- function(v) {

    resp <- mnar_c0 + mnar_c1 * log(v + mnar_off)   ## logit(p_mnar) ~ c0 + c1 * log(intensity)
    p_mnar = exp(resp) / (1 + exp(resp))            ## inverse logit

    i_mnar <- as.logical(stats::rbinom(length(v), 1, p_mnar))
    v[i_mnar] <- NA
    
    return(v)
  }
  
  return(apply(mat, 2, f))
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
#' @param n_obs Scalar number of observations to simulate, where \code{n_obs >= 1}.
#' @param n_feats Scalar number of features to simulate, where \code{n_feats >= 1}.
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
  
  ## mnar:
  mat <- f.mnar(mat, mnar_c0=mnar_c0, mnar_c1=mnar_c1)
  
  ## mcar:
  i_mcar <- as.logical(stats::rbinom(length(mat), 1, mcar_p))
  i_mcar <- matrix(i_mcar, nrow=nrow(mat), ncol=ncol(mat))
  mat[i_mcar] <- NA

  return(list(mat=round(mat), feat_mean=m, feat_cv=cv))
}

#' Simulate a two-condition dataset
#' @description 
#'   Simulates non-negative data matrix with \code{n_genes * peps_per_gene} 
#'   rows and \code{(n_samps1 + n_samps2) * reps_per_sample} columns. 
#' @details 
#'   For \code{n_genes_signif} features out of \code{n_genes} will have a 50:50 chance 
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
#' @param n_samps1 Scalar number of samples to simulate for condition1, with \code{n_samps1 >= 1}. 
#' @param n_samps2 Scalar number of samples to simulate for condition2, with \code{n_samps2 >= 1}. 
#' @param n_genes Scalar number of gene groups to simulate, with \code{n_genes >= 2}. 
#' @param n_genes_signif Number of significant (affected) genes. Scalar numeric, 
#'   with \code{0 <= n_genes_signif <= n_genes}.
#' @param fold_change Fold change for \code{n_sig} features. Scalar numeric.
#' @param peps_per_gene Number of peptides to simulate per gene. 
#' @param reps_per_sample Number of technical replicates to simulate per sample.
#' @param cv_reps Coeffient of variation between technical replicates.
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
#' rslt <- h0testr::f.sim2(n_samps1=3, n_samps2=3, n_genes=8, n_genes_signif=2)
#' print(rslt)
#' 
#' ## no missing values:
#' rslt <- h0testr::f.sim2(n_samps1=3, n_samps2=3, n_genes=8, n_genes_signif=2, mnar_c0=-Inf, mnar_c1=0, mcar_p=0)
#' print(rslt)
#' 
#' ## peptide-level, with technical replicate observations:
#' rslt <- h0testr::f.sim2(n_samps1=2, n_samps2=2, n_genes=8, reps_per_sample=2, peps_per_gene=3)
#' print(rslt)
#'
#' ## minimal \code{2x2} matrix with no effects:
#' rslt <- h0testr::f.sim2(n_samps1=1, n_samps2=1, n_genes=2)
#' print(rslt)

f.sim2 <- function(n_samps1, n_samps2, n_genes, n_genes_signif=0, 
    fold_change=0.5, peps_per_gene=1, reps_per_sample=1, cv_reps=0.1, 
    log_m_mean=11, log_m_sd=2.7, log_cv_mean=-0.75, log_cv_sd=0.5, 
    mnar_c0=4.65, mnar_c1=-0.5, mnar_off=0.0001, mcar_p=0.002) {
    
  if(n_samps1 < 1 || n_samps2 < 1) stop("n_samps1 < 1 || n_samps2 < 1")
  if(n_genes < 2 || n_genes_signif < 0) stop("n_genes < 2 || n_genes_signif < 0")
  
  ## feature mean, cv, and sds:
  cv <- exp(stats::rnorm(n=n_genes * peps_per_gene, mean=log_cv_mean, sd=log_cv_sd))
  m1 <- exp(stats::rnorm(n=n_genes * peps_per_gene, mean=log_m_mean, sd=log_m_sd))
  s1 <- m1 * cv
  
  ## genes, significance, peptides, and labels:
  gene <- paste0("gene", 1:n_genes)
  sig_genes <- sample(gene, n_genes_signif)
  if(peps_per_gene >= 2) {
    pep_labels <- do.call(paste0, expand.grid(paste0("_pep", 1:peps_per_gene), gene)[, c(2, 1)])
    gene_labels <- sub("\\_.*", "", pep_labels)
  } else {
    gene_labels <- pep_labels <- gene
  }
  
  ## group "grp1":
  mat1 <- f.sim0(n_obs=n_samps1, feat_means=m1, feat_sds=s1)
  colnames(mat1) <- paste0("grp1_samp", 1:n_samps1)
  
  ## effects; equal chance of up or down:
  m2 <- m1
  for(gene in sig_genes) {
    i_gene <- gene_labels %in% gene
    if(as.logical(stats::rbinom(1, 1, prob=0.5))) {
      m2[i_gene] <- m2[i_gene] * 2^(fold_change)
    } else {
      m2[i_gene] <- m2[i_gene] * 2^(-fold_change)
    }
  }
  s2 <- m2 * cv
  
  ## group "grp2":
  mat2 <- f.sim0(n_obs=n_samps2, feat_means=m2, feat_sds=s2)
  colnames(mat2) <- paste0("grp2_samp", 1:n_samps2)
  
  mat <- cbind(mat1, mat2)
  rownames(mat) <- pep_labels
  
  ## technical replication:
  if(reps_per_sample >= 2) {
    
    f <- function(v) {
      f0 <- function(val) {
        if(is.na(val)) {
          return(rep(NA, reps_per_sample))
        } else {
          return(stats::rnorm(reps_per_sample, mean=val, sd=cv_reps * val))
        }
      }
      return(list(t(sapply(v, f0))))
    }
    tmp_list <- apply(mat, 2, f)
    
    for(nom in names(tmp_list)) {
      tmp_list[[nom]] <- tmp_list[[nom]][[1]]
      colnames(tmp_list[[nom]]) <- paste0(nom, "_rep", 1:reps_per_sample)
    }
    mat <- do.call(cbind, tmp_list)
  }
  
  ## mnar:
  mat <- f.mnar(mat, mnar_c0=mnar_c0, mnar_c1=mnar_c1)
  
  ## mcar:
  i_mcar <- as.logical(stats::rbinom(length(mat), 1, mcar_p))
  i_mcar <- matrix(i_mcar, nrow=nrow(mat), ncol=ncol(mat))
  mat[i_mcar] <- NA
  
  return(list(mat=round(mat), feat_mean1=m1, feat_mean2=m2, feat_cv=cv))
}
