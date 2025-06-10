#' Impute missing values from between 0 and global LOD.
#' @description
#' `f.impute_unif_global_lod` imputes missing values by randomly drawing 
#'   from uniform distribution below an estimated global LOD. 
#' @details Imputed values are random draws from uniform distribution over the 
#'   interval [0, LOD], where LOD is an estimate global limit of 
#'   detection. Only NA values are considered as missing, so if you want `0` 
#'   to be considered missing, and have `0` in the data, do something like 
#'   `state$expression[state$expression %in% 0] <- NA` prior to imputing.
#' @param state A list with elements like that returned by `f.read_data()`:
#' \itemize{
#'   \item{expression}{Numeric matrix with feature rows and observation columns.}
#'   \item{features}{A data.frame with feature rows corresponding to rows of expression.}
#'   \item{samples}{A data.frame with observation rows corresponding to columns of expression.}
#' }
#' @param config List with configuration settings.
#' @param impute_quantile Numeric between 0 and 1 specifying 
#'   minimum non-zero/NA expression value for each feature to use as global
#'   LOD (limit of detection). If NULL, config$impute_quantile used.
#' @return An updated `state` list with the following elements:
#' \itemize{
#'   \item{expression}{Numeric matrix with feature rows and observation columns.}
#'   \item{features}{A data.frame with feature rows corresponding to rows of expression.}
#'   \item{samples}{A data.frame with observation rows corresponding to columns of expression.}
#' }
#' Returned `state$xpression` matrix contains strictly positive numeric values.
#' @examples
#' exprs[exprs %in% 0] <- NA
#' state <- list(expression=exprs)
#' config <- list(impute_quantile=0.05, log_file="")
#' state2 <- f.impute_unif_global_lod(state, config)
#' exprs2 <- state2$expression
#'
#' config <- list(log_file="")
#' state2 <- f.impute_unif_global_lod(state, config, impute_quantile=0.05)
#' exprs2 <- state2$expression

f.impute_unif_global_lod <- function(state, config, impute_quantile=NULL) {

  if(!is.matrix(state$expression)) {
     f.err("f.impute_unif_global_lod: !is.matrix(state$expression)", config=config)
  }
  if(is.null(impute_quantile)) impute_quantile <- config$impute_quantile
  
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
#' `f.impute_unif_sample_lod` imputes missing values by randomly drawing 
#'   from uniform distribution below an estimated observation-specific LOD. 
#' @details Imputed values are random draws from uniform distribution over the 
#'   interval [0, LOD], where LOD is an estimate observation limit of 
#'   detection. Only NA values are considered as missing, so if you want `0` 
#'   to be considered missing, and have `0` in the data, do something like 
#'   `exprs[exprs %in% 0] <- NA` prior to imputing.
#' @param state A list with elements like that returned by `f.read_data()`:
#' \itemize{
#'   \item{expression}{Numeric matrix with feature rows and observation columns.}
#'   \item{features}{A data.frame with feature rows corresponding to rows of expression.}
#'   \item{samples}{A data.frame with observation rows corresponding to columns of expression.}
#' }
#' @param config List with configuration settings.
#' @param impute_quantile Numeric between 0 and 1 specifying 
#'   minimum non-zero/NA expression value in each observation to use as the
#'   observation-specific LOD (limit of detection).
#' @return An updated `state` list with the following elements:
#' \itemize{
#'   \item{expression}{Numeric matrix with feature rows and observation columns.}
#'   \item{features}{A data.frame with feature rows corresponding to rows of expression.}
#'   \item{samples}{A data.frame with observation rows corresponding to columns of expression.}
#' }
#' @examples
#' exprs[exprs %in% 0] <- NA
#' state <- list(expression=exprs)
#' config <- list(log_file="")
#' state2 <- f.impute_unif_sample_lod(state, config, impute_quantile=0.01)
#' exprs2 <- state2$expression
#'
#' config <- list(impute_quantile=0.01, log_file="")
#' state2 <- f.impute_unif_sample_lod(state, config)
#' exprs2 <- state2$expression

f.impute_unif_sample_lod <- function(state, config, impute_quantile=0) {

  if(!is.matrix(state$expression)) {
    f.err("f.impute_unif_sample_lod: !is.matrix(state$expression)", config=config)
  }
  if(is.null(impute_quantile)) impute_quantile <- config$impute_quantile
  
  f <- function(v) {
    i <- is.na(v)
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
#' `f.impute_sample_lod` imputes missing values as the minimum observed
#'    value in the corresponding observation.
#' @details Missing values are set to an estimate observation limit of 
#'   detection (LOD). Only NA values are considered as missing, so if you 
#'   want `0` to be considered missing, and have `0` in the data, do something 
#'   like `exprs[exprs %in% 0] <- NA` prior to imputing.
#' @param state A list with elements like that returned by `f.read_data()`:
#' \itemize{
#'   \item{expression}{Numeric matrix with feature rows and observation columns.}
#'   \item{features}{A data.frame with feature rows corresponding to rows of expression.}
#'   \item{samples}{A data.frame with observation rows corresponding to columns of expression.}
#' }
#' @param config List with configuration settings.
#' @return An updated `state` list with the following elements:
#' \itemize{
#'   \item{expression}{Numeric matrix with feature rows and observation columns.}
#'   \item{features}{A data.frame with feature rows corresponding to rows of expression.}
#'   \item{samples}{A data.frame with observation rows corresponding to columns of expression.}
#' }
#' @examples
#' exprs[exprs %in% 0] <- NA
#' state <- list(expression=exprs)
#' config <- list(log_file="")
#' state2 <- f.impute_sample_lod(state, config)
#' exprs2 <- state2$expression

f.impute_sample_lod <- function(state, config) {

  if(!is.matrix(state$expression)) {
    f.err("f.impute_sample_lod: !is.matrix(state$expression)", config=config)
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
#' `f.impute_rnorm_feature` imputes missing values normally distributed 
#'   around the feature mean.
#' @details Missing values are draws from 
#'   `normal(mean=mean(exprs[feature, ]), sd=scale.*sd(exprs[feature, ]))`. 
#'   Only non-NA values are used in calculation of mean and sd. If you  
#'   want `0` to be considered missing, and have `0` in the data, do something 
#'   like `exprs[exprs %in% 0] <- NA` prior to imputing. All values are
#'   gauranteed non-negative; negative draws are repeated till non-negative
#'   drawn.
#' @param state A list with elements like that returned by `f.read_data()`:
#' \itemize{
#'   \item{expression}{Numeric matrix with feature rows and observation columns.}
#'   \item{features}{A data.frame with feature rows corresponding to rows of expression.}
#'   \item{samples}{A data.frame with observation rows corresponding to columns of expression.}
#' }
#' @param config List with configuration settings.
#' @param scale. Numeric greater than zero, linearly scaling the
#'   dispersion around the feature mean.
#' @return An updated `state` list with the following elements:
#' \itemize{
#'   \item{expression}{Numeric matrix with feature rows and observation columns.}
#'   \item{features}{A data.frame with feature rows corresponding to rows of expression.}
#'   \item{samples}{A data.frame with observation rows corresponding to columns of expression.}
#' }
#' @examples
#' exprs[exprs %in% 0] <- NA
#' state <- list(expression=exprs)
#' config <- list(log_file="")
#' state2 <- f.impute_rnorm_feature(state, config, scale.=1)
#' exprs2 <- state2$expression
#'
#' config <- list(impute_scale=1, log_file="")
#' state2 <- f.impute_rnorm_feature(state, config)
#' exprs2 <- state2$expression

f.impute_rnorm_feature <- function(state, config, scale.=NULL) {

  if(!is.matrix(state$expression)) {
    f.err("f.impute_rnorm_feature: !is.matrix(state$expression)", config=config)
  }
  if(is.null(scale.)) scale. <- config$scale
  
  i <- c(state$expression) < 0
  i[is.na(i)] <- F
  if(any(i)) {
    f.err("f.impute_rnorm_feature: state$expression contains negative values", config=config)
  }
  
  f <- function(v) {
    i <- is.na(v)
    if(any(i)) {
      m <- mean(v, na.rm=T)
      s <- sd(v, na.rm=T) * scale.
      v[i] <- stats::rnorm(sum(i), mean=m, sd=s)
      i <- v < 0
      while(any(i)) {
        v[i] <- stats::rnorm(sum(i), mean=m, sd=s)
        i <- v < 0
      }
    }
    return(v)
  }
  state$expression <- t(apply(state$expression, 1, f))
  
  return(state)
}

#' Impute by drawing from `p(missing|intensity)` from binomial glm.
#' @description
#' `f.impute_glm_binom` imputes missing values by randomly drawing
#'   from an estimated density of `p(missing|intensity)`.
#' @details Impute random draws by drawing from the estimated density of 
#'   `p(missing|intensity)`. Estimates `p(missing|intensity)` using a 
#'   binomial glm with formula `cbind(n.missing, n.found) ~ intensity`,
##   where intensity is assumed to be previously log transformed. If you  
#'   want `0` to be considered missing, and have `0` in the data, do 
#'   something like `exprs[exprs %in% 0] <- NA` prior to imputing. 
#' @param state A list with elements like that returned by `f.read_data()`:
#' \itemize{
#'   \item{expression}{Numeric matrix with feature rows and observation columns.}
#'   \item{features}{A data.frame with feature rows corresponding to rows of expression.}
#'   \item{samples}{A data.frame with observation rows corresponding to columns of expression.}
#' }
#' It is assumed that state$expression has been previously `log2(x+1)` transformed.
#' @param config List with configuration settings.
#' @param gran Numeric greater than zero. Granularity of prediction 
#'   grid. Smaller values lead to less chance of duplicate imputed values.
#'   Larger values require more compute time and memory.
#' @param off Numeric offset for calculating 
#'   `p.missing = (n.missing + off) / (n.total + off)`.
#' @param f_mid Function to use for calculating central tendency of 
#'   feature expression across samples. 
#' @return An updated `state` list with the following elements:
#' \itemize{
#'   \item{expression}{Numeric matrix with feature rows and observation columns.}
#'   \item{features}{A data.frame with feature rows corresponding to rows of expression.}
#'   \item{samples}{A data.frame with observation rows corresponding to columns of expression.}
#' }
#' @examples
#' exprs <- log2(exprs + 1)
#' exprs[exprs %in% 0] <- NA
#' state <- list(expression=exprs)
#' config <- list(log_file="")
#' state2 <- f.impute_glm_binom(state, config, impute_granularity=0.0001)
#' exprs2 <- state2$expression
#'
#' config <- list(gran=0.0001, log_file="")
#' state2 <- f.impute_glm_binom(state, config)
#' exprs2 <- state2$expression

f.impute_glm_binom <- function(state, config, gran=NULL, off=1, f_mid=median) {

  if(!is.matrix(state$expression)) {
    f.err("!is.matrix(state$expression)", config=config)
  }
  if(is.null(gran)) gran <- config$impute_granularity

  m <- apply(state$expression, 1, f_mid, na.rm=T)
  m[is.na(m)] <- 0
  n0 <- apply(state$expression, 1, function(v) sum(is.na(v) | v %in% 0))  ## n.missing
  n1 <- ncol(state$expression) - n0                                     
  p <- (n0 + off) / (ncol(state$expression) + off)
  dat <- data.frame(n0=n0, n1=n1, p=p, m=m)

  fit <- stats::glm(cbind(n0, n1) ~ m, data=dat, family="binomial")
  m <- seq(from=gran, to=max(c(state$expression), na.rm=T), by=gran)
  p <- predict(fit, newdata=data.frame(m=m), type="response")
  p[is.na(p)] <- min(p, na.rm=T)
  i_na <- is.na(c(state$expression))
  state$expression[i_na] <- sample(m, sum(i_na), replace=T, prob=p)

  return(state)
}

#' Impute by drawing from `p(missing|intensity)` estimated with loess.
#' @description
#' `f.impute_loess_logit` imputes missing values by randomly drawing
#'   from an estimated density of `p(missing|intensity)`.
#' @details Impute random draws by drawing from the estimated density of 
#'   `p(missing|intensity)`. Estimates `p(missing|intensity)` using a loess 
#'   fit with form `log(p_missing/p_found) ~ intensity`, where intensity is 
#'   assumed to be previously log transformed. If you want `0` to be 
#'   considered missing, and have `0` in the data, do something like 
#'   `exprs[exprs %in% 0] <- NA` prior to imputing. 
#' @param state A list with elements like that returned by `f.read_data()`:
#' \itemize{
#'   \item{expression}{Numeric matrix with feature rows and observation columns.}
#'   \item{features}{A data.frame with feature rows corresponding to rows of expression.}
#'   \item{samples}{A data.frame with observation rows corresponding to columns of expression.}
#' }
#' It is assumed that state$expression has been previously `log2(x+1)` transformed.
#' @param config List with configuration settings.
#' @param span Numeric greater than 0 and less than 1. 
#'   Span for loess fit.
#' @param gran Numeric greater than zero. Granularity of prediction 
#'   grid. Smaller values lead to less chance of duplicate imputed values.
#'   Larger values require more compute time and memory.
#' @param off Numeric offset for calculating 
#'   `p.missing = (n.missing + off) / (n.total + off)`.
#' @param f_mid Function to use for calculating central tendency of 
#'   feature expression across samples. 
#' @param degree Numeric in `c(1, 2)`. Degree for loess fit.
#' @param fam Character in `c("symmetric", "gaussian")`. Family for loess fit.
#' @return An updated `state` list with the following elements:
#' \itemize{
#'   \item{expression}{Numeric matrix with feature rows and observation columns.}
#'   \item{features}{A data.frame with feature rows corresponding to rows of expression.}
#'   \item{samples}{A data.frame with observation rows corresponding to columns of expression.}
#' }
#' @examples
#' exprs <- log2(exprs + 1)
#' exprs[exprs %in% 0] <- NA
#' state <- list(expression=exprs)
#' config <- list(log_file="")
#' state2 <- f.impute_loess_logit(state, config, span.=1, gran=0.0001)
#' exprs2 <- state2$expression
#'
#' config <- list(span.=1, gran=0.0001, log_file="")
#' state2 <- f.impute_loess_logit(state, config)
#' exprs2 <- state2$expression

f.impute_loess_logit <- function(state, config, span.=NULL, gran=NULL, 
    off=1, f_mid=median, degree=1, fam="symmetric") {

  if(!is.matrix(state$expression)) {
    f.err("!is.matrix(state$expression)", config=config)
  }
  if(is.null(span)) span <- config$impute_span
  if(is.null(gran)) gran <- config$impute_granularity

  m <- apply(state$expression, 1, f_mid, na.rm=T)
  m[is.na(m)] <- 0
  f <- function(v) sum(is.na(v) | v %in% 0)
  n0 <- apply(state$expression, 1, f)                     ## n.missing
  n1 <- ncol(state$expression) - n0                       ## n.found
  p <- (n0 + off) / (ncol(exprs) + off)                   ## p.missing
  dat <- data.frame(n0=n0, n1=n1, p=p, m=m)

  m_new <- seq(from=gran, to=max(c(state$expression), na.rm=T), by=gran)
  
  fit <- stats::loess(log(p/(1-p)) ~ m, data=dat, span=span, 
    degree=degree, family=fam)
    
  p_hat <- predict(fit, newdata=data.frame(m=m_new))      ## on logit scale
  p_hat[is.na(p_hat)] <- min(p_hat[p_hat > 0], na.rm=T)   ## is.na -> low p
  p_hat = exp(p_hat) / (1 + exp(p_hat))                   ## inverse logit
  i_na <- is.na(c(state$expression))                      ## values to be imputed
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

#' Impute missing values using `randomForest`.
#' @description
#' `f.impute_rf` imputes missing values using `randomForest`.
#' @details Imputes missing values using a `randomForest` model trained using 
#'   observations in which the feature was expressed. Iterates through 
#'   features beginning with those having fewest missing values. If 
#'   `aug_steps > 0`, augments observations with affine transormed versions.
#'   This is meant to enable extrapolation outside of observed intensity 
#'   range. If `aug_steps > 0`, assumes expression data have been previously
#'   `log(x+1)` transformed. If you want `0` to be considered missing, and 
#'   have `0` in the data, do something like `exprs[exprs %in% 0] <- NA` 
#'   prior to imputing. Augmentation uses `f.aug_mult()`.
#' @param state A list with elements like that returned by `f.read_data()`:
#' \itemize{
#'   \item{expression}{Numeric matrix with feature rows and observation columns.}
#'   \item{features}{A data.frame with feature rows corresponding to rows of expression.}
#'   \item{samples}{A data.frame with observation rows corresponding to columns of expression.}
#' }
#' It is assumed that state$expression has been previously `log2(x+1)` transformed.
#' @param config List with configuration settings.
#' @param f_imp Function to use for initial rough imputation.
#' @param ntree Numeric (greater than 0) number of trees in random forest.
#' @param mtry Numeric (greater than 0) number of features to sample at each 
#'   node in each tree. Random forest mtry.
#' @param aug_mult Numeric affine shift for training data augmentation. 
#' @param aug_add Numeric affine shift for training data augmentation.
#' @param aug_steps Numeric (non-negative) number of augmentation steps. Set to
#'   `0` to skip augmentation.
#' @param verbose Logical if TRUE, emits progress messages.
#' @return A list with the following elements:
#' \itemize{
#'   \item{exprs} \itemize{
#'     \item{expression}{Numeric matrix with feature rows and observation columns.}
#'     \item{features}{A data.frame with feature rows corresponding to rows of returned expression.}
#'     \item{samples}{A data.frame with observation rows corresponding to columns of returned expression.}
#'   }
#'   \item{log}{A data.frame logging statistics for each fit.}
#' }
#' @examples
#' exprs <- log2(exprs + 1)
#' exprs[exprs %in% 0] <- NA
#' state <- list(expression=exprs)
#' config <- list(log_file="")
#' state2 <- f.impute_rf(state, config)
#' exprs2 <- state2$expression
#'
#' config <- list(ntree=250, aug_steps=2, log_file="")
#' state2 <- f.impute_rf(state, config)
#' exprs2 <- state2$expression

f.impute_rf <- function(state, config, f_imp=f.impute_loess_logit, ntree=100, 
    mtry=NULL, aug_mult=0.33, aug_add=0, aug_steps=3, verbose=T) {

  if(!is.matrix(state$expression)) f.err("!is.matrix(exprs)", config=config)

  n_miss <- apply(state$expression, 1, function(v) sum(is.na(v) | v %in% 0))
  o <- order(n_miss, decreasing=F)
  tbl <- NULL
  
  for(idx_feat in o) {
  
    if(n_miss[idx_feat] < 1) next
    f.msg("processing", rownames(state$expression)[idx_feat], config=config)

    x_train <- f_imp(state$expression)
    x_train <- 2^x_train - 1
    x_train <- f.augment_affine(x_train, mult=aug_mult, add=aug_add, steps=aug_steps)
    x_train <- log2(x_train + 1)

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
    y_miss <- predict(fit, newdata=t(x_miss), type="response")
    if(any(y_miss < 0)) f.err("ERROR: y_miss < 0; y_miss:", y_miss, config=config)
    exprs[idx_feat, i_miss] <- y_miss

    tm_stmp <- format(Sys.time(), format='%Y%m%d%H%M%S')
    tbl_i <- data.frame(time=tm_stmp, feat=rownames(exprs)[idx_feat], 
      mtry=mtry0, n_miss=sum(i_miss), mean0=mean(y[!i_miss]), 
      mean1=mean(y_miss), pvar0=fit$rsq[1], pvar1=fit$rsq[ntree], 
      mse0=fit$mse[1], mse1=fit$mse[ntree]
    )
    
    if(verbose) f.msg(tbl_i, config=config)
    tbl <- rbind(tbl, tbl_i)
  }

  return(list(state=state, log=tbl))
}

#' Impute missing values using `glmnet`.
#' @description
#' `f.impute_rf` imputes missing values using `glmnet`.
#' @details Imputes missing values using a `glmnet::cv.glmnet` model trained 
#'   using observations in which the feature was expressed. Iterates through 
#'   features beginning with those having fewest missing values. If 
#'   `aug_steps > 0`, augments observations with affine transormed versions.
#'   This is meant to enable extrapolation outside of observed intensity 
#'   range. If `aug_steps > 0`, assumes expression data have been previously
#'   `log(x+1)` transformed. If you want `0` to be considered missing, and 
#'   have `0` in the data, do something like `exprs[exprs %in% 0] <- NA` 
#'   prior to imputing. Augmentation uses `f.aug_mult()`.
#' @param state A list with elements like that returned by `f.read_data()`:
#' \itemize{
#'   \item{expression}{Numeric matrix with feature rows and observation columns.}
#'   \item{features}{A data.frame with feature rows corresponding to rows of expression.}
#'   \item{samples}{A data.frame with observation rows corresponding to columns of expression.}
#' }
#' It is assumed that state$expression has been previously `log2(x+1)` transformed.
#' @param config List with configuration settings.
#' @param f_imp Function to use for initial rough imputation.
#' @param alpha Numeric (between 0 and 1) number of trees in random forest.
#' @param nfolds Numeric (greater than or equal to 2) number of folds for 
#'   cross-validation tuning of `lambda` value for `glmnet`.
#' @param measure Character in `c("mae", "mse")`. Loss function used during 
#'   cross-validation tuning of `lambda` value for `glmnet`.
#' @param aug_mult Numeric affine shift for training data augmentation. 
#' @param aug_add Numeric affine shift for training data augmentation.
#' @param aug_steps Numeric (non-negative) number of augmentation steps. Set to
#'   `0` to skip augmentation.
#' @param verbose Logical if TRUE, emits progress messages.
#' @return A list with the following elements:
#' \itemize{
#'   \item{exprs} \itemize{
#'     \item{expression}{Numeric matrix with feature rows and observation columns.}
#'     \item{features}{A data.frame with feature rows corresponding to rows of returned expression.}
#'     \item{samples}{A data.frame with observation rows corresponding to columns of returned expression.}
#'   }
#'   \item{log}{A data.frame logging statistics for each fit.}
#' }
#' @examples
#' exprs[exprs %in% 0] <- NA
#' state <- list(expression=exprs)
#' config <- list(log_file="")
#' state2 <- f.impute_glmnet(state, config, impute_quantile=0.01)
#' exprs2 <- state2$expression
#'
#' state2 <- f.impute_glmnet(state, config, nfolds=3, alpha=0.5, aug_steps=2)
#' state2 <- f.impute_glmnet(state, config)
#' exprs2 <- state2$expression

f.impute_glmnet <- function(state, config, f_imp=f.impute_loess_logit, 
    nfolds=5, alpha=1, measure="mae", aug_mult=0.33, aug_add=0, 
    aug_steps=3, verbose=T) {

  if(!is.matrix(state$expression)) {
    f.err("!is.matrix(state$expression)", config=config)
  }

  n_miss <- apply(state$expression, 1, function(v) sum(is.na(v) | v %in% 0))
  o <- order(n_miss, decreasing=F)
  tbl <- NULL
  
  for(idx_feat in o) {
  
    if(n_miss[idx_feat] < 1) next
    f.msg("processing", rownames(state$expression)[idx_feat], config=config)

    x_train <- f_imp(state$expression)
    x_train <- 2^x_train - 1
    x_train <- f.augment_affine(x_train, mult=aug_mult, add=aug_add, steps=aug_steps)
    x_train <- log2(x_train + 1)

    y <- state$expression[idx_feat, , drop=T]
    i_miss <- is.na(y) | y %in% 0
    col_names_miss <- colnames(state$expression)[i_miss]
    i_miss_train <- colnames(x_train) %in% col_names_miss
    x_train_i <- x_train[-idx_feat, !i_miss_train, drop=F]
    y_train_i <- x_train[idx_feat, !i_miss_train, drop=T]

    fit <- glmnet::cv.glmnet(x=t(x_train_i), y=y_train_i, family="gaussian", alpha=alpha, 
      type.measure=measure, nfolds=nfolds, parallel=F)

    x_miss <- x_train[-idx_feat, which(i_miss), drop=F]
    y_miss <- predict(fit, newx=t(x_miss), s=fit$lambda.1se, type="response")
    if(any(y_miss < 0)) f.err("y_miss < 0; y_miss:", y_miss, config=config)
    state$expression[idx_feat, i_miss] <- y_miss

    idx_lambda <- which(fit$lambda %in% fit$lambda.1se)[1] 
    tm_stmp <- format(Sys.time(), format='%Y%m%d%H%M%S')
    tbl_i <- data.frame(time=tm_stmp, feat=rownames(exprs)[idx_feat], alpha=alpha, nfolds=nfolds, 
      n_miss=sum(i_miss), mean0=mean(y[!i_miss]), mean1=mean(y_miss), 
      cvm0=max(fit$cvm), cvm1=fit$cvm[idx_lambda], cvup0=max(fit$cvup), cvup1=fit$cvup[idx_lambda]
    )
    
    if(verbose) f.msg(tbl_i, config=config)
    tbl <- rbind(tbl, tbl_i)
  }

  return(list(state=state, log=tbl))
}

