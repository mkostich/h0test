###############################################################################
## notes:

## contains some useful functions for proteomics
##   meant to be sourced after defining LOG_FILE:
##     LOG_FILE <- "my.log"        ## log to file
##     LOG_FILE <- ""              ## log to console
##     source("/path/to/h0test_frm.R")

## uses limma and edgeR packages

## usual flow:
##   1) prefilter for all NA or all 0 rows (features) or columns (observations)
##   2) normalize, and if not part of normalization, transform (e.g. log2)
##   3) combine technical replicates (e.g. take median)
##   4) filter features and observations based on number of values > 0
##   5) impute missing values
##   6) test hypotheses, etc.

## f.msg(...)
## 
## f.log(...)
## 
## f.err(...)
## 
## f.save_tsv(dat, file_out, row.names=T, col.names=T)
## 
## f.filter_features(exprs, samps, feats, n_samples_min=2)
##   return(list(exprs=exprs, samps=samps, feats=feats))
##
## f.filter_samples(exprs, samps, feats, n_features_min=1000)
##   return(list(exprs=exprs, samps=samps, feats=feats))
## 
## f.samples_per_feature(exprs) return(n)
##
## f.feature_median_expression(exprs) return(m)
##
## f.features_per_sample(exprs) return(n)
##
## f.normalize_edger(exprs, method="TMM") return(exprs)
##   Inter-sample edger (tmm) normalization;
##     normalized cpm values using edger norm factors
##   method in c("TMM", "TMMwsp", "RLE", "upperquartile", "none")
##
## f.normalize_quantile <- function(exprs, norm_quantile=0.5, multiplier=1e6) return(exprs)
##   Inter-sample quantile normalization; same as limma scale normalization 
##     when norm_quantile == 0.5 (median); multiplier like with cpm:
## 
## Counts per million (if multiplier is 1e6):
##   f.normalize_cpm(exprs, multiplier=1e6) return(exprs)
## 
## f.normalize_vsn(exprs) return(exprs)
##   Variance stabilizing normalization: results on log2 scale:
## 
## f.normalize_loess <- function(exprs, span=0.7, method="affy") return(exprs)
##   Cyclic loess normalization; 
##     method in c("fast", "affy", "pairs"); 
##       slow; fast is linear, others are quadratic:
##
## f.impute_unif_global_lod(exprs, impute_quantile=0.01) return(exprs)
##   Impute missing values assuming left-censored only.
##     Returns expression matrix same as exprs but w/ missing values imputed.
##       ONLY NA VALUES CONSIDERED MISSING, so may need to first:
##         exprs[exprs == 0] <- NA
##   Imputed values are random draws from uniform in [0, LOD].
##   impute_quantile: quantile of minimum non-zero expression value for each 
##     gene to use as global LOD. 
## 
## f.impute_sample_lod(exprs) return(exprs)
##   for each sample, convert NAs to minimum non-NA value in that sample:
## 
## f.impute_unif_sample_lod(exprs, impute_quantile=0.01) return(exprs)
##   for each sample, impute random draws from:
##     uniform in (0, quantile(exprs[, sample], probs=impute_quantile))
##
## f.impute_rnorm_feature(exprs, scale.=1) return(exprs)
##   impute random draws from:
##     normal(mean=mean(exprs[feature, ]), sd=scale.*sd(exprs[feature, ])):
##
## f.impute_glm_binom(x, gran=0.0001, off=1, f_mid=median) return(x)
##     impute random draws by estimating p(missing|intensity) using a 
##   binomial glm of form cbind(n.missing, n.found) ~ intensity,
##     where intensity is assumed to be previously log transformed.
##   x: numeric matrix with feature rows and observation columns;
##   gran: granularity of predictions (smaller = less chance of dups); numeric > 0
##   off: offset for p.missing = (n.miss + off)/(n.total + off); numeric;
##   f_mid: function for aggregating intensities into a central estimate; function
##
## f.impute_loess_logit <- function(x, n=1e5, off=1, f_mid=median, 
##   span=0.25, degree=1, fam="symmetric") return(x)
##   impute random draws by estimating p(missing|intensity) using a 
##     loess fit of form log(p_missing/p_found) ~ intensity,
##       where intensity is assumed to be previously log transformed.
##   x: numeric matrix with feature rows and observation columns;
##   n: granularity of predictions (larger = less chance of dups); integer > 0
##   off: offset for p.missing = (n.miss + off)/(n.total + off); numeric >= 0;
##   f_mid: function for aggregating intensities into a central estimate; function
##   span: span for loess fit; numeric in closed interval (0, 1)
##   degree: degree for loess fit; integer in c(1, 2)
##   fam: family for loess fit; character in c("gaussian", "symmetric")
##
## f.impute_rf(x, f_imp=f.impute_loess_logit, ntree=100, mtry=NULL, 
##     aug_mult=0.33, aug_add=0, aug_steps=3)
##       return(list(x, log)); where log is data.frame w/ logging info;
##   f_imp: method to use for initial rough imputation;
##   ntree: number of trees in random forest
##   mtry: random forest mtry
##   aug_mult: multiplier for x_train = f.augment_affine(x)
##   aug_add: shift for x_train = f.augment_affine(x)
##   aug_steps: steps for x_train = f.aug_mult(x)
##
## f.impute_glmnet(x, f_imp=f.impute_loess_logit, alpha=1, nfolds=5, 
##     measure="mae", aug_mult=0.33, aug_add=0, aug_steps=3); 
##       return(list(x, log)) where tbl is data.frame w/ logging info;
##   measure: loss in c("mae", "mse")
##   nfolds: for CV tuning integer in 2:ncol(x)
##   aug_steps: for augmenting data matrix integer in 0:Inf; 
##   alpha numeric in 0:1; 
## 
## f.augment_affine <- function(x, mult=1, add=0, steps=1) return(x)
##     augments x with affine transformed versions of itself;
##   e.g., for steps=3: x = cbind(x, mult*x+add, (mult^2)*x+2*add, (mult^3)*x+3*add) 
##   mult: scale for affine transform
##   add: shift for affine transform
##   steps: number of times to apply affine transform
## 
## f.test_voom <- function(exprs, samps, frm, test_var, normalize.method="none") 
##   return(tbl)
## 
## f.test_trend <- function(exprs, samps, frm, test_var) return(tbl)


###############################################################################
## logging:

if(!exists("LOG_FILE")) {
  stop("!exists('LOG_FILE'): define variable 'LOG_FILE' before sourcing")
} else if(file.exists(LOG_FILE)) {
  stop("file.exists(LOG_FILE): LOG_FILE must not exist")
}

library("limma")
library("edgeR")

f.msg <- function(...) {
  cat(..., "\n", file=LOG_FILE, append=T)
  flush.console()
}

f.log <- function(...) {
  cat(..., "at:", format(Sys.time(), format='%Y%m%d%H%M%S'), "\n", file=LOG_FILE, append=T)
  flush.console()
}

f.err <- function(...) {
  f.log("ERROR:", ...)
  stop("Stopping", call.=F)
}


###############################################################################
## save files:

f.save_tsv <- function(dat, file_out, row.names=T, col.names=T) {
  tryCatch(
    write.table(dat, file=file_out, quote=F, sep="\t", 
      row.names=row.names, col.names=col.names),
    error=function(msg) f.err("write.table() error: writing to ", file_out, ": ", msg$message),
    warning=function(msg) f.err("write.table() warning: writing to ", file_out, ": ", msg$message)
  )
}


###############################################################################
## filtering:

## filter features based on number of expressing samples:

f.filter_features <- function(exprs, samps, feats, n_samples_min=2) {
  if(!is.matrix(exprs)) f.err("f.filter_features: !is.matrix(exprs)")

  f <- function(v) {
    i <- sum(v > 0, na.rm=T) >= n_samples_min
    i[is.na(i)] <- F
    return(i)
  }
  i <- apply(exprs, 1, f)

  f.msg("filtering", sum(!i), "features, keeping", sum(i))
  exprs <- exprs[i, ]
  feats <- feats[i, ]
  f.msg("features filtered: nrow(exprs):", nrow(exprs), "; ncol(exprs):", ncol(exprs))

  return(list(exprs=exprs, samps=samps, feats=feats))
}

f.filter_samples <- function(exprs, samps, feats, n_features_min=1000) {
  if(!is.matrix(exprs)) f.err("f.filter_samples: !is.matrix(exprs)")

  f <- function(v) {
    i <- sum(v > 0, na.rm=T) >= n_features_min
    i[is.na(i)] <- F
    return(i)
  }
  i <- apply(exprs, 2, f)

  f.msg("filtering", sum(!i), "samples, keeping", sum(i))
  exprs <- exprs[, i]
  samps <- samps[i, ]
  f.msg("samples filtered: nrow(exprs):", nrow(exprs), "; ncol(exprs):", ncol(exprs))

  return(list(exprs=exprs, samps=samps, feats=feats))
}


###############################################################################
## filter-related stats:

f.samples_per_feature <- function(exprs) {
  if(!is.matrix(exprs)) f.err("f.samples_per_feature: !is.matrix(exprs)")

  f <- function(v) {
    n <- sum(v > 0, na.rm=T)
    n[is.na(n)] <- 0      ## only if all(is.na(v))
    return(n)
  }
  n <- apply(exprs, 1, f)
  return(n)
}

f.feature_median_expression <- function(exprs) {
  if(!is.matrix(exprs)) f.err("f.feature_median_expression: !is.matrix(exprs)")
  f <- function(v) {
    m <- median(v, na.rm=T)
    m[is.na(m)] <- 0      ## only if all(is.na(v))
    return(m)
  }
  m <- apply(exprs, 1, f)
  return(m)
}

f.features_per_sample <- function(exprs) {
  if(!is.matrix(exprs)) f.err("f.features_per_sample: !is.matrix(exprs)")
  f <- function(v) {
    n <- sum(v > 0, na.rm=T)
    n[is.na(n)] <- 0      ## only if all(is.na(v))
    return(n)
  }
  n <- apply(exprs, 2, f)
  return(n)
}


###############################################################################
## normalization:

## inter-sample edger (tmm) normalization;
##   normalized cpm values using edger norm factors
##   method in c("TMM", "TMMwsp", "RLE", "upperquartile", "none")
##   returns list with elements: exprs, samps, feats

f.normalize_edger <- function(exprs, method="TMM", p=0.75) {
  if(!is.matrix(exprs)) f.err("f.normalize_edger: !is.matrix(exprs)")

  f.log("convert NA to zero")
  exprs[is.na(exprs)] <- 0

  f.log("making edgeR object")
  obj <- edgeR::DGEList(exprs)

  f.log("edgeR::calcNormFactors")
  obj <- edgeR::calcNormFactors(obj, method=method, p=p)

  f.log("making normalized expression")
  exprs <- edgeR::cpm(obj, normalized.lib.sizes=T, log=F, prior.count=1)

  f.log("convert zero back to NA")
  exprs[exprs == 0] <- NA

  return(exprs)
}

## inter-sample quantile normalization; same as limma scale normalization 
##   when norm_quantile == 0.5 (median):

f.normalize_quantile <- function(exprs, norm_quantile=0.5, multiplier=1e3) {
  if(!is.matrix(exprs)) f.err("f.normalize_quantile: !is.matrix(exprs)")
  f <- function(v) multiplier * v / quantile(v, probs=norm_quantile, na.rm=T)
  exprs <- apply(exprs, 2, f)
  return(exprs)
}

## counts per million (if multiplier is 1e6):

f.normalize_cpm <- function(exprs, multiplier=1e6) {
  if(!is.matrix(exprs)) f.err("f.normalize_cpm: !is.matrix(exprs)")
  f <- function(v) multiplier * (v / sum(v, na.rm=T))
  exprs <- apply(exprs, 2, f)
  return(exprs)
}

## variance stabilizing normalization; 
##   results on log2 scale:

f.normalize_vsn <- function(exprs) {
  if(!is.matrix(exprs)) f.err("f.normalize_vsn: !is.matrix(exprs)")
  exprs <- limma::normalizeVSN(exprs)
  return(exprs)
}

## cyclic loess normalization; 
##   method in c("fast", "affy", "pairs"); 
##   slow! fast is linear, others are quadratic:
f.normalize_loess <- function(exprs, span=0.7, method="affy") {
  if(!is.matrix(exprs)) f.err("f.normalize_loess: !is.matrix(exprs)")
  exprs <- limma::normalizeCyclicLoess(exprs, span=span, method=method)
  return(exprs)
}

## Old-school 'quantile' normalization
##    matches (all quantiles of) signal distribution across samples;
f.normalize_qquantile <- function(exprs) {
  if(!is.matrix(exprs)) f.err("f.normalize_qquantile: !is.matrix(exprs)")
  exprs <- limma::normalizeQuantiles(exprs)
  return(exprs)
}


###############################################################################
## imputation:

## Impute missing values assuming left-censored only.
##   impute_quantile: quantile of minimum non-zero expression value for each 
##     gene to use as global LOD. 
##   Returns expression matrix same as exprs but w/ missing values imputed.
##     ONLY NA VALUES CONSIDERED MISSING, so may need to first:
##       exprs[exprs == 0] <- NA
##   Imputed values are random draws from uniform in [0, LOD].

f.impute_unif_global_lod <- function(exprs, impute_quantile=0) {
  if(!is.matrix(exprs)) f.err("f.impute_unif_global_lod: !is.matrix(exprs)")
  v <- apply(exprs, 1, min, na.rm=T)
  v <- v[!is.na(v)]
  v <- v[v > 0]        ## minimum non-zero values for each gene (with non-0 vals)
  max_val <- quantile(v, probs=impute_quantile[1], na.rm=T)    ## quantile of min vals
  if(is.na(max_val)) f.err("f.impute_unif_global_lod: is.na(max_val)")
  names(max_val) <- NULL
  f.msg("f.impute_unif_global_lod: inpute_quantile: ", impute_quantile, "; max_val:", max_val)
  f <- function(v, max_val) {
    i <- is.na(v)
    if(any(i)) {
      v[i] <- runif(sum(i), 0, max_val)
    }
    v
  }
  exprs <- t(apply(exprs, 1, f, max_val))
  return(exprs)
}

## for each sample, convert NAs to minimum non-NA value in that sample:
f.impute_sample_lod <- function(exprs) {
  if(!is.matrix(exprs)) f.err("f.impute_sample_lod: !is.matrix(exprs)")
  f <- function(v) {
    v[is.na(v)] <- min(v, na.rm=T)
    v
  }
  exprs <- apply(exprs, 2, f)
  return(exprs)
}

## for each sample, impute random draws from:
##   uniform in (0, quantile(exprs[, sample], probs=impute_quantile))
f.impute_unif_sample_lod <- function(exprs, impute_quantile=0) {
  if(!is.matrix(exprs)) f.err("f.impute_unif_sample_lod: !is.matrix(exprs)")
  f <- function(v) {
    i <- is.na(v)
    if(any(i)) {
      max_val <- quantile(v, probs=impute_quantile, na.rm=T)
      v[i] <- runif(sum(i), 0, max_val)
    }
    v
  }
  exprs <- apply(exprs, 2, f)
  return(exprs)
}

## impute random draws from:
##   normal(mean=mean(exprs[feature, ]), sd=scale.*sd(exprs[feature, ])),
##     ensuring all imputed values are non-negative (redraw on negative):

f.impute_rnorm_feature <- function(exprs, scale.=1) {
  if(!is.matrix(exprs)) f.err("f.impute_rnorm_feature: !is.matrix(exprs)")
  i <- c(exprs) < 0
  i[is.na(i)] <- F
  if(any(i)) f.err("f.impute_rnorm_feature: exprs contains negative values")
  f <- function(v) {
    i <- is.na(v)
    if(any(i)) {
      m <- mean(v, na.rm=T)
      s <- sd(v, na.rm=T) * scale.
      v[i] <- rnorm(sum(i), mean=m, sd=s)
      i <- v < 0
      while(any(i)) {
        v[i] <- rnorm(sum(i), mean=m, sd=s)
        i <- v < 0
      }
    }
    return(v)
  }
  exprs <- t(apply(exprs, 1, f))
  return(exprs)
}

## impute random draws by estimating p(missing|intensity) using a 
##   binomial glm of form cbind(n.missing, n.found) ~ intensity,
##     where intensity is assumed to be previously log transformed.
##   f.impute_glm_binom(x, gran=0.0001, off=1, f_mid=median) return(x)
##     x: numeric matrix with feature rows and observation columns;
##     gran: granularity of predictions (smaller = less chance of dups); numeric > 0
##     off: offset for p.missing = (n.miss + off)/(n.total + off); numeric;
##     f_mid: function for aggregating intensities into a central estimate; function

f.impute_glm_binom <- function(x, gran=0.0001, off=1, f_mid=median) {

  if(!is.matrix(x)) stop("!is.matrix(x)")
  if(!is.numeric(x)) stop("!is.numeric(x)")

  m <- apply(x, 1, f_mid, na.rm=T)
  m[is.na(m)] <- 0
  n0 <- apply(x, 1, function(v) sum(is.na(v) | v %in% 0))  ## n.missing
  n1 <- ncol(x) - n0                                     
  p <- (n0 + off) / (ncol(x) + off)
  dat <- data.frame(n0=n0, n1=n1, p=p, m=m)

  fit <- glm(cbind(n0, n1) ~ m, data=dat, family="binomial")
  m <- seq(from=gran, to=max(c(x), na.rm=T), by=gran)
  p <- predict(fit, newdata=data.frame(m=m), type="response")
  p[is.na(p)] <- min(p, na.rm=T)
  x[is.na(c(x))] <- sample(m, sum(is.na(c(x))), replace=T, prob=p)

  return(x)
}

## impute random draws by estimating p(missing|intensity) using a 
##   loess fit of form log(p_missing/p_found) ~ intensity,
##     where intensity is assumed to be previously log transformed.
##   f.impute_loess_logit(x, gran=0.0001, off=1, f_mid=median, 
##       span=0.25, degree=1, fam="symmetric") return(x)
##     x: numeric matrix with feature rows and observation columns;
##     gran: granularity of predictions (smaller: less chance of duplicates); numeric > 0
##     off: offset for probability.missing = (n.miss + off)/(n.total + off); numeric;
##     f_mid: function for aggregating intensities into a central estimate; function
##     span: span for loess fit; numeric in closed interval (0, 1)
##     degree: degree for loess fit; integer in c(1, 2)
##     fam: family for loess fit; character in c("gaussian", "symmetric")

f.impute_loess_logit <- function(x, n=1e5, off=1, f_mid=median, 
    span=0.25, degree=1, fam="symmetric") {

  if(!is.matrix(x)) stop("!is.matrix(x)")
  if(!is.numeric(x)) stop("!is.numeric(x)")

  m <- apply(x, 1, f_mid, na.rm=T)
  m[is.na(m)] <- 0
  n0 <- apply(x, 1, function(v) sum(is.na(v) | v %in% 0))  ## n.missing
  n1 <- ncol(x) - n0                                       ## n.found
  p <- (n0 + off) / (ncol(x) + off)                        ## p.missing
  dat <- data.frame(n0=n0, n1=n1, p=p, m=m)

  m_new <- seq(from=1/n, to=max(c(x), na.rm=T), length=n)  ## newdata
  fit <- loess(log(p/(1-p)) ~ m, data=dat, span=span, degree=degree, family=fam)
  p_hat <- predict(fit, newdata=data.frame(m=m_new))       ## on logit scale
  p_hat[is.na(p_hat)] <- min(p_hat[p_hat > 0], na.rm=T)
  p_hat = exp(p_hat) / (1 + exp(p_hat))                    ## inverse logit
  x[is.na(c(x))] <- sample(m, sum(is.na(c(x))), replace=T, prob=p)

  return(x)
}

## helper for augmenting training data with linearly transformed
##   (multiplied by mult) versions of itself;
##   mult: numeric multiplier > 0;
##   steps: number of augmentatons, each one x' = x * mult^step:

f.augment_affine <- function(x, mult=1, add=0, steps=1) {

  xs <- list()
  xs[[length(xs) + 1]] <- x

  while(steps > 0) {
    x <- (x * mult) + add
    xs[[length(xs) + 1]] <- x
    steps = steps - 1
  }
  x <- do.call(cbind, xs)
  return(x)
}

## f_imp: method to use for initial rough imputation;
## ntree: number of trees in random forest
## mtry: random forest mtry
## aug_mult: multiplier for x_train = f.aug_mult(x)
## aug_steps: steps for x_train = f.aug_mult(x)

f.impute_rf <- function(x, f_imp=f.impute_loess_logit, ntree=100, mtry=NULL, aug_mult=0.33, aug_add=0, aug_steps=3) {

  if(!is.matrix(x)) stop("!is.matrix(x)")
  if(!is.numeric(x)) stop("!is.numeric(x)")

  n_miss <- apply(x, 1, function(v) sum(is.na(v) | v %in% 0))
  o <- order(n_miss, decreasing=F)
  tbl <- NULL
  for(idx_feat in o) {
    if(n_miss[idx_feat] < 1) next
    cat("processing", rownames(x)[idx_feat], "\n"); flush.console()

    x_train <- f_imp(x)
    x_train <- 2^x_train - 1
    x_train <- f.augment_affine(x_train, mult=aug_mult, add=aug_add, steps=aug_steps)
    x_train <- log2(x_train + 1)

    y <- x[idx_feat, , drop=T]
    i_miss <- is.na(y) | y %in% 0
    col_names_miss <- colnames(x)[i_miss]
    i_miss_train <- colnames(x_train) %in% col_names_miss
    x_train_i <- x_train[-idx_feat, !i_miss_train, drop=F]
    y_train_i <- x_train[idx_feat, !i_miss_train, drop=T]
    if(is.null(mtry)) { mtry0 <- round(sqrt(nrow(x_train_i))) } else { mtry0 <- mtry }
    fit <- randomForest::randomForest(x=t(x_train_i), y=y_train_i, mtry=mtry0, 
      ntree=ntree, importance=F, replace=T)

    x_miss <- x_train[-idx_feat, which(i_miss), drop=F]
    y_miss <- predict(fit, newdata=t(x_miss), type="response")
    if(any(y_miss < 0)) stop("ERROR: y_miss < 0; y_miss:", y_miss)
    x[idx_feat, i_miss] <- y_miss

    tm_stmp <- format(Sys.time(), format='%Y%m%d%H%M%S')
    tbl_i <- data.frame(time=tm_stmp, feat=rownames(x)[idx_feat], mtry=mtry0, n_miss=sum(i_miss), 
      mean0=mean(y[!i_miss]), mean1=mean(y_miss), pvar0=fit$rsq[1], pvar1=fit$rsq[ntree], 
      mse0=fit$mse[1], mse1=fit$mse[ntree]
    )
    print(tbl_i); flush.console()
    tbl <- rbind(tbl, tbl_i)
  }

  return(list(x=x, log=tbl))
}

## measure in c("mae", "mse"); nfolds integer in 2:ncol(x); aug_steps integer in 0:Inf; 
##   alpha numeric in 0:1; return(list(x=x, log=tbl))

f.impute_glmnet <- function(x, f_imp=f.impute_loess_logit, alpha=1, nfolds=5, 
    measure="mae", aug_mult=0.33, aug_add=0, aug_steps=3) {

  if(!is.matrix(x)) stop("!is.matrix(x)")
  if(!is.numeric(x)) stop("!is.numeric(x)")

  n_miss <- apply(x, 1, function(v) sum(is.na(v) | v %in% 0))
  o <- order(n_miss, decreasing=F)
  tbl <- NULL
  for(idx_feat in o) {
    if(n_miss[idx_feat] < 1) next
    cat("processing", rownames(x)[idx_feat], "\n"); flush.console()

    x_train <- f_imp(x)
    x_train <- 2^x_train - 1
    x_train <- f.augment_affine(x_train, mult=aug_mult, add=aug_add, steps=aug_steps)
    x_train <- log2(x_train + 1)

    y <- x[idx_feat, , drop=T]
    i_miss <- is.na(y) | y %in% 0
    col_names_miss <- colnames(x)[i_miss]
    i_miss_train <- colnames(x_train) %in% col_names_miss
    x_train_i <- x_train[-idx_feat, !i_miss_train, drop=F]
    y_train_i <- x_train[idx_feat, !i_miss_train, drop=T]

    fit <- glmnet::cv.glmnet(x=t(x_train_i), y=y_train_i, family="gaussian", alpha=alpha, 
      type.measure=measure, nfolds=nfolds, parallel=F)

    x_miss <- x_train[-idx_feat, which(i_miss), drop=F]
    y_miss <- predict(fit, newx=t(x_miss), s=fit$lambda.1se, type="response")
    if(any(y_miss < 0)) stop("ERROR: y_miss < 0; y_miss:", y_miss)
    x[idx_feat, i_miss] <- y_miss

    idx_lambda <- which(fit$lambda %in% fit$lambda.1se)[1] 
    tm_stmp <- format(Sys.time(), format='%Y%m%d%H%M%S')
    tbl_i <- data.frame(time=tm_stmp, feat=rownames(x)[idx_feat], alpha=alpha, nfolds=nfolds, 
      n_miss=sum(i_miss), mean0=mean(y[!i_miss]), mean1=mean(y_miss), 
      cvm0=max(fit$cvm), cvm1=fit$cvm[idx_lambda], cvup0=max(fit$cvup), cvup1=fit$cvup[idx_lambda]
    )
    print(tbl_i); flush.console()
    tbl <- rbind(tbl, tbl_i)
  }

  return(list(x=x, log=tbl))
}


###############################################################################
## hypothesis testing:
 
## can use quantile normalization here:
f.test_voom <- function(exprs, samps, frm, test_term, normalize.method="none") {
  if(!is.matrix(exprs)) f.err("f.test_voom: !is.matrix(exprs)")
  i <- apply(exprs, 1, function(v) any(is.na(v)))
  exprs <- exprs[!i, ]
  design <- model.matrix(frm, data=samps)
  obj <- limma::voom(exprs, design, plot=F, normalize.method=normalize.method)
  fit <- limma::lmFit(obj, design)
  fit <- limma::eBayes(fit, trend=F)
  cols <- colnames(model.matrix(as.formula(paste("~0 + ", test_term)), data=samps))
  cols <- cols[cols %in% colnames(design)]
  i <- colnames(design) %in% cols
  tbl <- limma::topTable(fit, coef=which(i), number=Inf)
  f.msg("tested", nrow(exprs), "genes")
  f.msg("found", sum(tbl$adj.P.Val < 0.05, na.rm=T), "hits")
  return(tbl)
}

f.test_trend <- function(exprs, samps, frm, test_term) {
  if(!is.matrix(exprs)) f.err("f.test_trend: !is.matrix(exprs)")
  design <- model.matrix(frm, data=samps)
  fit <- limma::lmFit(exprs, design)
  fit <- limma::eBayes(fit, trend=T)
  cols <- colnames(model.matrix(as.formula(paste("~0 + ", test_term)), data=samps))
  cols <- cols[cols %in% colnames(design)]
  i <- colnames(design) %in% cols
  tbl <- limma::topTable(fit, coef=which(i), number=Inf)
  f.msg("tested", nrow(exprs), "genes")
  f.msg("found", sum(tbl$adj.P.Val < 0.05, na.rm=T), "hits")
  return(tbl)
}

