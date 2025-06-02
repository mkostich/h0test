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
## f.samples_per_feature(exprs)
##   return(n)
##
## f.feature_median_expression(exprs)
##   return(m)
##
## f.features_per_sample(exprs)
##   return(n)
##
## Inter-sample edger (tmm) normalization;
##   normalized cpm values using edger norm factors
##   method in c("TMM", "TMMwsp", "RLE", "upperquartile", "none")
## f.normalize_edger(exprs, method="TMM") 
##   return(exprs)
## 
## Inter-sample quantile normalization; same as limma scale normalization 
##   when norm_quantile == 0.5 (median); multiplier like with cpm:
## f.normalize_quantile <- function(exprs, norm_quantile=0.5, multiplier=1e6)
##   return(exprs)
## 
## Counts per million (if multiplier is 1e6):
## f.normalize_cpm(exprs, multiplier=1e6)
##   return(exprs)
## 
## Variance stabilizing normalization:
##   results on log2 scale:
## f.normalize_vsn(exprs)
##   return(exprs)
## 
## Cyclic loess normalization; 
##   method in c("fast", "affy", "pairs"); 
##   slow; fast is linear, others are quadratic:
## f.normalize_loess <- function(exprs, span=0.7, method="affy")
##   return(exprs)
## 
## Impute missing values assuming left-censored only.
##   impute_quantile: quantile of minimum non-zero expression value for each 
##     gene to use as global LOD. 
##   Returns expression matrix same as exprs but w/ missing values imputed.
##     ONLY NA VALUES CONSIDERED MISSING, so may need to first:
##       exprs[exprs == 0] <- NA
##   Imputed values are random draws from uniform in [0, LOD].
## f.impute_unif_global_lod(exprs, impute_quantile=0.01)
##   return(exprs)
## 
## for each sample, convert NAs to minimum non-NA value in that sample:
## f.impute_sample_lod(exprs)
##   return(exprs)
## 
## for each sample, impute random draws from:
##   uniform in (0, quantile(exprs[, sample], probs=impute_quantile))
## f.impute_unif_sample_lod(exprs, impute_quantile=0.01)
##   return(exprs)
## 
## impute random draws from:
##   normal(mean=mean(exprs[feature, ]), sd=scale.*sd(exprs[feature, ])):
## f.impute_rnorm_feature(exprs, scale.=1)
##   return(exprs)
## 
## f.test_voom <- function(exprs, samps, frm, test_var, normalize.method="none")
##   return(tbl)
## 
## f.test_trend <- function(exprs, samps, frm, test_var)
##   return(tbl)


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

## impute random draws by ...:




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

