#' Normalize expression using \code{edgeR}
#' @description
#'   Normalize expression using functionality from \code{edgeR} package.
#' @details Inter-observation normalization. Uses \code{edgeR::calcNormFactors()}. 
#'   Returned values on a counts-per-million scale.
#' @param state List with elements formatted like the list returned by \code{f.read_data()}:
#'   \tabular{ll}{
#'     \code{expression} \cr \tab Numeric matrix with non-negative expression values. \cr
#'     \code{features}   \cr \tab A data.frame with feature meta-data for rows of expression. \cr
#'     \code{samples}    \cr \tab A data.frame with observation meta-data for columns of expression. \cr
#'   } 
#' @param config List with configuration values. Uses the following keys:
#'   \tabular{ll}{
#'     \code{norm_method}     \cr \tab Normalization method in \code{c("TMM", "TMMwsp", "RLE", "upperquartile")}. \cr
#'     \code{norm_quantile}   \cr \tab Quantile to use for \code{norm_method \%in\% c("quantile", "upperquartile")}. \cr
#'   }
#' @param method Character in set
#'   \code{c("TMM", "TMMwsp", "RLE", "upperquartile", "none")}. Can be set with
#'   \code{config$norm_method}.
#' @param p Numeric in closed interval \code{[0, 1]} specifying quantile to use for method 
#'   \code{upperquartile}. Can be set with \code{config$norm_quantile}.
#' @return A list (the processed state) with the following elements:
#'   \tabular{ll}{
#'     \code{expression} \cr \tab Numeric matrix with non-negative expression values. \cr
#'     \code{features}   \cr \tab A data.frame with feature meta-data for rows of expression. \cr
#'     \code{samples}    \cr \tab A data.frame with observation meta-data for columns of expression. \cr
#'   } 
#' @examples
#' set.seed(101)
#' exprs <- h0testr::f.sim1(n_obs=6, n_feats=12)$mat
#' feats <- data.frame(feature_id=rownames(exprs))
#' samps <- data.frame(observation_id=colnames(exprs))
#' state <- list(expression=exprs, features=feats, samples=samps)
#'
#' config <- list(norm_method="upperquartile", norm_quantile=0.75)
#' state2 <- h0testr::f.normalize_edger(state, config)
#' print(state$expression)
#' print(state2$expression)
#'
#' ## norm_quantile ignored for norm_method="TMM":
#' config <- list(norm_method="TMM", norm_quantile=0)
#' state2 <- h0testr::f.normalize_edger(state, config)
#' print(state$expression)
#' print(state2$expression)

f.normalize_edger <- function(state, config, method=NULL, p=NULL) {

  f.check_config(config)

  if(!is.matrix(state$expression)) {
    f.err("f.normalize_edger: !is.matrix(state$expression)", 
      config=config)
  }
  if(is.null(method)) method <- config$norm_method
  if(is.null(method)) {
    f.err("f.normalize_edger: method and config$norm_method both unset", 
      config=config)
  }
  
  if(method %in% "upperquartile") {
    if(is.null(p)) p <- config$norm_quantile
    if(is.null(p) && method %in% "upperquartile") {
      f.err("f.normalize_edger: p and config$norm_quantile both unset", 
        config=config)
    }
    if(!is.numeric(p)) {
        f.err("f.normalize_edger: p and config$norm_quantile need to be numeric", 
          config=config)
    }
  }

  f.log("convert NA to zero", config=config)
  state$expression[is.na(state$expression)] <- 0

  f.log("making edgeR object", config=config)
  obj <- edgeR::DGEList(state$expression)

  f.log("edgeR::calcNormFactors", config=config)
  obj <- edgeR::calcNormFactors(obj, method=method, p=p)

  f.log("making normalized expression", config=config)
  state$expression <- edgeR::cpm(obj, normalized.lib.sizes=T, 
    log=F, prior.count=1)

  f.log("convert zero back to NA", config=config)
  state$expression[state$expression == 0] <- NA

  return(state)
}

#' Normalize expression using selected quantile as scaling factor.
#' @description
#'   Normalize expression using variant of quantile normalization which 
#'     excludes missing values.
#' @details 
#'   Inter-observation normalization. Setting \code{p <- 0.5} is equivalent 
#'     to median scaling with median calculated after exclusion of missing 
#'     values. Similarly, setting \code{p < 0.75} is upperquartile normalization 
#'     ignoring missing values. 
#' @param state List with elements formatted like the list returned by \code{f.read_data()}:
#'   \tabular{ll}{
#'     \code{expression} \cr \tab Numeric matrix with non-negative expression values. \cr
#'     \code{features}   \cr \tab A data.frame with feature meta-data for rows of expression. \cr
#'     \code{samples}    \cr \tab A data.frame with observation meta-data for columns of expression. \cr
#'   } 
#' @param config List with configuration values. Uses the following keys:
#'   \tabular{ll}{
#'     \code{norm_quantile}   \cr \tab Quantile to use where \code{norm_method \%in\% c("quantile", "upperquartile")}. \cr
#'   }
#' @param norm_quantile Numeric in closed interval \code{[0, 1]} specifying quantile to use.
#' @param multiplier Numeric used to scale returned values after 
#'   dividing by selected quantile. 
#' @return A list (the processed state) with the following elements:
#'   \tabular{ll}{
#'     \code{expression} \cr \tab Numeric matrix with non-negative expression values. \cr
#'     \code{features}   \cr \tab A data.frame with feature meta-data for rows of expression. \cr
#'     \code{samples}    \cr \tab A data.frame with observation meta-data for columns of expression. \cr
#'   } 
#' @examples
#' set.seed(101)
#' exprs <- h0testr::f.sim1(n_obs=6, n_feats=12)$mat
#' feats <- data.frame(feature_id=rownames(exprs))
#' samps <- data.frame(observation_id=colnames(exprs))
#' state <- list(expression=exprs, features=feats, samples=samps)
#' config <- list(norm_quantile=0.75)
#' state2 <- h0testr::f.normalize_quantile(state, config)
#' print(state$expression)
#' print(state2$expression)
#' apply(state$expression, 2, quantile, probs=c(0.5, 0.75, 0.9), na.rm=TRUE)
#' apply(state2$expression, 2, quantile, probs=c(0.5, 0.75, 0.9), na.rm=TRUE)

f.normalize_quantile <- function(state, config, norm_quantile=NULL, multiplier=1e3) {

  f.check_config(config)

  if(!is.matrix(state$expression)) {
    f.err("f.normalize_quantile: !is.matrix(state$expression)", config=config)
  }
  if(is.null(norm_quantile)) norm_quantile <- config$norm_quantile
  if(is.null(norm_quantile)) {
    f.err("f.normalize_quantile: is.null(norm_quantile)", config=config)
  }
  if(!is.numeric(norm_quantile)) {
    f.err("f.normalize_quantile: !is.numeric(norm_quantile)", config=config)
  }
  
  f <- function(v) {
    multiplier * v / stats::quantile(v, probs=norm_quantile, na.rm=T)
  }
  state$expression <- apply(state$expression, 2, f)
  
  return(state)
}

#' Normalize expression using total expression in each sample.
#' @description
#'   Normalize expression using variant of CPM normalization which excludes missing values.
#' @details 
#'   Inter-observation normalization, based on total expression
#'     in each observation. Makes total expression (excluding missing values)
#'     equal in each observation.
#' @param state List with elements formatted like the list returned by \code{f.read_data()}:
#'   \tabular{ll}{
#'     \code{expression} \cr \tab Numeric matrix with non-negative expression values. \cr
#'     \code{features}   \cr \tab A data.frame with feature meta-data for rows of expression. \cr
#'     \code{samples}    \cr \tab A data.frame with observation meta-data for columns of expression. \cr
#'   } 
#' @param config List with configuration values. Does not use any keys, so can pass empty list.
#' @param multiplier Numeric greater than zero used to scale returned values 
#'   after dividing by total counts in observation. For example, the default 
#'   \code{multiplier=1e6} yields normalized expression as CPM (counts per million). 
#' @return A list (the processed state) with the following elements:
#'   \tabular{ll}{
#'     \code{expression} \cr \tab Numeric matrix with non-negative expression values. \cr
#'     \code{features}   \cr \tab A data.frame with feature meta-data for rows of expression. \cr
#'     \code{samples}    \cr \tab A data.frame with observation meta-data for columns of expression. \cr
#'   } 
#' @examples
#' set.seed(101)
#' exprs <- h0testr::f.sim1(n_obs=6, n_feats=12)$mat
#' feats <- data.frame(feature_id=rownames(exprs))
#' samps <- data.frame(observation_id=colnames(exprs))
#' state <- list(expression=exprs, features=feats, samples=samps)
#' config <- list()
#' state2 <- h0testr::f.normalize_cpm(state, config)
#' print(state$expression)
#' print(state2$expression)
#' apply(state$expression, 2, function(v) sum(v, na.rm=TRUE))
#' apply(state2$expression, 2, function(v) sum(v, na.rm=TRUE))

f.normalize_cpm <- function(state, config, multiplier=1e6) {

  if(!is.matrix(state$expression)) {
    f.err("f.normalize_cpm: !is.matrix(state$expression)", config=config)
  }
  
  f <- function(v) multiplier * (v / sum(v, na.rm=T))
  state$expression <- apply(state$expression, 2, f)
  return(state)
}

#' Variance stabilizing transformation.
#' @description
#'   Normalize expression using variance stabilizing transformation.
#' @details 
#'   Inter-observation normalization. Uses \code{limma::normalizeVSN()}. 
#'     Unlike most other normalization methods, results are returned on a 
#'     log2-like scale.
#' @param state List with elements formatted like the list returned by \code{f.read_data()}:
#'   \tabular{ll}{
#'     \code{expression} \cr \tab Numeric matrix with non-negative expression values. \cr
#'     \code{features}   \cr \tab A data.frame with feature meta-data for rows of expression. \cr
#'     \code{samples}    \cr \tab A data.frame with observation meta-data for columns of expression. \cr
#'   } 
#' @param config List with configuration values. Does not use any keys, so can pass empty list.
#' @return A list (the processed state) with the following elements:
#'   \tabular{ll}{
#'     \code{expression} \cr \tab Numeric matrix with non-negative expression values. \cr
#'     \code{features}   \cr \tab A data.frame with feature meta-data for rows of expression. \cr
#'     \code{samples}    \cr \tab A data.frame with observation meta-data for columns of expression. \cr
#'   } 
#' @examples
#' set.seed(101)
#' exprs <- h0testr::f.sim1(n_obs=6, n_feats=50, mnar_c0=-Inf)$mat
#' feats <- data.frame(feature_id=rownames(exprs))
#' samps <- data.frame(observation_id=colnames(exprs))
#' state <- list(expression=exprs, features=feats, samples=samps)
#' config <- list()
#' state2 <- h0testr::f.normalize_vsn(state, config)
#' head(state$expression)
#' head(state2$expression)
#' summary(apply(state$expression, 2, sd, na.rm=TRUE))
#' summary(apply(log2(state$expression+1), 2, sd, na.rm=TRUE))
#' summary(apply(state2$expression, 2, sd, na.rm=TRUE))
#' summary(apply(state$expression, 1, sd, na.rm=TRUE))
#' summary(apply(log2(state$expression+1), 1, sd, na.rm=TRUE))
#' summary(apply(state2$expression, 1, sd, na.rm=TRUE))

f.normalize_vsn <- function(state, config) {

  if(!is.matrix(state$expression)) {
    f.err("f.normalize_vsn: !is.matrix(state$expression)", config=config)
  }
  
  state$expression <- limma::normalizeVSN(state$expression)
  return(state)
}

#' Cyclic loess normalization.
#' @description
#'   Normalize expression using the cyclic-loess algorithm.
#' @details 
#'   Inter-observation normalization using cyclic-loess results
#'     in similar signal distributions across all samples, similar to 
#'     \code{f.normalize_qquantile()}. This is a slow method, especially if
#'     \code{method \%in\% c("affy", "pairs")}, which scale quadratically. Calls
#'     \code{limma::normalizeCyclicLoess()} under the hood.
#' @param state List with elements formatted like the list returned by `f.read_data()`:
#'   \tabular{ll}{
#'     \code{expression} \cr \tab Numeric matrix with non-negative expression values. \cr
#'     \code{features}   \cr \tab A data.frame with feature meta-data for rows of expression. \cr
#'     \code{samples}    \cr \tab A data.frame with observation meta-data for columns of expression. \cr
#'   } 
#' @param config List with configuration values. Does not use any keys, so can pass empty list.
#' @param span Numeric between 0 and 1 specifying span for loess fit.
#'   Higher numbers result in smoother (less localized) fit.
#' @param method Character in \code{c("fast", "affy", "pairs")}.
#' @return A list (the processed state) with the following elements:
#'   \tabular{ll}{
#'     \code{expression} \cr \tab Numeric matrix with non-negative expression values. \cr
#'     \code{features}   \cr \tab A data.frame with feature meta-data for rows of expression. \cr
#'     \code{samples}    \cr \tab A data.frame with observation meta-data for columns of expression. \cr
#'   } 
#' @examples
#' set.seed(101)
#' exprs <- h0testr::f.sim1(n_obs=6, n_feats=30, mnar_c0=-Inf)$mat
#' feats <- data.frame(feature_id=rownames(exprs))
#' samps <- data.frame(observation_id=colnames(exprs))
#' state <- list(expression=exprs, features=feats, samples=samps)
#' config <- list()
#' state2 <- h0testr::f.normalize_loess(state, config)
#' head(state$expression)
#' head(state2$expression)
#' summary(apply(state$expression, 2, sd, na.rm=TRUE))
#' summary(apply(state2$expression, 2, sd, na.rm=TRUE))
#' summary(apply(state$expression, 1, sd, na.rm=TRUE))
#' summary(apply(state2$expression, 1, sd, na.rm=TRUE))

f.normalize_loess <- function(state, config, span=0.7, method="affy") {

  if(!is.matrix(state$expression)) {
    f.err("f.normalize_loess: !is.matrix(state$expression)", config=config)
  }
  
  state$expression <- limma::normalizeCyclicLoess(state$expression, 
    span=span, method=method)
    
  return(state)
}

#' Old-school 'quantile' normalization.
#' @description
#'   Normalize expression values using the 'quantile normalization' algorithm.
#' @details 
#'   Inter-observation normalization resulting in nearly identical
#'     signal distributions across all samples, so all quantiles in \code{0:1} match
#'     across all samples. Calls \code{limma::normalizeQuantiles()} under the hood.
#' @param state List with elements formatted like the list returned by \code{f.read_data()}:
#'   \tabular{ll}{
#'     \code{expression} \cr \tab Numeric matrix with non-negative expression values. \cr
#'     \code{features}   \cr \tab A data.frame with feature meta-data for rows of expression. \cr
#'     \code{samples}    \cr \tab A data.frame with observation meta-data for columns of expression. \cr
#'   } 
#' @param config List with configuration values. Does not use any keys, so can pass empty list.
#' @return A list (the processed state) with the following elements:
#'   \tabular{ll}{
#'     \code{expression} \cr \tab Numeric matrix with non-negative expression values. \cr
#'     \code{features}   \cr \tab A data.frame with feature meta-data for rows of expression. \cr
#'     \code{samples}    \cr \tab A data.frame with observation meta-data for columns of expression. \cr
#'   } 
#' @examples
#' set.seed(101)
#' exprs <- h0testr::f.sim1(n_obs=6, n_feats=12, mnar_c0=-Inf)$mat
#' feats <- data.frame(feature_id=rownames(exprs))
#' samps <- data.frame(observation_id=colnames(exprs))
#' state <- list(expression=exprs, features=feats, samples=samps)
#' config <- list()
#' state2 <- h0testr::f.normalize_qquantile(state, config)
#' print(state$expression)
#' print(state2$expression)
#' 
#' summary(apply(state$expression, 2, sd, na.rm=TRUE))
#' summary(apply(state2$expression, 2, sd, na.rm=TRUE))
#'
#' ## afterwards, all quantiles line up:
#' apply(state$expression, 2, quantile, probs=c(0.1, 0.25, 0.5, 0.75, 0.9), na.rm=TRUE)
#' apply(state2$expression, 2, quantile, probs=c(0.1, 0.25, 0.5, 0.75, 0.9), na.rm=TRUE)

f.normalize_qquantile <- function(state, config) {

  if(!is.matrix(state$expression)) {
    f.err("f.normalize_qquantile: !is.matrix(state$expression)", config=config)
  }
  
  state$expression <- limma::normalizeQuantiles(state$expression)
  
  return(state)
}

#' Inter-sample normalization
#' @description
#'   Normalize expression data to reduce effects of technical differences 
#'     between samples.
#' @details 
#'   Inter-observation normalization using any of the methods available in this package.
#'     See individual methods for more details.
#' @param state List with elements formatted like the list returned by \code{f.read_data()}:
#'   \tabular{ll}{
#'     \code{expression} \cr \tab Numeric matrix with non-negative expression values. \cr
#'     \code{features}   \cr \tab A data.frame with feature meta-data for rows of expression. \cr
#'     \code{samples}    \cr \tab A data.frame with observation meta-data for columns of expression. \cr
#'   } 
#' @param config List with configuration values. Uses the following keys:
#'   \tabular{ll}{
#'     \code{norm_method}   \cr \tab Character scalar in \code{c("vsn","cpm","quantile","qquantile","TMM","TMMwsp","RLE","upperquartile")}. \cr
#'     \code{norm_quantile} \cr \tab Quantile (numeric between 0 and 1) for \code{norm_method \%in\% c("quantile", "upperquartile")}. \cr
#'     \code{feat_col}      \cr \tab Column of \code{state$features} matching \code{rownames(state$expression)}. \cr
#'     \code{obs_col}       \cr \tab Column in \code{state$samples} matching \code{colnames(state$expression)}. \cr
#'   }
#' @return A list (the processed state) with the following elements:
#'   \tabular{ll}{
#'     \code{expression} \cr \tab Numeric matrix with non-negative expression values. \cr
#'     \code{features}   \cr \tab A data.frame with feature meta-data for rows of expression. \cr
#'     \code{samples}    \cr \tab A data.frame with observation meta-data for columns of expression. \cr
#'   } 
#' @examples
#' set.seed(101)
#' exprs <- h0testr::f.sim1(n_obs=6, n_feats=12)$mat
#' exprs[, 4:6] <- exprs[, 4:6] * 2
#' feats <- data.frame(feature_id=rownames(exprs))
#' samps <- data.frame(observation_id=colnames(exprs))
#' state <- list(expression=exprs, features=feats, samples=samps)
#' 
#' config <- list(save_state=FALSE, norm_method="RLE", 
#'   feat_col="feature_id", obs_col="observation_id")
#' out <- h0testr::f.normalize(state, config)
#' summary(apply(state$expression, 2, sd, na.rm=TRUE))
#' summary(apply(out$state$expression, 2, sd, na.rm=TRUE))
#' summary(apply(state$expression, 1, sd, na.rm=TRUE))
#' summary(apply(out$state$expression, 1, sd, na.rm=TRUE))
#' 
#' config$norm_method <- "quantile"
#' config$norm_quantile <- 0.5
#' out <- h0testr::f.normalize(state, config)
#' summary(apply(state$expression, 2, sd, na.rm=TRUE))
#' summary(apply(out$state$expression, 2, sd, na.rm=TRUE))
#' summary(apply(state$expression, 1, sd, na.rm=TRUE))
#' summary(apply(out$state$expression, 1, sd, na.rm=TRUE))

f.normalize <- function(state, config) {

  f.check_config(config)
  f.check_state(state, config)
  
  if(config$norm_method %in% c("TMM", "TMMwsp", "RLE", "upperquartile")) {
    state <- f.normalize_edger(state, config)
  } else if(config$norm_method %in% "quantile") {
    state <- f.normalize_quantile(state, config)
  } else if(config$norm_method %in% "cpm") {
    state <- f.normalize_cpm(state, config)
  } else if(config$norm_method %in% "vsn") {
    state <- f.normalize_vsn(state, config)
  } else if(config$norm_method %in% "qquantile") {
    state <- f.normalize_qquantile(state, config)
  } else if(config$norm_method %in% "log2") {
    f.msg("config$norm_method %in% 'log2'", config=config)
  } else if(config$norm_method %in% "none") {
    f.msg("skipping normalization: config$norm_method %in% 'none'", 
      config=config)
  } else {
    f.err("f.normalize: unexpected config$norm_method:", config$norm_method, 
      config=config)
  }
  
  if(!(config$norm_method %in% c("vsn", "none"))) {
    f.log("transforming data", config=config)
    state$expression <- log2(state$expression + 1)
  } 

  f.check_state(state, config)
  f.report_state(state, config)
  i <- config$run_order %in% "normalize"
  prfx <- paste0(which(i)[1] + 2, ".normalized")
  f.save_state(state, config, prefix=prfx)
  
  return(list(state=state, config=config))
}


