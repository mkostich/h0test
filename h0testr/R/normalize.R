#' Normalize expression using \code{edgeR}
#' @description
#'   Normalize expression using functionality from \code{edgeR} package.
#' @details Inter-observation normalization. Uses \code{edgeR::calcNormFactors()}. 
#'   Returned values on a counts-per-million scale.
#'   See documentation for \code{h0testr::new_config()} 
#'     for more detailed description of configuration parameters. 
#' @param state List with elements formatted like the list returned by \code{read_data()}:
#'   \tabular{ll}{
#'     \code{expression} \cr \tab Numeric matrix with non-negative expression values. \cr
#'     \code{features}   \cr \tab A data.frame with feature meta-data for rows of expression. \cr
#'     \code{samples}    \cr \tab A data.frame with observation meta-data for columns of expression. \cr
#'   } 
#' @param config List with configuration values. None required. Uses the following keys:
#'   \tabular{ll}{
#'     \code{normalization_method}     \cr \tab Used if \code{is.null(method)}. \cr
#'     \code{normalization_quantile}   \cr \tab Used if \code{is.null(normalization_quantile)}. \cr
#'   }
#' @param method Character in set
#'   \code{c("RLE", "upperquartile", "TMM", "TMMwsp", "none")}. Can be set with
#'   \code{config$normalization_method}.
#' @param normalization_quantile Numeric where \code{0 <= normalization_quantile <= 1.0}; quantile for method 
#'   \code{upperquartile}. Can set with \code{config$normalization_quantile}.
#' @return A list (the processed state) with the following elements:
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
#'
#' ## example configured via parameters:
#' config <- list()
#' state2 <- h0testr::normalize_edger(state, config, 
#'   method="upperquartile", normalization_quantile=0.75)
#' print(state$expression)
#' print(state2$expression)
#'
#' ## example configured with settings in config:
#' config <- list(normalization_method="upperquartile", normalization_quantile=0.75)
#' state2 <- h0testr::normalize_edger(state, config)
#' print(state$expression)
#' print(state2$expression)
#'
#' ## normalizaiton_quantile not needed for e.g. method="RLE":
#' state2 <- h0testr::normalize_edger(state, config, method="RLE")
#' print(state$expression)
#' print(state2$expression)

normalize_edger <- function(state, config, method=NULL, normalization_quantile=NULL) {

  check_config(config)

  if(!is.matrix(state$expression)) {
    f.err("normalize_edger: !is.matrix(state$expression)", "\n",
      "class(state$expression):", class(state$expression), config=config)
  }
  if(is.null(method)) method <- config$normalization_method
  if(is.null(method)) method <- "RLE"
  allowed <- c("RLE", "upperquartile", "TMM", "TMMwsp", "none")
  if(!(method %in% allowed)) {
    f.err("normalize_edger: !(method %in% allowed); method:", method, "\n",
      "allowed:", allowed, config=config)
  }
  
  if(is.null(normalization_quantile)) normalization_quantile <- config$normalization_quantile
  if(is.null(normalization_quantile)) normalization_quantile <- 0.75
  if(!is.numeric(normalization_quantile)) {
      f.err("normalize_edger: !is.numeric(normalization_quantile); normalization_quantile:", 
        normalization_quantile, config=config)
  }
  if(normalization_quantile < 0 || normalization_quantile > 1) {
    f.err(
      "normalize_edger:", 
      "normalization_quantile < 0 || normalization_quantile > 1\n",
      "normalization_quantile:", normalization_quantile, config=config
    )
  }
  
  f.log("convert NA to zero", config=config)
  state$expression[is.na(state$expression)] <- 0
  
  f.log("making edgeR object", config=config)
  obj <- edgeR::DGEList(state$expression)
  
  f.log("edgeR::calcNormFactors", config=config)
  obj <- edgeR::calcNormFactors(obj, method=method, p=normalization_quantile)
  
  f.log("making normalized expression", config=config)
  state$expression <- edgeR::cpm(obj, normalized.lib.sizes=T, 
    log=F, prior.count=1)
  
  f.log("convert zero back to NA", config=config)
  state$expression[state$expression == 0] <- NA
  
  return(state)
}

#' Normalize expression using quantile of expression intensity
#' @description
#'   Normalize expression using variant of quantile normalization which 
#'     excludes missing values.
#' @details 
#'   Inter-observation normalization. Setting \code{p <- 0.5} is equivalent 
#'     to median scaling with median calculated after exclusion of missing 
#'     values. Similarly, setting \code{p < 0.75} is upperquartile normalization 
#'     ignoring missing values. 
#'   See documentation for \code{h0testr::new_config()} 
#'     for more detailed description of configuration parameters. 
#' @param state List with elements formatted like the list returned by \code{read_data()}:
#'   \tabular{ll}{
#'     \code{expression} \cr \tab Numeric matrix with non-negative expression values. \cr
#'     \code{features}   \cr \tab A data.frame with feature meta-data for rows of expression. \cr
#'     \code{samples}    \cr \tab A data.frame with observation meta-data for columns of expression. \cr
#'   } 
#' @param config List with configuration values. Uses the following keys:
#'   \tabular{ll}{
#'     \code{normalization_quantile}   \cr \tab Quantile to use where \code{normalization_method \%in\% c("quantile", "upperquartile")}. \cr
#'   }
#' @param normalization_quantile Numeric in closed interval \code{[0, 1]} specifying quantile to use.
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
#' exprs <- h0testr::sim1(n_obs=6, n_feats=12)$mat
#' feats <- data.frame(feature_id=rownames(exprs))
#' samps <- data.frame(observation_id=colnames(exprs))
#' state <- list(expression=exprs, features=feats, samples=samps)
#' config <- list(normalization_quantile=0.75)
#' state2 <- h0testr::normalize_quantile(state, config)
#' print(state$expression)
#' print(state2$expression)
#' apply(state$expression, 2, quantile, probs=c(0.5, 0.75, 0.9), na.rm=TRUE)
#' apply(state2$expression, 2, quantile, probs=c(0.5, 0.75, 0.9), na.rm=TRUE)

normalize_quantile <- function(state, config, normalization_quantile=NULL, multiplier=1e3) {
  
  check_config(config)
  
  if(!is.matrix(state$expression)) {
    f.err("normalize_quantile: !is.matrix(state$expression)", "\n",
      "class(state$expression):", class(state$expression), config=config)
  }
  if(is.null(normalization_quantile)) normalization_quantile <- config$normalization_quantile
  if(is.null(normalization_quantile)) normalization_quantile <- 0.75
  if(!is.numeric(normalization_quantile)) {
    f.err(
      "normalize_quantile: !is.numeric(normalization_quantile);", 
      "normalization_quantile:", normalization_quantile, ";", 
      "typeof(normalization_quantile):", typeof(normalization_quantile), 
      config=config
    )
  }
  if(normalization_quantile < 0 || normalization_quantile > 1) {
    f.err(
      "normalize_quantile:", 
      "normalization_quantile < 0 || normalization_quantile > 1\n",
      "normalization_quantile:", normalization_quantile, config=config
    )
  }
  
  f <- function(v) {
    multiplier * v / stats::quantile(v, probs=normalization_quantile, na.rm=T)
  }
  state$expression <- apply(state$expression, 2, f)
  
  return(state)
}

#' Normalize expression using total expression in each sample
#' @description
#'   Normalize expression using variant of CPM normalization which excludes missing values.
#' @details 
#'   Inter-observation normalization, based on dividing each expression value 
#'     by total expression in each observation. Makes total expression 
#'     (excluding missing values) equal in each observation.
#' @param state List with elements formatted like the list returned by \code{read_data()}:
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
#' exprs <- h0testr::sim1(n_obs=6, n_feats=12)$mat
#' feats <- data.frame(feature_id=rownames(exprs))
#' samps <- data.frame(observation_id=colnames(exprs))
#' state <- list(expression=exprs, features=feats, samples=samps)
#' config <- list()
#' state2 <- h0testr::normalize_cpm(state, config)
#' print(state$expression)
#' print(state2$expression)
#' apply(state$expression, 2, function(v) sum(v, na.rm=TRUE))
#' apply(state2$expression, 2, function(v) sum(v, na.rm=TRUE))

normalize_cpm <- function(state, config, multiplier=1e6) {

  if(!is.matrix(state$expression)) {
    f.err("normalize_cpm: !is.matrix(state$expression)", "\n",
      "class(state$expression):", class(state$expression), config=config)
  }
  
  f <- function(v) multiplier * (v / sum(v, na.rm=T))
  state$expression <- apply(state$expression, 2, f)
  
  return(state)
}

#' Variance stabilizing transformation
#' @description
#'   Normalize expression using variance stabilizing transformation.
#' @details 
#'   Inter-observation normalization. Uses \code{limma::normalizeVSN()}, which 
#'     wraps \code{vsn::vsn2}. 
#'     Unlike most other normalization methods, results are returned on a 
#'     log2-like scale.
#' @param state List formatted like the list returned by \code{read_data()}:
#'   \tabular{ll}{
#'     \code{expression} \cr \tab Numeric matrix with non-negative expression values. \cr
#'     \code{features}   \cr \tab A data.frame with feature meta-data for rows of expression. \cr
#'     \code{samples}    \cr \tab A data.frame with observation meta-data for columns of expression. \cr
#'   } 
#' @param config List with configuration values. Does not use any keys, so can pass empty list.
#' @param n_pts Scalar minimum number of data points per stratum. See \code{vsn::vsn2}. Default: 42L.
#' @return A list (the processed state) with the following elements:
#'   \tabular{ll}{
#'     \code{expression} \cr \tab Numeric matrix with non-negative expression values. \cr
#'     \code{features}   \cr \tab A data.frame with feature meta-data for rows of expression. \cr
#'     \code{samples}    \cr \tab A data.frame with observation meta-data for columns of expression. \cr
#'   } 
#' @examples
#' set.seed(101)
#' exprs <- h0testr::sim1(n_obs=6, n_feats=50, mnar_c0=-Inf)$mat
#' feats <- data.frame(feature_id=rownames(exprs))
#' samps <- data.frame(observation_id=colnames(exprs))
#' state <- list(expression=exprs, features=feats, samples=samps)
#' config <- list()
#' state2 <- h0testr::normalize_vsn(state, config)
#' head(state$expression)
#' head(state2$expression)
#' summary(apply(state$expression, 2, sd, na.rm=TRUE))
#' summary(apply(log2(state$expression+1), 2, sd, na.rm=TRUE))
#' summary(apply(state2$expression, 2, sd, na.rm=TRUE))
#' summary(apply(state$expression, 1, sd, na.rm=TRUE))
#' summary(apply(log2(state$expression+1), 1, sd, na.rm=TRUE))
#' summary(apply(state2$expression, 1, sd, na.rm=TRUE))

normalize_vsn <- function(state, config, n_pts=42L) {
  
  if(!is.matrix(state$expression)) {
    f.err("normalize_vsn: !is.matrix(state$expression)", "\n",
      "class(state$expression):", class(state$expression), config=config)
  }
  
  n_pts_max <- round(sqrt(nrow(state$expression)))
  if(n_pts > n_pts_max) {
    f.msg("normalize_vsn: n_pts > n_pts_max; n_pts:", n_pts, "\n",
      "Setting n_pts to n_pts_max:", n_pts_max, config=config)
    n_pts <- n_pts_max
  }
  
  state$expression <- limma::normalizeVSN(state$expression, 
    minDataPointsPerStratum=n_pts)

  return(state)
}

#' Cyclic loess normalization
#' @description
#'   Normalize expression using the cyclic-loess algorithm.
#' @details 
#'   Inter-observation normalization using cyclic-loess results
#'     in similar signal distributions across all samples, similar to 
#'     \code{normalize_qquantile()}. This is a slow method, especially if
#'     \code{method \%in\% c("affy", "pairs")}, which scale quadratically. Calls
#'     \code{limma::normalizeCyclicLoess()} under the hood.
#' @param state List with elements formatted like the list returned by `read_data()`:
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
#' exprs <- h0testr::sim1(n_obs=6, n_feats=30, mnar_c0=-Inf)$mat
#' feats <- data.frame(feature_id=rownames(exprs))
#' samps <- data.frame(observation_id=colnames(exprs))
#' state <- list(expression=exprs, features=feats, samples=samps)
#' config <- list()
#' state2 <- h0testr::normalize_loess(state, config)
#' head(state$expression)
#' head(state2$expression)
#' summary(apply(state$expression, 2, sd, na.rm=TRUE))
#' summary(apply(state2$expression, 2, sd, na.rm=TRUE))
#' summary(apply(state$expression, 1, sd, na.rm=TRUE))
#' summary(apply(state2$expression, 1, sd, na.rm=TRUE))

normalize_loess <- function(state, config, span=NULL, method="affy") {

  if(!is.matrix(state$expression)) {
    f.err("normalize_loess: !is.matrix(state$expression)", "\n",
      "class(state$expression):", class(state$expression), config=config)
  }
  
  if(is.null(span)) span <- config$normalization_span
  if(is.null(span)) span <- 0.7
  
  if(!is.numeric(span)) {
    f.err("normalize_loess: !is.numeric(span);", "span:", span, ";", 
      "typeof(span):", typeof(span), config=config)
  }
  if(span < 0 || span > 1) {
    f.err("normalize_loess: span < 0 || span > 1; span:", span, config=config)
  }
  
  state$expression <- limma::normalizeCyclicLoess(state$expression, 
    span=span, method=method)
    
  return(state)
}

#' Old-school 'quantile' normalization
#' @description
#'   Normalize expression values using the 'quantile normalization' algorithm.
#' @details 
#'   Inter-observation normalization resulting in nearly identical
#'     signal distributions across all samples, so all quantiles in \code{0:1} match
#'     across all samples. Calls \code{limma::normalizeQuantiles()} under the hood.
#' @param state List with elements formatted like the list returned by \code{read_data()}:
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
#' exprs <- h0testr::sim1(n_obs=6, n_feats=12, mnar_c0=-Inf)$mat
#' feats <- data.frame(feature_id=rownames(exprs))
#' samps <- data.frame(observation_id=colnames(exprs))
#' state <- list(expression=exprs, features=feats, samples=samps)
#' config <- list()
#' state2 <- h0testr::normalize_qquantile(state, config)
#' print(state$expression)
#' print(state2$expression)
#' 
#' summary(apply(state$expression, 2, sd, na.rm=TRUE))
#' summary(apply(state2$expression, 2, sd, na.rm=TRUE))
#'
#' ## afterwards, all quantiles line up:
#' apply(state$expression, 2, quantile, probs=c(0.1, 0.25, 0.5, 0.75, 0.9), na.rm=TRUE)
#' apply(state2$expression, 2, quantile, probs=c(0.1, 0.25, 0.5, 0.75, 0.9), na.rm=TRUE)

normalize_qquantile <- function(state, config) {

  if(!is.matrix(state$expression)) {
    f.err("normalize_qquantile: !is.matrix(state$expression)", "\n",
      "class(state$expression):", class(state$expression), config=config)
  }
  
  state$expression <- limma::normalizeQuantiles(state$expression)
  
  return(state)
}

#' Normalize expression data using \code{MsCoreUtils} package
#' @description
#'   Normalization using the \code{MsCoreUtils::normalize_matrix()} function.
#' @details
#'   Uses the \code{MsCoreUtils::normalize_matrix()} function to normalize
#'     expression data in \code{state$expression}. Does not affect 
#'     \code{state$features} or \code{state$samples}. Potential values 
#'     for \code{method} are restricted to those that always yield non-negative 
#'     values from non-negative inputs. Also excluded 
#'     \code{c("vsn", "quantiles")} because already available elsewhere 
#'     (as \code{h0testr::normalize_vsn()} and 
#'     \code{h0testr::normalize_qquantile()}). 
#'     Acceptable \code{method} values include:
#'       \code{c("sum", "max", "div.mean", "div.median", "quantiles.robust")}.
#'   See documentation for \code{MsCoreUtils::normalize_matrix()} for details of each method.
#'   See documentation for \code{h0testr::new_config()} 
#'     for more detailed description of configuration parameters. 
#' @param state List with elements like those returned by \code{read_data()}:
#'   \tabular{ll}{
#'     \code{expression} \cr \tab Numeric matrix with non-negative expression values. \cr
#'     \code{features}   \cr \tab A data.frame with feature meta-data for rows of expression. \cr
#'     \code{samples}    \cr \tab A data.frame with observation meta-data for columns of expression. \cr
#'   } 
#' @param config List with configuration values. The only optional setting used is
#'   \code{normalization_method}, which is used to set \code{method} if \code{is.null(method)}.
#' @param method Name (character scalar) of method to use for normalization. See details.
#' @return A list with elements: 
#'   \tabular{ll}{
#'     \code{expression} \cr \tab Numeric matrix with normalized expression values. \cr
#'     \code{features}   \cr \tab A data.frame with feature meta-data for rows of expression. \cr
#'     \code{samples}    \cr \tab A data.frame with observation meta-data for columns of expression. \cr
#'   }
#' @examples
#' set.seed(101)
#' exprs <- h0testr::sim1(n_obs=6, n_feats=12, mnar_c0=-Inf)$mat
#' feats <- data.frame(feature_id=rownames(exprs))
#' samps <- data.frame(observation_id=colnames(exprs))
#' state <- list(expression=exprs, features=feats, samples=samps)
#' config <- list()
#' state2 <- h0testr::normalize_mscoreutils(state, config, method="quantiles.robust")
#' print(state$expression)
#' print(state2$expression)
#'
#' ## afterwards, quantiles more similar across observations:
#' apply(state$expression, 2, quantile, probs=c(0.1, 0.25, 0.5, 0.75, 0.9), na.rm=TRUE)
#' apply(state2$expression, 2, quantile, probs=c(0.1, 0.25, 0.5, 0.75, 0.9), na.rm=TRUE)

normalize_mscoreutils <- function(state, config, method=NULL) {
  
  if(is.null(method) || method %in% "") {
    if(is.null(config$normalization_method) || config$normalization_method %in% "") {
      f.err("normalize_mscoreutils: method and config$normalization_method both unset", 
        config=config)
    }
    method <- config$normalization_method
  }
  
  ## restricted to subset which yield non-negative if fed non-negative:
  allowed <- c("sum", "max", "div.mean", "div.median", "quantiles.robust")
  
  if(!(method %in% allowed)) {
    f.err("normalize_mscoreutils: unrecognized method:", method, "\n", 
      "allowed values:", allowed, config=config)
  }
  
  if(!is.matrix(state$expression)) {
    f.err("normalize_mscoreutils: !is.matrix(state$expression)", "\n",
      "class(state$expression):", class(state$expression), config=config)
  }
  
  ## does not work with integer:
  rnames <- rownames(state$expression)
  state$expression <- apply(state$expression, 2, as.numeric)
  state$expression <- MsCoreUtils::normalize_matrix(state$expression, method=method)
  rownames(state$expression) <- rnames
  
  return(state)
}

#' Get vector of normalization method names
#' @description
#'   Get a vector with acceptable values of \code{method} parameter for \code{h0testr::normalize}.
#' @return
#'   Character vector with names of acceptable values for \code{h0testr::normalize(..., method=)}.
#' @examples
#' normalization_methods <- h0testr::normalize_methods()
#' cat("Available normalization methods:\n")
#' for(method in normalization_methods) {
#'   cat("method:", method, "\n")
#' }

normalize_methods <- function() {
  return(
    c("quantile", "cpm", "vsn", "loess", "qquantile",
      "RLE", "upperquartile", "TMM", "TMMwsp", 
      "sum", "max", "div.mean", "div.median", "quantiles.robust", 
      "log2", "none")
  )
}

#' Inter-sample normalization
#' @description
#'   Normalize expression data to reduce effects of technical differences 
#'     between samples.
#' @details 
#'   Inter-observation normalization using any of the methods available in 
#'     the \code{h0testr} package. See individual methods for more details. 
#'     Normalizes \code{state$expression}. Does not affect 
#'     \code{state$features} or \code{state$samples}.
#'   See documentation for \code{h0testr::new_config()} 
#'     for more detailed description of configuration parameters. 
#' @param state List with elements formatted like the list returned by \code{read_data()}:
#'   \tabular{ll}{
#'     \code{expression} \cr \tab Numeric matrix with non-negative expression values. \cr
#'     \code{features}   \cr \tab A data.frame with feature meta-data for rows of expression. \cr
#'     \code{samples}    \cr \tab A data.frame with observation meta-data for columns of expression. \cr
#'   } 
#' @param config List with configuration values. Uses the following keys:
#'   \tabular{ll}{
#'     \code{normalization_method}   \cr \tab Character scalar in \code{c("vsn","cpm","quantile","qquantile","TMM","TMMwsp","RLE","upperquartile")}. \cr
#'     \code{normalization_quantile} \cr \tab Quantile (numeric between 0 and 1) for \code{normalization_method \%in\% c("quantile", "upperquartile")}. \cr
#'     \code{feat_col}      \cr \tab Column of \code{state$features} matching \code{rownames(state$expression)}. \cr
#'     \code{obs_col}       \cr \tab Column in \code{state$samples} matching \code{colnames(state$expression)}. \cr
#'   }
#' @param method Name of method to use, where scalar \code{method \%in\% h0testr::normalization_methods()}.
#' @param normalization_quantile Quantile for methods \code{c("quantile", "upperquartile")}, where \code{0 <= normalization_quantile <= 1}.
#' @param span Span for method \code{"loess"}, where \code{0 < span < 1}.
#' @return A list (the processed state) with the following elements:
#'   \tabular{ll}{
#'     \code{expression} \cr \tab Numeric matrix with normalized expression values. \cr
#'     \code{features}   \cr \tab Feature meta-data \code{data.frame} corresponding to rows of \code{expression}. \cr
#'     \code{samples}    \cr \tab Observation meta-data \code{data.frame} corresponding to columns of \code{expression}. \cr
#'   } 
#' @examples
#' ## some toy data:
#' set.seed(101)
#' exprs <- h0testr::sim1(n_obs=6, n_feats=12)$mat
#' exprs[, 4:6] <- exprs[, 4:6] * 2
#' feats <- data.frame(feature_id=rownames(exprs))
#' samps <- data.frame(observation_id=colnames(exprs))
#' state <- list(expression=exprs, features=feats, samples=samps)
#'
#' normalization_methods <- h0testr::normalize_methods()
#' cat("Available normalization methods:\n")
#' for(method in normalization_methods) {
#'   cat("method:", method, "\n")
#' }
#' 
#' ## example configured using parameters:
#' config <- list(feat_col="feature_id", obs_col="observation_id")
#' out <- h0testr::normalize(state, config, method="RLE")
#' summary(apply(state$expression, 2, sd, na.rm=TRUE))
#' summary(apply(out$state$expression, 2, sd, na.rm=TRUE))
#' summary(apply(state$expression, 1, sd, na.rm=TRUE))
#' summary(apply(out$state$expression, 1, sd, na.rm=TRUE))
#' 
#' ## example configured using config:
#' config$normalization_method <- "quantile"
#' config$normalization_quantile <- 0.5
#' out <- h0testr::normalize(state, config)
#' summary(apply(state$expression, 2, sd, na.rm=TRUE))
#' summary(apply(out$state$expression, 2, sd, na.rm=TRUE))
#' summary(apply(state$expression, 1, sd, na.rm=TRUE))
#' summary(apply(out$state$expression, 1, sd, na.rm=TRUE))

normalize <- function(state, config, method=NULL, 
    normalization_quantile=NULL, span=NULL) {

  check_config(config)
  f.check_state(state, config)
  
  if(is.null(method) || method %in% "") method <- config$normalization_method
  if(is.null(method) || method %in% "") {
    f.err("normalize: method and config$normalization_method both unset", 
      config=config)
  }
  if(is.null(normalization_quantile)) normalization_quantile <- config$normalization_quantile
  if(is.null(span)) span <- config$normalization_span
  
  f.msg("normalize: normalization_method:", method, 
    "; normalization_quantile:", normalization_quantile, 
    "; normalization_span:", span, config=config)
  
  if(method %in% c("TMM", "TMMwsp", "RLE", "upperquartile")) {
    state <- normalize_edger(state, config, method=method, normalization_quantile=normalization_quantile)
  } else if(method %in% c("sum", "max", "div.mean", "div.median", "quantiles.robust")) {
    state <- normalize_mscoreutils(state, config, method=method)
  } else if(method %in% "quantile") {
    state <- normalize_quantile(state, config, normalization_quantile=normalization_quantile)
  } else if(method %in% "cpm") {
    state <- normalize_cpm(state, config)
  } else if(method %in% "vsn") {
    state <- normalize_vsn(state, config)
  } else if(method %in% "loess") {
    state <- normalize_loess(state, config, span=span)
  } else if(method %in% "qquantile") {
    state <- normalize_qquantile(state, config)
  } else if(method %in% "log2") {
    f.msg("normalization: log2 transformation", config=config)
  } else if(method %in% "none") {
    f.msg("skipping normalization: config$normalization_method %in% 'none'", 
      config=config)
  } else {
    allowed <- normalize_methods()
    f.err("normalize: unexpected config$normalization_method:", 
      config$normalization_method, "\n",
      "allowed:", allowed, config=config
    )
  }
    
  if(!(method %in% c("vsn", "none"))) {
    f.log("transforming data", config=config)
    state$expression <- log2(state$expression + 1)
  } 

  f.check_state(state, config)
  f.report_state(state, config)
  
  if(!is.null(config$run_order)) {
    i <- config$run_order %in% "normalize"
    if(any(i)) {
      prfx <- paste0(which(i)[1] + 2, ".normalized")
    } else {
      prfx <- "normalized"
    }
  }
  f.save_state(state, config, prefix=prfx)
  
  return(list(state=state, config=config))
}
