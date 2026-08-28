## observation with no measured value has nothing for 
##   inter-observation method to work with: 

f.obs_measured <- function(state, config, method, fn, refuse=FALSE) {

  i <- apply(state$expression, 2, function(v) all(is.na(v)))

  if(!any(i)) return(invisible(NULL))

  nom <- colnames(state$expression)[i]
  if(is.null(nom)) nom <- as.character(which(i))

  msg <- list(sum(i), "of", length(i), "observations have no measured value,",
    "so there is nothing in them for normalization method", method,
    "to work with;", "\n",
    "  observations:", paste(utils::head(nom, 5), collapse=", "), "\n",
    "  drop them before normalizing: filter_observations() keeps an",
    "observation only if it measures at least config$n_features_min features",
    config=config)

  if(refuse) {
    do.call(f.err, c(paste0(fn, ": cannot normalize:"), msg))
  } else {
    do.call(f.msg, c(paste0(fn, ": WARNING:"), msg))
  }

  return(invisible(NULL))
}

## the methods that come back with the same value for every observation when
##   state$expression has a single feature:

f.one_feature_methods <- function() {
  return(
    c("quantile", "cpm", "loess", "RLE", "upperquartile", "TMM", "TMMwsp",
      "div.mean", "div.median", "quantiles.robust")
  )
}

## reported rather than refused: a one-feature matrix:

f.one_feature <- function(state, config, method, fn) {

  if(nrow(state$expression) > 1) return(invisible(NULL))
  if(!(method %in% f.one_feature_methods())) return(invisible(NULL))

  f.msg(paste0(fn, ": WARNING:"), "state$expression has",
    nrow(state$expression), "feature, and normalization method", method,
    "divides each observation by a statistic of that observation, which with",
    "one feature is the value itself;", "\n",
    "  every observation will come back with the same value, so every",
    "difference between observations is erased;", "\n",
    "  of the other methods, sum and max come through, MsCoreUtils",
    "normalizing those per feature rather than per observation, log2 and none",
    "do not normalize, and vsn and qquantile refuse a single feature",
    config=config)

  return(invisible(NULL))
}

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
#' @param method Character scalar in set
#'   \code{c("RLE", "upperquartile", "TMM", "TMMwsp", "none")}. Can be set with
#'   \code{config$normalization_method}.
#' @param normalization_quantile Single number where \code{0 <= normalization_quantile <= 1.0}; quantile for method 
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
#' ## normalization_quantile not needed for e.g. method="RLE":
#' state2 <- h0testr::normalize_edger(state, config, method="RLE")
#' print(state$expression)
#' print(state2$expression)

normalize_edger <- function(state, config, method=NULL, normalization_quantile=NULL) {

  check_config(config)

  if(!is.matrix(state$expression)) {
    f.err("normalize_edger: !is.matrix(state$expression)", "\n",
      "class(state$expression):", class(state$expression), config=config)
  }
  if(length(method) %in% 0) method <- config$normalization_method
  if(length(method) %in% 0) method <- "RLE"
  allowed <- c("RLE", "upperquartile", "TMM", "TMMwsp", "none")

  if(!(length(method) %in% 1 && is.character(method) && !is.na(method) &&
      method %in% allowed)) {
    f.err("normalize_edger: method has to be one of",
      paste(allowed, collapse=" "), ";", "\n", " method:", method,
      "; class:", class(method), "; length:", length(method), config=config)
  }
  
  if(length(normalization_quantile) %in% 0) {
    normalization_quantile <- config$normalization_quantile
  }
  if(length(normalization_quantile) %in% 0) normalization_quantile <- 0.75

  if(!(length(normalization_quantile) %in% 1 &&
      is.numeric(normalization_quantile) && !is.na(normalization_quantile))) {
    f.err("normalize_edger: normalization_quantile has to be a single number;",
      "\n", " normalization_quantile:", normalization_quantile, "; class:",
      class(normalization_quantile), "; length:",
      length(normalization_quantile), config=config)
  }
  if(normalization_quantile < 0 || normalization_quantile > 1) {
    f.err(
      "normalize_edger:", 
      "normalization_quantile < 0 || normalization_quantile > 1\n",
      "normalization_quantile:", normalization_quantile, config=config
    )
  }
  
  f.obs_measured(state, config, method, "normalize_edger", refuse=TRUE)

  ## edgeR's RLE takes reference from features measured in every observation: 

  if(method %in% "RLE") {
    n_complete <- sum(apply(state$expression, 1, function(v) !anyNA(v)))
    if(n_complete %in% 0) {
      f.err("normalize_edger: method RLE cannot be used on this matrix:", "\n",
        "  no feature is measured in every observation:",
        n_complete, "of", nrow(state$expression), "features;", "\n",
        "  edgeR's RLE takes its reference from the features that are, so",
        "without one edgeR::calcNormFactors() returns NA normalization factors",
        "and edgeR::cpm() then stops with 'library sizes should be finite and",
        "non-negative';", "\n",
        "  method 'TMM', 'TMMwsp' and 'upperquartile' do not need a feature",
        "measured in every observation. See h0testr::normalize_methods()",
        config=config)
    }
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
#' @param normalization_quantile Single number in closed interval \code{[0, 1]} specifying quantile to use.
#' @param multiplier Single finite number greater than zero, used to scale returned values after 
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
  if(length(normalization_quantile) %in% 0) {
    normalization_quantile <- config$normalization_quantile
  }
  if(length(normalization_quantile) %in% 0) normalization_quantile <- 0.75

  if(!(length(normalization_quantile) %in% 1 &&
      is.numeric(normalization_quantile) && !is.na(normalization_quantile))) {
    f.err("normalize_quantile: normalization_quantile has to be a single",
      "number;", "\n", " normalization_quantile:", normalization_quantile,
      "; class:", class(normalization_quantile), "; length:",
      length(normalization_quantile), config=config)
  }
  if(normalization_quantile < 0 || normalization_quantile > 1) {
    f.err(
      "normalize_quantile:", 
      "normalization_quantile < 0 || normalization_quantile > 1\n",
      "normalization_quantile:", normalization_quantile, config=config
    )
  }

  if(!(length(multiplier) %in% 1 && is.numeric(multiplier) &&
      is.finite(multiplier) && multiplier > 0)) {
    f.err("normalize_quantile: multiplier has to be a single finite number greater than zero;",
      "\n", " multiplier:", multiplier, "; class:", class(multiplier),
      "; length:", length(multiplier), config=config)
  }
  f <- function(v) {
    multiplier * v / stats::quantile(v, probs=normalization_quantile, na.rm=T)
  }
  
  for(j in seq_len(ncol(state$expression))) {
    state$expression[, j] <- f(state$expression[, j])
  }
  
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
#' @param config List with configuration values. Does not use any keys, so can pass an empty list; checked with \code{h0testr::check_config()}.
#' @param multiplier Single finite number greater than zero, used to scale returned values 
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

  check_config(config)

  if(!is.matrix(state$expression)) {
    f.err("normalize_cpm: !is.matrix(state$expression)", "\n",
      "class(state$expression):", class(state$expression), config=config)
  }

  if(!(length(multiplier) %in% 1 && is.numeric(multiplier) &&
      is.finite(multiplier) && multiplier > 0)) {
    f.err("normalize_cpm: multiplier has to be a single finite number greater than zero;",
      "\n", " multiplier:", multiplier, "; class:", class(multiplier),
      "; length:", length(multiplier), config=config)
  }
  f <- function(v) multiplier * (v / sum(v, na.rm=T))
  
  for(j in seq_len(ncol(state$expression))) {
    state$expression[, j] <- f(state$expression[, j])
  }
  
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
#'   A single feature is refused: \code{vsn::vsn2} fits a transformation per
#'     stratum; one feature is not enough to fit, so is refused.
#' @param state List formatted like the list returned by \code{read_data()}:
#'   \tabular{ll}{
#'     \code{expression} \cr \tab Numeric matrix with non-negative expression values. \cr
#'     \code{features}   \cr \tab A data.frame with feature meta-data for rows of expression. \cr
#'     \code{samples}    \cr \tab A data.frame with observation meta-data for columns of expression. \cr
#'   } 
#' @param config List with configuration values. Does not use any keys, so can pass an empty list; checked with \code{h0testr::check_config()}.
#' @param n_pts Scalar minimum number of data points per stratum, as a whole number of 1 or more.
#'   Reduced to \code{round(sqrt(nrow(state$expression)))} with a message if it exceeds that.
#'   See \code{vsn::vsn2}. Default: 42L.
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

  check_config(config)

  f.need_pkgs("vsn", "normalize_vsn", config)
  
  if(!is.matrix(state$expression)) {
    f.err("normalize_vsn: !is.matrix(state$expression)", "\n",
      "class(state$expression):", class(state$expression), config=config)
  }
  
  ## n_pts is handed to vsn2 as minDataPointsPerStratum: 

  if(!(length(n_pts) %in% 1 && is.numeric(n_pts) && !is.na(n_pts))) {
    f.err("normalize_vsn: n_pts has to be a single number;", "\n",
      " n_pts:", n_pts, "; class:", class(n_pts), "; length:", length(n_pts),
      config=config)
  }
  if(n_pts != round(n_pts) || n_pts < 1) {
    f.err("normalize_vsn: n_pts has to be a whole number of 1 or more;",
      "n_pts:", n_pts, config=config)
  }

  n_pts_max <- round(sqrt(nrow(state$expression)))
  if(n_pts > n_pts_max) {
    f.msg("normalize_vsn: n_pts > n_pts_max; n_pts:", n_pts, "\n",
      "Setting n_pts to n_pts_max:", n_pts_max, config=config)
    n_pts <- n_pts_max
  }
  
  ## single feature refused: 

  if(nrow(state$expression) < 2) {
    f.err("normalize_vsn: needs at least 2 features;", "\n",
      " nrow(state$expression):", nrow(state$expression), "\n",
      "vsn fits a transformation per stratum, and one feature is not enough",
      "to fit one to: where it returns anything at all it returns the same",
      "value for every observation, and otherwise it reports 'L-BFGS-B needs",
      "finite values of fn' or 'Elements of argument Sstrat must be in",
      "ascending order'", config=config)
  }

  state$expression <- limma::normalizeVSN(state$expression, 
    minDataPointsPerStratum=n_pts)

  return(state)
}

#' Cyclic loess normalization
#' @description
#'   Normalize expression using cyclic-loess.
#' @details
#'   Inter-observation normalization using cyclic-loess results
#'     in similar signal distributions across all samples, similar to
#'     \code{normalize_qquantile()}. This is a slow method, especially if
#'     \code{method \%in\% c("affy", "pairs")}, which scale quadratically. Calls
#'     \code{limma::normalizeCyclicLoess()} under the hood.
#'   Only \code{method="fast"} can be used on data with missing values, so it is
#'     the default here, as it is in \code{limma}. The other two normalize each
#'     pair of observations on their difference, which is missing for any feature
#'     the pair does not both measure; since the adjustments of all pairs are
#'     summed, one missing value removes that feature from every observation
#'     rather than from the one it was missing from. Rather than return that
#'     quietly, \code{method \%in\% c("affy", "pairs")} is an error when
#'     \code{state$expression} has any missing value.
#'   Cyclic loess is an \emph{additive} correction, so \code{state$expression} is
#'     expected on a log scale, but raw intensities are not refused.
#' @param state List with elements formatted like the list returned by `read_data()`:
#'   \tabular{ll}{
#'     \code{expression} \cr \tab Numeric matrix with non-negative expression values. \cr
#'     \code{features}   \cr \tab A data.frame with feature meta-data for rows of expression. \cr
#'     \code{samples}    \cr \tab A data.frame with observation meta-data for columns of expression. \cr
#'   } 
#' @param config List with configuration values, checked with \code{h0testr::check_config()}. Uses the following key:
#'   \tabular{ll}{
#'     \code{normalization_span} \cr \tab Span (numeric between 0 and 1) for the loess fit; read when \code{span} is \code{NULL}. \cr
#'   }
#' @param span Single number between 0 and 1 specifying span for loess fit.
#'   Higher numbers result in smoother (less localized) fit. Defaults to
#'   \code{config$normalization_span}, or to \code{0.7} when that is unset as well.
#' @param method Character in \code{c("fast", "affy", "pairs")}. Default:
#'   \code{"fast"}, the only one usable on data with missing values; see details.
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

normalize_loess <- function(state, config, span=NULL, method="fast") {

  check_config(config)

  if(!is.matrix(state$expression)) {
    f.err("normalize_loess: !is.matrix(state$expression)", "\n",
      "class(state$expression):", class(state$expression), config=config)
  }

  allowed <- c("fast", "affy", "pairs")

  if(!(length(method) %in% 1 && is.character(method) && !is.na(method) &&
      method %in% allowed)) {
    f.err("normalize_loess: method has to be one of",
      paste(allowed, collapse=" "), ";", "\n", " method:", method,
      "; class:", class(method), "; length:", length(method), config=config)
  }

  ## limma's "affy" and "pairs" normalize each pair of observations on
  ##   x[, j] - x[, i], which is missing for any feature the pair does not both
  ##   measure:

  n_na <- sum(is.na(state$expression))

  if(method %in% c("affy", "pairs") && n_na > 0) {
    f.err("normalize_loess: method", method, "cannot be used on data with",
      "missing values;", "\n", " missing values in state$expression:", n_na,
      "of", length(state$expression), "\n",
      "it normalizes each pair of observations on their difference, so a",
      "feature the pair does not both measure is dropped from every",
      "observation, not just from the one it was missing from;", "\n",
      "use method='fast' (the default), which fits each observation against the",
      "mean of the values each feature does have", config=config)
  }

  if(length(span) %in% 0) span <- config$normalization_span
  if(length(span) %in% 0) span <- 0.7

  ## as with method above: !is.numeric() let a span of length 2 and an NA
  ##   through to the range check, which is where R reported them:

  if(!(length(span) %in% 1 && is.numeric(span) && !is.na(span))) {
    f.err("normalize_loess: span has to be a single number;", "\n",
      " span:", span, "; class:", class(span), "; length:", length(span),
      config=config)
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
#'   Normalize expression values using 'quantile normalization'.
#' @details 
#'   Inter-observation normalization resulting in nearly identical
#'     signal distributions across all samples, so all quantiles in \code{0:1} match
#'     across all samples. Calls \code{limma::normalizeQuantiles()} under the hood.
#'   An observation with no measured value is refused. Drop such observations with
#'     \code{h0testr::filter_observations()} first.
#'   A single feature is refused as well: quantiles are matched by
#'     interpolating between the values of each observation, and a single
#'     value has nothing to interpolate between, which \code{limma} reports
#'     as \code{"need at least two non-NA values to interpolate"}.
#' @param state List with elements formatted like the list returned by \code{read_data()}:
#'   \tabular{ll}{
#'     \code{expression} \cr \tab Numeric matrix with non-negative expression values. \cr
#'     \code{features}   \cr \tab A data.frame with feature meta-data for rows of expression. \cr
#'     \code{samples}    \cr \tab A data.frame with observation meta-data for columns of expression. \cr
#'   } 
#' @param config List with configuration values. Does not use any keys, so can pass an empty list; checked with \code{h0testr::check_config()}.
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

  check_config(config)

  if(!is.matrix(state$expression)) {
    f.err("normalize_qquantile: !is.matrix(state$expression)", "\n",
      "class(state$expression):", class(state$expression), config=config)
  }

  f.obs_measured(state, config, "qquantile", "normalize_qquantile", refuse=TRUE)

  ## limma::normalizeQuantiles() matches quantiles across observations by
  ##   interpolating between the values of each one; needs 2+ values:
  
  if(nrow(state$expression) < 2) {
    f.err("normalize_qquantile: needs at least 2 features;", "\n",
      " nrow(state$expression):", nrow(state$expression), "\n",
      "quantiles are matched across observations by interpolating between",
      "the values of each one, and a single feature has nothing to",
      "interpolate between; limma reports this as 'need at least two",
      "non-NA values to interpolate'", config=config)
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
#'   \code{"div.mean"}, \code{"div.median"} and \code{"quantiles.robust"}
#'     work across observations, while \code{"sum"} and \code{"max"} work
#'     across features:
#'     \code{MsCoreUtils::normalize_matrix(m, method="sum")} is
#'     \code{m / rowSums(m)} and \code{method="max"} is
#'     \code{m / apply(m, 1, max)}, checked with \code{all.equal()} against
#'     both the per-row and the per-observation forms on a 60 x 6 matrix. So
#'     those two rescale each feature by its own total or maximum rather than
#'     equalizing anything between observations.
#'   \code{"quantiles.robust"} is an error when \code{state$expression} has
#'     any missing value. \code{preprocessCore} assigns values by rank within
#'     each observation, and a missing value takes a rank of its own, so it
#'     comes back at that rank in every observation, on whichever feature sorts
#'     there rather than on the one it was missing from. The other four methods 
#'     here return the missing values they were given, in the cells they were given them in;
#'     \code{h0testr::normalize_qquantile()} matches quantiles across
#'     observations without the robust trimming and without moving anything.
#'   See documentation for \code{MsCoreUtils::normalize_matrix()} for details of each method.
#'   See documentation for \code{h0testr::new_config()} 
#'     for more detailed description of configuration parameters. 
#' @param state List with elements like those returned by \code{read_data()}:
#'   \tabular{ll}{
#'     \code{expression} \cr \tab Numeric matrix with non-negative expression values. \cr
#'     \code{features}   \cr \tab A data.frame with feature meta-data for rows of expression. \cr
#'     \code{samples}    \cr \tab A data.frame with observation meta-data for columns of expression. \cr
#'   } 
#' @param config List with configuration values, checked with \code{h0testr::check_config()}. The only optional setting used is
#'   \code{normalization_method}, which is used to set \code{method} if \code{is.null(method)}.
#' @param method Name (character scalar) of method to use for normalization.
#'   \code{"quantiles.robust"} cannot be used on data with any missing value.
#'   See details.
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

  check_config(config)

  if(length(method) %in% 0 || isTRUE(method %in% "")) {
    if(length(config$normalization_method) %in% 0 ||
        isTRUE(config$normalization_method %in% "")) {
      f.err("normalize_mscoreutils: method and config$normalization_method",
        "both unset", config=config)
    }
    method <- config$normalization_method
  }
  
  ## restricted to subset which yield non-negative if fed non-negative:
  allowed <- c("sum", "max", "div.mean", "div.median", "quantiles.robust")
  
  ## the method can arrive from config as well as from the argument:

  if(!(length(method) %in% 1 && is.character(method) && !is.na(method) &&
      method %in% allowed)) {
    f.err("normalize_mscoreutils: method has to be one of",
      paste(allowed, collapse=" "), ";", "\n", " method:", method,
      "; class:", class(method), "; length:", length(method), config=config)
  }
  
  if(!is.matrix(state$expression)) {
    f.err("normalize_mscoreutils: !is.matrix(state$expression)", "\n",
      "class(state$expression):", class(state$expression), config=config)
  }

  ## missing value takes a rank of its own:

  n_na <- sum(is.na(state$expression))

  if(method %in% "quantiles.robust" && n_na > 0) {
    f.err("normalize_mscoreutils: method", method, "cannot be used on data",
      "with missing values;", "\n",
      " missing values in state$expression:", n_na, "of",
      length(state$expression), "\n",
      "it assigns values by rank within each observation, so a missing value",
      "comes back at that rank in every observation, in whatever feature sorts",
      "there rather than in the one it was missing from;", "\n",
      "the other methods here (sum, max, div.mean, div.median) leave missing",
      "values where they are, as does h0testr::normalize_qquantile(), which",
      "matches quantiles across observations without the robust trimming",
      config=config)
  }

  ## MsCoreUtils::normalize_matrix() does not work with integer input:

  storage.mode(state$expression) <- "double"
  state$expression <- MsCoreUtils::normalize_matrix(state$expression, method=method)
  
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
#'     between obervations.
#' @details 
#'   Inter-observation normalization using any of the methods available in 
#'     the \code{h0testr} package. See individual methods for more details. 
#'     Normalizes \code{state$expression}. Does not affect 
#'     \code{state$features} or \code{state$samples}.
#'   Every method other than \code{"none"} leaves data log transformed, so
#'     \code{config$is_log_transformed} is set to \code{TRUE} on return and the
#'     rest of the workflow reads the scale from there rather than tracking it.
#'     For that reason a method other than \code{"none"} is an error when
#'     \code{config$is_log_transformed} is already \code{TRUE}. \code{"vsn"}
#'     counts as transforming: its output is arsinh-scaled rather than log2, but
#'     is not to be transformed again either.
#'   \code{"loess"} is transformed before normalizing rather than after it, since
#'     cyclic loess is an additive correction: on raw intensities it returns
#'     negative fitted values for the smallest measurements, and \code{log2(x +
#'     1)} of anything at or below \code{-1} is \code{NaN}, so the old order
#'     deleted those measurements. See \code{h0testr::normalize_loess()}.
#'   \code{"quantiles.robust"} is an error when \code{state$expression} has
#'     any missing value, since it assigns values by rank within each
#'     observation and so returns a missing value at that rank in every
#'     observation, on whichever feature sorts there rather than on the one it
#'     was missing from. Normalization comes before imputation, so that is most
#'     matrices; \code{"qquantile"} matches quantiles across observations
#'     without changing anything. See \code{h0testr::normalize_mscoreutils()}.
#'   An observation with no measured value is reported whatever the method, since
#'     there is nothing for an inter-observation method to work with. Such observations
#'     should be dropped beforehand by, \code{h0testr::filter_observations()}.
#'   A single feature is normalized by whichever method can do it, with a
#'     \code{WARNING} from here: ten of the sixteen methods return the same
#'     value for every observation, so every difference between observations is
#'     erased. Each of the ten divides an observation by a statistic of that
#'     observation, and with one feature that statistic is the value itself.
#'     Measured on a 1 x 6 matrix: \code{"quantile"}, \code{"cpm"},
#'     \code{"loess"}, \code{"RLE"}, \code{"upperquartile"}, \code{"TMM"},
#'     \code{"TMMwsp"}, \code{"div.mean"}, \code{"div.median"} and
#'     \code{"quantiles.robust"}. \code{"sum"} and \code{"max"} come through
#'     instead, \code{MsCoreUtils} normalizing those per feature rather than
#'     per observation, which is also the exception to the description of these
#'     methods as inter-observation:
#'     \code{MsCoreUtils::normalize_matrix(m, "sum")} is
#'     \code{m / rowSums(m)}. \code{"qquantile"} needs two features to
#'     interpolate between and \code{"vsn"} needs two to fit a transformation
#'     to; both refuse one, in \code{h0testr::normalize_qquantile()} and
#'     \code{h0testr::normalize_vsn()} respectively.
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
#'     \code{normalization_method}   \cr \tab Character scalar naming one of the methods \code{h0testr::normalize_methods()} returns. \cr
#'     \code{normalization_quantile} \cr \tab Quantile (numeric between 0 and 1) for \code{normalization_method \%in\% c("quantile", "upperquartile")}. \cr
#'     \code{is_log_transformed} \cr \tab Optional logical; if \code{TRUE}, any method other than \code{"none"} is an error, since it would transform the data a second time. \cr
#'     \code{feat_col}      \cr \tab Column of \code{state$features} matching \code{rownames(state$expression)}. \cr
#'     \code{obs_col}       \cr \tab Column in \code{state$samples} matching \code{colnames(state$expression)}. \cr
#'   }
#' @param method Name of method to use, where scalar \code{method \%in\% h0testr::normalize_methods()}.
#'   Checked before anything else happens, so a name no method goes by is
#'   refused without the data having been touched.
#' @param normalization_quantile Quantile for methods \code{c("quantile", "upperquartile")}, where \code{0 <= normalization_quantile <= 1}.
#' @param span Span for method \code{"loess"}, where \code{0 < span < 1}.
#' @return A list with the following two elements:
#'   \tabular{ll}{
#'     \code{state}  \cr \tab The processed state; see below. \cr
#'     \code{config} \cr \tab The configuration, with \code{is_log_transformed} set to \code{TRUE} unless the method was \code{"none"}. \cr
#'   }
#'   The element \code{state} is a list with the following three elements:
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

  if(length(method) %in% 0 || isTRUE(method %in% "")) {
    method <- config$normalization_method
  }
  if(length(method) %in% 0 || isTRUE(method %in% "")) {
    f.err("normalize: method and config$normalization_method both unset",
      config=config)
  }

  ## method named here rather than at end of if/else chain below:

  allowed <- normalize_methods()

  if(!(length(method) %in% 1 && is.character(method) && !is.na(method) &&
      method %in% allowed)) {
    f.err("normalize: method has to be one of", paste(allowed, collapse=" "),
      ";", "\n", " method:", method, "; class:", class(method),
      "; length:", length(method), config=config)
  }

  if(is.null(normalization_quantile)) normalization_quantile <- config$normalization_quantile
  if(is.null(span)) span <- config$normalization_span
  
  f.msg("normalize: normalization_method:", method,
    "; normalization_quantile:", normalization_quantile,
    "; normalization_span:", span, config=config)

  ## every method other than "none" ends by transforming data:

  if(!(method %in% "none") && isTRUE(config$is_log_transformed)) {
    f.err("normalize: config$is_log_transformed is TRUE, so state$expression is",
      "already on a log-like scale, but normalization method is",
      method, ", which would transform it again;", "\n",
      "set method (or config$normalization_method) to 'none', or set",
      "config$is_log_transformed to FALSE if the input really is raw",
      config=config)
  }

  ## observation with nothing measured is named in the log method:

  if(!(method %in% c("TMM", "TMMwsp", "RLE", "upperquartile", "qquantile",
      "none"))) {
    f.obs_measured(state, config, method, "normalize")
  }

  ## a single feature is the other thing there is nothing to work with, and it
  ##   is the quiet one: 

  f.one_feature(state, config, method, "normalize")

  ## cyclic loess is an *additive* correction; limma documents log-expression
  ##   input for it, so the transformation belongs before the fit:

  log_first <- method %in% "loess"

  if(log_first) {
    f.log("transforming data before the loess fit", config=config)
    state$expression <- log2(state$expression + 1)
  }

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
    ## not reachable from outside: 
    f.err("normalize: no branch for method:", method, ";", "\n",
      "it is one of the names h0testr::normalize_methods() returns, but",
      "nothing in h0testr::normalize() handles it", config=config)
  }
    
  if(!log_first && !(method %in% c("vsn", "none"))) {
    f.log("transforming data", config=config)
    state$expression <- log2(state$expression + 1)
  }

  ## record scale so that no caller has to track it: 
  if(!(method %in% "none")) config$is_log_transformed <- TRUE

  ## config$log_from_raw is narrower, and says that an exact 0 in this matrix
  ##   cannot be a measurement, which is what f.check_state() then watches for:

  if(!(method %in% "none")) config$log_from_raw <- FALSE

  if(!(method %in% c("vsn", "none"))) {
    n_zero <- sum(state$expression %in% 0)     ## a count; NAs are not counted
    if(n_zero %in% 0) {
      config$log_from_raw <- TRUE
    } else {
      f.msg("normalize: WARNING:", n_zero, "exact zeros in the transformed",
        "data, which log2(x + 1) of a positive value cannot produce, so",
        method, "must have emitted them;", "\n",
        "  they cannot be told apart from a missing value written as a",
        "measurement, so that check is disabled for this run", config=config)
    }
  }

  f.check_state(state, config)
  f.report_state(state, config)
  
  prfx <- "normalized"
  if(!is.null(config$run_order)) {
    i <- config$run_order %in% "normalize"
    if(any(i)) {
      prfx <- paste0(which(i)[1] + 2, ".normalized")
    }
  }
  f.save_state(state, config, prefix=prfx)
  
  return(list(state=state, config=config))
}
