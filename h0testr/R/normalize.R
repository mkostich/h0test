#' Normalize expression using edgeR
#' @description
#' `f.normalize_edger` normalizes expression using functionality from edgeR.
#' @details Inter-observation normalization. `edgeR::calcNormFactors()` 
#'   called under the hood. Returned values on a counts-per-million scale.
#' @param state List with elements formatted like the list returned by `f.read_data()`:
#'   \tabular{ll}{
#'     \code{expression} \cr \tab Numeric matrix with non-negative expression values. \cr
#'     \code{features}   \cr \tab A data.frame with feature meta-data for rows of expression. \cr
#'     \code{samples}    \cr \tab A data.frame with observation meta-data for columns of expression. \cr
#'   } 
#' @param config List with configuration values.
#' @param method Character in 
#'   `c("TMM", "TMMwsp", "RLE", "upperquartile", "none")`.
#' @param p Numeric between `0` and `1` specifying quantile to use for method 
#'   `upperquartile`.
#' @return A list (the processed state) with the following elements:
#'   \tabular{ll}{
#'     \code{expression} \cr \tab Numeric matrix with non-negative expression values. \cr
#'     \code{features}   \cr \tab A data.frame with feature meta-data for rows of expression. \cr
#'     \code{samples}    \cr \tab A data.frame with observation meta-data for columns of expression. \cr
#'   } 
#' @examples
#'   \dontrun{
#'     state <- list(expression=exprs)
#'     config <- list(norm_method="TMM", log_file="")
#'     state <- f.normalize_edger(state, config)
#'     norm_exprs <- state$expression
#'
#'     config <- list(norm_method="upperquartile", norm_quantile=0.75, log_file="")
#'     state <- f.normalize_edger(state, config)
#'     norm_exprs <- state$expression
#' 
#'     config <- list(log_file="")
#'     state <- f.normalize_edger(state, config, method="upperquartile", p=0.75)
#'     norm_exprs <- state$expression
#'   }

f.normalize_edger <- function(state, config, method=NULL, p=NULL) {

  if(!is.matrix(state$expression)) {
    f.err("f.normalize_edger: !is.matrix(state$expression)", config=config)
  }
  if(is.null(method)) method <- config$norm_method
  if(is.null(p)) p <- config$norm_quantile

  f.log("convert NA to zero", config=config)
  state$expression[is.na(state$expression)] <- 0

  f.log("making edgeR object", config=config)
  obj <- edgeR::DGEList(state$expression)

  f.log("edgeR::calcNormFactors", config=config)
  obj <- edgeR::calcNormFactors(obj, method=method, p=p)

  f.log("making normalized expression", config=config)
  state$expression <- edgeR::cpm(obj, normalized.lib.sizes=T, log=F, prior.count=1)

  f.log("convert zero back to NA", config=config)
  state$expression[state$expression == 0] <- NA

  return(state)
}

#' Normalize expression using selected quantile as scaling factor.
#' @description
#' `f.normalize_quantile` normalizes expression using variant of quantile 
#'   normalization which excludes missing values.
#' @details Inter-observation normalization. Setting p to 0.5 is equivalent 
#'   to median scaling with median calculated after exclusion of missing 
#'   values. Similarly, setting p to 0.75 is upperquartile normalization 
#'   ignoring missing values. 
#' @param state List with elements formatted like the list returned by `f.read_data()`:
#'   \tabular{ll}{
#'     \code{expression} \cr \tab Numeric matrix with non-negative expression values. \cr
#'     \code{features}   \cr \tab A data.frame with feature meta-data for rows of expression. \cr
#'     \code{samples}    \cr \tab A data.frame with observation meta-data for columns of expression. \cr
#'   } 
#' @param config List with configuration values.
#' @param norm_quantile Numeric in `0:1` specifying quantile to use.
#' @param multiplier Numeric used to scale returned values after 
#'   dividing by selected quantile. 
#' @return A list (the processed state) with the following elements:
#'   \tabular{ll}{
#'     \code{expression} \cr \tab Numeric matrix with non-negative expression values. \cr
#'     \code{features}   \cr \tab A data.frame with feature meta-data for rows of expression. \cr
#'     \code{samples}    \cr \tab A data.frame with observation meta-data for columns of expression. \cr
#'   } 
#' @examples
#'   \dontrun{
#'     state <- list(expression=exprs)
#'     config <- list(log_file="", norm_quantile=0.75)
#'     state <- f.normalize_quantile(state, config)
#'     norm_exprs <- state$expression
#'
#'     config <- list(log_file="")
#'     state <- f.normalize_quantile(state, config, norm_quantile=0.75)
#'     norm_exprs <- state$expression
#'  }

f.normalize_quantile <- function(state, config, norm_quantile=NULL, multiplier=1e3) {

  if(!is.matrix(state$expression)) {
    f.err("f.normalize_edger: !is.matrix(state$expression)", config=config)
  }
  if(is.null(norm_quantile)) norm_quantile <- config$norm_quantile
  
  f <- function(v) {
    multiplier * v / stats::quantile(v, probs=norm_quantile, na.rm=T)
  }
  state$expression <- apply(state$expression, 2, f)
  
  return(state)
}

#' Normalize expression using total expression in each sample.
#' @description
#' `f.normalize_cpm` normalizes expression using variant of CPM 
#'   normalization which excludes missing values.
#' @details Inter-observation normalization, based on total expression
#'   in each observation. Makes total expression (excluding missing values)
#'   equal in each observation.
#' @param state List with elements formatted like the list returned by `f.read_data()`:
#'   \tabular{ll}{
#'     \code{expression} \cr \tab Numeric matrix with non-negative expression values. \cr
#'     \code{features}   \cr \tab A data.frame with feature meta-data for rows of expression. \cr
#'     \code{samples}    \cr \tab A data.frame with observation meta-data for columns of expression. \cr
#'   } 
#' @param config List with configuration values.
#' @param multiplier Numeric greater than zero used to scale returned values 
#'   after dividing by total counts in observation. For example, the default 
#'   `multiplier=1e6` yields normalized expression as CPM (counts per million). 
#' @return A list (the processed state) with the following elements:
#'   \tabular{ll}{
#'     \code{expression} \cr \tab Numeric matrix with non-negative expression values. \cr
#'     \code{features}   \cr \tab A data.frame with feature meta-data for rows of expression. \cr
#'     \code{samples}    \cr \tab A data.frame with observation meta-data for columns of expression. \cr
#'   } 
#' @examples
#'   \dontrun{
#'     state <- list(expression=exprs)
#'     config <- list(log_file="")
#'     state <- f.normalize_cpm(state, config)
#'     norm_exprs <- state$expression
#'   }

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
#'   `f.normalize_vsn` normalizes expression using variance stabilizing 
#'     transformation.
#' @details Inter-observation normalization. Under the hood, it calls 
#'   `limma::normalizeVSN()`. Unlike most other normalization methods,
#'   results are returned on a log2-like scale.
#' @param state List with elements formatted like the list returned by `f.read_data()`:
#'   \tabular{ll}{
#'     \code{expression} \cr \tab Numeric matrix with non-negative expression values. \cr
#'     \code{features}   \cr \tab A data.frame with feature meta-data for rows of expression. \cr
#'     \code{samples}    \cr \tab A data.frame with observation meta-data for columns of expression. \cr
#'   } 
#' @param config List with configuration values.
#' @return A list (the processed state) with the following elements:
#'   \tabular{ll}{
#'     \code{expression} \cr \tab Numeric matrix with non-negative expression values. \cr
#'     \code{features}   \cr \tab A data.frame with feature meta-data for rows of expression. \cr
#'     \code{samples}    \cr \tab A data.frame with observation meta-data for columns of expression. \cr
#'   } 
#' @examples
#'   \dontrun{
#'     state <- list(expression=exprs)
#'     config <- list(log_file="")
#'     state <- f.normalize_vsn(state, config)
#'     norm_exprs <- state$expression
#'   }

f.normalize_vsn <- function(state, config) {

  if(!is.matrix(state$expression)) {
    f.err("f.normalize_vsn: !is.matrix(state$expression)", config=config)
  }
  
  state$expression <- limma::normalizeVSN(state$expression)
  return(state)
}

#' Cyclic loess normalization.
#' @description
#' `f.normalize_loess` normalizes expression using the cyclic-loess
#'   algorithm.
#' @details Inter-observation normalization using cyclic-loess results
#'   in similar signal distributions across all samples, similar to 
#'   `f.normalize_qquantile()`. This is a slow method, especially if
#'   `method %in% c("affy", "pairs")`, which scale quadratically. Calls
#'   `limma::normalizeCyclicLoess()` under the hood.
#' @param state List with elements formatted like the list returned by `f.read_data()`:
#'   \tabular{ll}{
#'     \code{expression} \cr \tab Numeric matrix with non-negative expression values. \cr
#'     \code{features}   \cr \tab A data.frame with feature meta-data for rows of expression. \cr
#'     \code{samples}    \cr \tab A data.frame with observation meta-data for columns of expression. \cr
#'   } 
#' @param config List with configuration values.
#' @param span Numeric between 0 and 1 specifying span for loess fit.
#'   Higher numbers result in smoother (less localized) fit.
#' @param method Character in `c("fast", "affy", "pairs")`.
#' @return A list (the processed state) with the following elements:
#'   \tabular{ll}{
#'     \code{expression} \cr \tab Numeric matrix with non-negative expression values. \cr
#'     \code{features}   \cr \tab A data.frame with feature meta-data for rows of expression. \cr
#'     \code{samples}    \cr \tab A data.frame with observation meta-data for columns of expression. \cr
#'   } 
#' @examples
#'   \dontrun{
#'     state <- list(expression=exprs)
#'     config <- list(log_file="")
#'     state <- f.normalize_loess(state, config)
#'     norm_exprs <- state$expression
#'   }

f.normalize_loess <- function(state, config, span=0.7, method="affy") {

  if(!is.matrix(state$expression)) {
    f.err("f.normalize_loess: !is.matrix(state$expression)")
  }
  
  state$expression <- limma::normalizeCyclicLoess(state$expression, span=span, method=method)
  return(state)
}

#' Old-school 'quantile' normalization.
#' @description
#' `f.normalize_qquantile` normalizes data using the 'quantile normalization'
#'   algorithm.
#' @details Inter-observation normalization resulting in nearly identical
#'   signal distributions across all samples, so all quantiles in 0:1 match
#'   across all samples. Calls `limma::normalizeQuantiles()` under the hood.
#' @param state List with elements formatted like the list returned by `f.read_data()`:
#'   \tabular{ll}{
#'     \code{expression} \cr \tab Numeric matrix with non-negative expression values. \cr
#'     \code{features}   \cr \tab A data.frame with feature meta-data for rows of expression. \cr
#'     \code{samples}    \cr \tab A data.frame with observation meta-data for columns of expression. \cr
#'   } 
#' @param config List with configuration values.
#' @return A list (the processed state) with the following elements:
#'   \tabular{ll}{
#'     \code{expression} \cr \tab Numeric matrix with non-negative expression values. \cr
#'     \code{features}   \cr \tab A data.frame with feature meta-data for rows of expression. \cr
#'     \code{samples}    \cr \tab A data.frame with observation meta-data for columns of expression. \cr
#'   } 
#' @examples
#'   \dontrun{
#'     state <- list(expression=exprs)
#'     config <- list(log_file="")
#'     state <- f.normalize_qquantile(state, config)
#'     norm_exprs <- state$expression
#'   }

f.normalize_qquantile <- function(state, config) {

  if(!is.matrix(state$expression)) {
    f.err("f.normalize_qquantile: !is.matrix(state$expression)", config=config)
  }
  
  state$expression <- limma::normalizeQuantiles(state$expression)
  return(state)
}

