#' Normalize expression using \code{edgeR}
#' @description
#'   Normalize expression using functionality from \code{edgeR} package.
#' @details Inter-observation normalization. \code{edgeR::calcNormFactors()} 
#'   called under the hood. Returned values on a counts-per-million scale.
#' @param state List with elements formatted like the list returned by \code{f.read_data()}:
#'   \tabular{ll}{
#'     \code{expression} \cr \tab Numeric matrix with non-negative expression values. \cr
#'     \code{features}   \cr \tab A data.frame with feature meta-data for rows of expression. \cr
#'     \code{samples}    \cr \tab A data.frame with observation meta-data for columns of expression. \cr
#'   } 
#' @param config List with configuration values.
#' @param method Character in set
#'   \code{c("TMM", "TMMwsp", "RLE", "upperquartile", "none")}.
#' @param p Numeric in closed interval \code{[0, 1]} specifying quantile to use for method 
#'   \code{upperquartile}.
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
#' @param config List with configuration values.
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
#' @param config List with configuration values.
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
#'   Normalize expression using variance stabilizing transformation.
#' @details 
#'   Inter-observation normalization. Under the hood, it calls 
#'     \code{limma::normalizeVSN()}. Unlike most other normalization methods,
#'     results are returned on a log2-like scale.
#' @param state List with elements formatted like the list returned by \code{f.read_data()}:
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
#' @param config List with configuration values.
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

#' Inter-sample normalization
#' @description
#'   Normalize expression data to reduce effects of technical differences 
#'     between samples.
#' @details 
#'   Inter-observation normalization using any of the methods available in this package.
#' @param state List with elements formatted like the list returned by \code{f.read_data()}:
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
#'     config <- list(norm_method="RLE", log_file="")
#'     out <- f.normalize(state, config)
#'     config2 <- out$config
#'     state2 <- out$state
#'     exprs2 <- state2$expression
#'
#'     state <- list(expression=exprs)
#'     config <- list(norm_method="quantile", norm_quantile=0.75, log_file="")
#'     out <- f.normalize(state, config)
#'     config2 <- out$config
#'     state2 <- out$state
#'     exprs2 <- state2$expression
#'
#'     state <- list(expression=exprs)
#'     config <- list(norm_method="cpm", multiplier=1e3, log_file="")
#'     out <- f.normalize(state, config)
#'     config2 <- out$config
#'     state2 <- out$state
#'     exprs2 <- state2$expression
#'   }

f.normalize <- function(state, config) {
  
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
  } else if(config$norm_method %in% "none") {
    f.msg("skipping normalization: config$norm_method %in% 'none'", config=config)
  } else f.err("unexpected config$norm_method:", config$norm_method, config=config)
  
  if(!(config$norm_method %in% "vsn")) {
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


#' Combine technical replicates
#' @description
#'   Combine expression signals from technical replicates.
#' @details Combines signals for each feature across technical replicates by taking median.
#' @param state List with elements formatted like the list returned by \code{f.read_data()}:
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
#'     config <- list(obs_id_col="observation_id", sample_id_col="sample_id", log_file="")
#'     state <- list(expression=exprs, samples=samps, features=feats)
#'     cat("config$obs_col:", config$obs_col, "\n")
#'     out <- f.combine_reps(state, config)
#'     config2 <- out$config
#'     cat("config2$obs_col:", config2$obs_col, "\n")
#'     state2 <- out$state
#'     exprs2 <- state2$expression
#'     samps2 <- state2$samples
#'     f.do_stuff(state2, config2)
#'   }

f.combine_reps <- function(state, config) {

  f <- function(v, s) tapply(v, s, stats::median, na.rm=T)
  sample_ids <- state$samples[, config$sample_id_col, drop=T]
  state$expression <- t(apply(state$expression, 1, f, sample_ids))
  state$samples <- state$samples[!duplicated(sample_ids), ]
  state$samples[, config$obs_id_col] <- NULL                  ## !!!!!
  sample_ids <- state$samples[, config$sample_id_col, drop=T]
  if(!all(sample_ids %in% colnames(state$expression))) {
    f.err("!all(samples[, config$sample_id_col] %in% colnames(expression))", config=config)
  }  
  state$expression <- state$expression[, sample_ids]
  config$obs_col <- config$sample_id_col

  f.check_state(state, config)
  f.report_state(state, config)
  i <- config$run_order %in% "combine_reps"
  prfx <- paste0(which(i)[1] + 2, ".combined")
  f.save_state(state, config, prefix=prfx)
  
  return(list(state=state, config=config))
}

