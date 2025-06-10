#' Report configuration
#' @description
#' `f.report_config` reports configuration used for run, and checks if
#'   `config$test_term` is compatible with `config$frm`.
#' @details Report written to config$log_file; if `config$log_file == ""`,
#'   written to standard out (console or terminal). Only supports non-lists
#'     and lists of non-lists (not lists of lists) as `config` values.
#'   Throws error if `config$test_term` is not compatible with `config$frm`.
#' @param config List with configuration values
#' @return NULL
#' @examples
#'   \dontrun{
#'     config <- list(a=1, b="accorn", c=~x1+x2+x1:x2, d=list(e=1, b=FALSE))
#'     f.report_config(config)
#'   }

f.report_config <- function(config) {
  for(k1 in names(config)) {
    v1 <- config[[k1]]
    if(is.list(v1)) {
      for(k2 in names(v1)) {
        f.msg(k1, ":", k2, ":", paste(as.character(v1), sep=", "), config=config)
      }
    } else f.msg(k1, ":", paste(as.character(v1), sep=", "), config=config)
  }
  
  ## check config$test_term compatible with config$frm; throws error if not, 
  ##   else returns NULL:

  trms <- as.character(config$frm)[2]
  trms <- gsub("[[:space:]]+", "", trms)
  trms <- unlist(strsplit(trms, split="\\+"))
  
  if(!(config$test_term %in% trms)) {
    f.err("config$test_term", config$test_term, "not in config$frm", 
      as.character(config$frm), config=config)
  }
}

###############################################################################

#' Read data
#' @description
#' `f.read_data` reads expression data, feature metadata, and observation 
#'   metadata.
#' @details 
#'   Expression data from `config$data_file_in` in `config$dir_in`.
#'   Feature metadata from `config$feature_file_in` in `config$dir_in`.
#'   Observation metadata from `config$sample_file_in` in `config$dir_in`.
#' @param config List with configuration values.
#' @return A list (the initial state) with the following elements:
#'   \tabular{ll}{
#'     \code{expression} \cr \tab Numeric matrix with non-negative expression values. \cr
#'     \code{features}   \cr \tab A data.frame with feature meta-data for rows of expression. \cr
#'     \code{samples}    \cr \tab A data.frame with observation meta-data for columns of expression. \cr
#'   }
#' @examples
#'   \dontrun{
#'     config <- list(dir_in=".", feature_file_in="feats.tsv", 
#'       sample_file_in="samps.tsv", data_file_in="exprs.tsv", log_file="")
#'     state <- f.read_data(config)
#'     exprs <- state$expression
#'     feats <- state$features
#'     samps <- state$samples
#'   }

f.read_data <- function(config) {

  f.log("reading data", config=config)
  
  file_path <- paste(config$dir_in, config$feature_file_in, sep="/")
  feats <- utils::read.table(file_path, header=T, sep="\t", quote="", as.is=T)
  
  file_path <- paste(config$dir_in, config$sample_file_in, sep="/")
  samps <- utils::read.table(file_path, header=T, sep="\t", quote="", as.is=T)
  
  file_path <- paste(config$dir_in, config$data_file_in, sep="/")
  exprs <- utils::read.table(file_path, header=T, sep="\t", quote="", as.is=T)
  exprs <- as.matrix(exprs)

  if(!(typeof(exprs) %in% "double")) {
    f.err("!(typeof(exprs) %in% 'double')", config=config)
  }
  
  state <- list(expression=exprs, features=feats, samples=samps)
  f.check_state(state, config)
  
  return(state)
}

###############################################################################

## set up types and levels of covariates; uses config$feat_id_col:

f.check_covariates <- function(state, config) {

  for(nom in c("n_samples_expr_col", "median_raw_col", "n_features_expr_col")) {
    if(nom %in% names(state$samples)) {
      f.err("nom %in% names(state$samples); nom:", nom, config=config)
    }
  }

  for(nom in c(config$obs_col, config$sample_col)) {
    if(!(nom %in% names(state$samples))) {
      f.err("!(nom %in% names(state$samples)); nom:", nom, config=config)
    }
  }

  if(!(config$feat_id_col %in% names(state$features))) {
    f.err("!(config$feat_id_col %in% names(state$features))", config=config)
  }
}

## uses config$sample_factors, config$obs_col, and config$sample_col:

f.subset_covariates <- function(state, config) {

  ## variables referred to in formula:
  vars <- as.character(config$frm)[2]
  vars <- gsub("[[:space:]]+", "", vars)
  vars <- gsub("[\\:\\*\\-]", "+", vars)
  vars <- unlist(strsplit(vars, split="\\+"))
  vars <- sort(unique(vars))

  ## make sure all needed variables in samps:
  if(!all(vars %in% names(state$samples))) {
    f.err("!all(vars %in% names(state$samples))", config=config)
  }
  
  if(!all(names(config$sample_factors) %in% vars)) {
    f.err("!all(names(config$sample_factors) %in% vars)", config=config)
  }
  
  f.msg("subsetting sample metadata", config=config)
  state$samples <- state$samples[, c(config$obs_col, config$sample_col, vars)]
  
  return(state)
}

## uses config$sample_factors

f.set_covariate_factor_levels <- function(state, config) {

  f.msg("setting factor levels", config=config)
  
  for(nom in names(config$sample_factors)) {
  
    ## check for potential misconfiguration first:
    if(!(nom %in% names(state$samples))) {
      f.err("!(nom %in% names(state$samples)); nom:", nom, config=config)
    }
    
    lvls1 <- config$sample_factors[[nom]]
    lvls2 <- sort(unique(as.character(state$samples[[nom]])))
    if(!all(lvls1 %in% lvls2)) {
      f.err("!all(lvls1 %in% lvls2) for nom:", nom, config=config)
    }
    if(!all(lvls2 %in% lvls1)) {
      f.err("!all(lvls2 %in% lvls1) for nom:", nom, config=config)
    }
    state$samples[[nom]] <- as.character(state$samples[[nom]])
    state$samples[[nom]] <- factor(
      state$samples[[nom]], 
      levels=config$sample_factors[[nom]]
    )
  }
  
  return(state)
}

#' Preprocess covariates
#' @description
#' `f.preprocess_covariates` checks covariates referred to in
#'   `config$frm` are in `state$samples`, subsets only the needed
#'   variables into `state$samples`, and sets covariate factor 
#'   levels according to `config$sample_factors`.
#' @details 
#'   Wrapper for `f.check_covariates()`, then `f.subset_covariates()`, 
#'     followed by `f.set_covariate_factor_levels()`.
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
#'     state <- list(expression=exprs, features=feats, samples=samps)
#'     config <- list(
#'       feat_id_col="gene_id", 
#'       sample_factors=list(age=c("Young", "Old"), gender=c("Male", "Female")),
#'       obs_col="observation_id",
#'       sample_col="sample_id",
#'       log_file=""
#'     )
#'     state <- f.preprocess_covariates(state, config)
#'     exprs <- state$expression
#'     feats <- state$features
#'     samps <- state$samples
#'   }

f.preprocess_covariates <- function(state, config) {

  f.log("preprocessing covariates", config=config)

  f.check_covariates(state, config)
  state <- f.subset_covariates(state, config)
  state <- f.set_covariate_factor_levels(state, config)
  
  return(state)
}

###############################################################################

#' Add filter statistics
#' @description
#' `f.add_filter_stats` adds filtering-related statistics to `state$features`, 
#'   and `state$samples`.
#' @details 
#'   Wrapper for `f.samples_per_feature()`, `f.feature_median_expression()`, 
#'     `f.features_per_sample()`. Also reports quantiles of distributions. 
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
#'     state <- list(expression=exprs, features=feats, samples=samps)
#'     config <- list(log_file="")
#'     state <- f.preprocess_covariates(state, config)
#'     exprs <- state$expression
#'     feats <- state$features
#'     samps <- state$samples
#'   }

f.add_filter_stats <- function(state, config) {

  n <- f.samples_per_feature(state, config)
  if(!all(names(n) == state$features[, config$feat_id_col, drop=T])) {
    f.err("!all(names(n) == state$features[, config$feat_id_col])", config=config)
  }
  state$features[, config$n_samples_expr_col] <- n

  m <- f.feature_median_expression(state, config)
  if(!all(names(m) == state$features[, config$feat_id_col, drop=T])) {
    f.err("!all(names(m) == state$features[, config$feat_id_col, drop=T])", config=config)
  }
  state$features[, config$median_raw_col] <- m

  n <- f.features_per_sample(state, config)
  if(!all(names(n) == state$samples[, config$obs_col, drop=T])) {
    f.err("!all(names(n) == state$samples[, config$obs_col])", config=config)
  } 
  state$samples[, config$n_features_expr_col] <- n
  
  n <- apply(state$expression, 1, function(v) sum(v > 0, na.rm=T))
  f.msg("samples per feature:", config=config)
  f.quantile(n, config, digits=0)

  n <- apply(state$expression, 2, function(v) sum(v > 0, na.rm=T))
  f.msg("features per sample", config=config)
  f.quantile(n, config, digits=0)
}

#' Prefilter data
#' @description
#' `f.prefilter` adds filtering-related statistics to `state$features`, 
#'   and `state$samples`.
#' @details 
#'   Wrapper for `f.samples_per_feature()`, `f.feature_median_expression()`, 
#'     `f.features_per_sample()`. Also reports quantiles of distributions. 
#' @param state List with elements formatted like the list returned by `f.read_data()`:
#'   \tabular{ll}{
#'     \code{expression} \cr \tab Numeric matrix with non-negative expression values. \cr
#'     \code{features}   \cr \tab A data.frame with feature meta-data for rows of expression. \cr
#'     \code{samples}    \cr \tab A data.frame with observation meta-data for columns of expression. \cr
#'   } 
#' @param config List with configuration values.
#' @param n_samples_min minimum number of samples per feature; numeric >= 1.
#' @param n_features_min minimum number of features per sample; numeric >= 1.
#' @return A list (the filtered state) with the following elements:
#'   \tabular{ll}{
#'     \code{expression} \cr \tab Numeric matrix with non-negative expression values. \cr
#'     \code{features}   \cr \tab A data.frame with feature meta-data for rows of expression. \cr
#'     \code{samples}    \cr \tab A data.frame with observation meta-data for columns of expression. \cr
#'   } 
#' @examples
#'   \dontrun{
#'     state <- list(expression=exprs, features=feats, samples=samps)
#'     config <- list(log_file="")
#'     state <- f.prefilter(state, config)
#'     exprs <- state$expression
#'     feats <- state$features
#'     samps <- state$samples
#'   }

f.prefilter <- function(state, config, n_samples_min=1, n_features_min=1) {
  
  f.log("prefilter features and samples", config=config)
  f.msg("before filtering features:", config=config)
  f.report_state(state, config)
  
  state <- f.filter_features(state, config, n_samples_min=n_samples_min)
  f.msg("after filtering features", config=config)
  f.report_state(state, config)

  state <- f.filter_samples(state, config, n_features_min=n_features_min)
  f.msg("after filtering samples", config=config)
  f.report_state(state, config)
  
  f.check_state(state, config)
  return(state)
}

#' Permute data
#' @description Permute 
#' @details If variable is NULL, uses config$permute_var instead. If variable is 
#'   NULL and config$permute_var == "", skips permutation (normal execution).
#' @param state List with elements formatted like the list returned by `f.read_data()`:
#'   \tabular{ll}{
#'     \code{expression} \cr \tab Numeric matrix with non-negative expression values. \cr
#'     \code{features}   \cr \tab A data.frame with feature meta-data for rows of expression. \cr
#'     \code{samples}    \cr \tab A data.frame with observation meta-data for columns of expression. \cr
#'   } 
#' @param config List with configuration values.
#' @param variable Character name of variable (column in samples) to permute.
#' @return A list (the filtered state) with the following elements:
#'   \tabular{ll}{
#'     \code{expression} \cr \tab Numeric matrix with non-negative expression values. \cr
#'     \code{features}   \cr \tab A data.frame with feature meta-data for rows of expression. \cr
#'     \code{samples}    \cr \tab A data.frame with observation meta-data for columns of expression. \cr
#'   } 
#' @examples
#'   \dontrun{
#'     state <- list(expression=exprs, samples=samps)
#'     config <- list(log_file="")
#'     state2 <- f.permute(state, config, permute_var="age")
#'     exprs2 <- state2$expression
#'     samps2 <- state2$samples
#'
#'     config <- list(log_file="")
#'     state2 <- f.permute(state, config, variable="age")
#'     exprs2 <- state2$expression
#'     samps2 <- state2$samples
#'   }

f.permute <- function(state, config, variable=NULL) {

  if(is.null(variable)) variable <- config$permute_var
  
  if(is.null(variable) || variable %in% "") {
    f.msg("skipping permutation", config=config)
  } else {
    if(!(variable %in% colnames(state$samples))) {
      f.err("!(variable %in% colnames(state$samples))", config=config)
    }
    f.msg("permuting", variable, config=config)
    tmp <- state$samples[
      !duplicated(state$samples[, config$sample_col]), 
      c(config$sample_col, variable), 
      drop=T
    ]
    rownames(tmp) <- tmp[, config$sample_col]
    tmp[, variable] <- sample(tmp[, variable, drop=T], nrow(tmp), replace=F)
    state$samples[, variable] <- tmp[state$samples[, config$sample_col], variable, drop=T]
  } 

  f.check_state(state, config)
  return(state)
}

## SAVE 2: prepped
## f.check_state(state, config)
## f.report_state(state, config)
## f.save_state(state, config, "2.prepped")

