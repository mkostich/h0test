#' Read data
#' @description
#' Reads expression data, feature metadata, and observation metadata files.
#' @details 
#'   Expression data from \code{config$data_file_in} in \code{config$dir_in}.
#'   Feature metadata from \code{config$feature_file_in} in \code{config$dir_in}.
#'   Observation metadata from \code{config$sample_file_in} in \code{config$dir_in}.
#' @param config List with configuration values. Requires the following keys:
#'   \tabular{ll}{
#'     \code{dir_in}    \cr \tab Path to directory where input files located (character). \cr
#'     \code{feature_file_in} \cr \tab Name of file with feature metadata; assumed in \code{dir_in}. \cr
#'     \code{sample_file_in}  \cr \tab Name of file with observation metadata; assumed in \code{dir_in}. \cr
#'     \code{data_file_in}    \cr \tab Name of file with signal data; assumed in \code{dir_in}. \cr
#'     \code{feat_id_col}     \cr \tab Name of column (character) in \code{feature_file_in} that corresponds to rows of \code{data_file_in}. \cr
#'     \code{obs_id_col}      \cr \tab Name of column (character) in \code{sample_file_in} that corresponds to columns of \code{data_file_in}. \cr
#'     \code{sample_id_col}   \cr \tab Name of column (character) in \code{sample_file_in} with unique sample labels. \cr
#'   }
#' @return A list (the initial state) with the following elements:
#'   \tabular{ll}{
#'     \code{expression} \cr \tab Numeric matrix with non-negative expression values. \cr
#'     \code{features}   \cr \tab A data.frame with feature meta-data for rows of expression. \cr
#'     \code{samples}    \cr \tab A data.frame with observation meta-data for columns of expression. \cr
#'   }
#' @examples
#' config <- h0testr::f.new_config()
#' config$dir_in <- system.file("extdata", package="h0testr")  ## where example data 
#' config$feature_file_in <- "features.tsv"
#' config$sample_file_in <- "samples.tsv"
#' config$data_file_in <- "expression.tsv" 
#' config$feat_id_col <- "feature_id"
#' config$obs_id_col <- "observation_id"
#' 
#' state <- h0testr::f.read_data(config)
#' 
#' names(state)
#' print(state$features)
#' print(state$samples)
#' print(state$expression[1:6, 1:6])

f.read_data <- function(config) {

  f.check_config(config)
  
  f.log("reading data", config=config)
  
  file_path <- paste(config$dir_in, config$feature_file_in, sep="/")
  feats <- utils::read.table(file_path, header=T, sep="\t", quote="", as.is=T)
  
  file_path <- paste(config$dir_in, config$sample_file_in, sep="/")
  samps <- utils::read.table(file_path, header=T, sep="\t", quote="", as.is=T)
  
  file_path <- paste(config$dir_in, config$data_file_in, sep="/")
  exprs <- utils::read.table(file_path, header=T, sep="\t", quote="", as.is=T)
  exprs <- as.matrix(exprs)
  
  if(!(typeof(exprs) %in% "double")) {
    f.err("f.read_data: !(typeof(exprs) %in% 'double')", config=config)
  }
  
  f.log("checking initial state", config=config)
  state <- list(expression=exprs, features=feats, samples=samps)
  f.check_state(state, config)
  
  return(state)
}

###############################################################################

## set up types and levels of covariates; uses config$feat_id_col:

f.check_covariates <- function(state, config, initialized=F) {

  f.check_config(config)
  
  reqd_params <- c("n_samples_expr_col", "median_raw_col", 
    "n_features_expr_col", "obs_id_col", "sample_id_col", "feat_id_col", 
    "sample_factors")
    
  for(param in reqd_params) {
    if(!(param %in% names(config))) {
      f.err("f.check_covariates: !(param %in% names(config)), for param:", 
        param, config=config)
    }
  }
  
  if(!is.null(config$save_state) && config$save_state) {
    reqd_params <- c("dir_out", "data_mid_out", "feature_mid_out", 
      "sample_mid_out", "suffix_out")
    for(param in reqd_params) {
      if(!(param %in% names(config))) {
        f.err("f.check_covariates: !(param %in% names(config)), for param:", 
          param, config=config)
      }
    }
  }
  
  if(!initialized) {
    noms <- c(config$n_samples_expr_col, config$median_raw_col, config$n_features_expr_col)
    for(nom in noms) {
      if(nom %in% names(state$samples)) {
        f.err("f.check_covariates: nom %in% names(state$samples); nom:", 
          nom, config=config)
      }
    }
  }
  
  for(nom in c(config$obs_id_col, config$sample_id_col)) {
    if(!(nom %in% names(state$samples))) {
      f.err("f.check_covariates: !(nom %in% names(state$samples)); nom:", 
        nom, config=config)
    }
  }
  
  if(!(config$feat_id_col %in% names(state$features))) {
    f.err("f.check_covariates: !(config$feat_id_col %in% names(state$features))", config=config)
  }
  
  if(any(duplicated(state$features[[config$feat_id_col]]))) {
    f.err("f.check_covariates: any(duplicated(state$features[[config$feat_id_col]]))")
  }
  
  if(any(duplicated(state$samples[[config$obs_id_col]]))) {
    f.err("f.check_covariates: any(duplicated(state$samples[[config$obs_id_col]]))")
  }
}

## uses config$sample_factors, config$obs_id_col, and config$sample_id_col:

f.subset_covariates <- function(state, config) {

  f.check_config(config)
  
  ## variables referred to in formula:
  vars <- as.character(config$frm)
  vars <- gsub("[[:space:]]+", "", vars)
  vars <- gsub("[\\~\\:\\*\\-]", "+", vars)
  vars <- unlist(strsplit(vars, split="\\+"))
  vars <- sort(unique(vars))
  n <- nchar(vars)
  n[is.na(n)] <- 0
  vars <- vars[n > 0]
  
  ## make sure all needed variables in samps:
  if(!all(vars %in% names(state$samples))) {
    f.err("f.subset_covariates: !all(vars %in% names(state$samples)); vars:", vars, config=config)
  }
  
  if(!all(names(config$sample_factors) %in% vars)) {
    f.err("f.subset_covariates: !all(names(config$sample_factors) %in% config$frm); vars:", vars, config=config)
  }
  
  f.msg("subsetting sample metadata", config=config)
  state$samples <- state$samples[, unique(c(config$obs_id_col, config$sample_id_col, vars))]
  
  return(state)
}

## uses config$sample_factors

f.set_covariate_factor_levels <- function(state, config) {

  f.check_config(config)
  
  f.msg("setting factor levels", config=config)
  
  for(nom in names(config$sample_factors)) {
    
    ## check for potential misconfiguration first:
    if(!(nom %in% names(state$samples))) {
      f.err("f.set_covariate_factor_levels: !(nom %in% names(state$samples)); nom:", nom, config=config)
    }
    
    lvls1 <- config$sample_factors[[nom]]
    lvls2 <- sort(unique(as.character(state$samples[[nom]])))
    if(!all(lvls1 %in% lvls2)) {
      f.err("f.set_covariate_factor_levels: !all(lvls1 %in% lvls2) for nom:", nom, config=config)
    }
    if(!all(lvls2 %in% lvls1)) {
      f.err("f.set_covariate_factor_levels: !all(lvls2 %in% lvls1) for nom:", nom, config=config)
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
#'   Check covariates referred to in
#'     \code{config$frm} are in \code{state$samples}, subsets only the needed
#'     variables into \code{state$samples}, and sets covariate factor 
#'     levels according to \code{config$sample_factors}.
#' @details 
#'   Wrapper for \code{f.check_covariates()}, then \code{f.subset_covariates()}, 
#'     followed by \code{f.set_covariate_factor_levels()}. If \code{initialized=FALSE}, 
#'     then checks if \code{state$features} has columns with names in 
#'     \code{c(config$n_samples_expr_col, config$median_raw_col)}; and if 
#'     \code{state$samples} has a column named \code{config$n_features_expr_col}. If either
#'     is \code{TRUE}, results in error.
#' @param state List with elements formatted like the list returned by \code{f.read_data()}:
#'   \tabular{ll}{
#'     \code{expression} \cr \tab Numeric matrix with non-negative expression values. \cr
#'     \code{features}   \cr \tab A data.frame with feature meta-data for rows of expression. \cr
#'     \code{samples}    \cr \tab A data.frame with observation meta-data for columns of expression. \cr
#'   } 
#' @param config List with configuration values. Requires the following keys:
#'   \tabular{ll}{
#'     \code{sample_factors}  \cr \tab List specifying levels of factor variables in \code{config$frm} (see examples). \cr
#'     \code{feat_id_col}     \cr \tab Name of column (character) in \code{feature_file_in} that corresponds to rows of \code{data_file_in}. \cr
#'     \code{obs_id_col}      \cr \tab Name of column (character) in \code{sample_file_in} that corresponds to columns of \code{data_file_in}. \cr
#'     \code{sample_id_col}   \cr \tab Name of column (character) in \code{sample_file_in} with unique sample labels. \cr
#'   }
#' @param initialized Logical scalar indicating if \code{state} has already has filter statistics initialized. 
#' @return A list (the processed state) with the following elements:
#'   \tabular{ll}{
#'     \code{expression} \cr \tab Numeric matrix with non-negative expression values. \cr
#'     \code{features}   \cr \tab A data.frame with feature meta-data for rows of expression. \cr
#'     \code{samples}    \cr \tab A data.frame with observation meta-data for columns of expression. \cr
#'   } 
#' @examples
#' exprs <- h0testr::f.sim1(n_obs=6, n_feats=12)$mat
#' feats <- data.frame(feature_id=rownames(exprs))
#' samps <- data.frame(observation_id=colnames(exprs), 
#'   condition=c(rep("ctl", 3), rep("trt", 3)))
#' state <- list(expression=exprs, features=feats, samples=samps)
#' config <- h0testr::f.new_config()
#' config$feat_id_col <- "feature_id"
#' config$obs_id_col <- "observation_id"
#' config$sample_id_col <- "observation_id"
#' config$frm <- ~condition
#' config$test_term <- "condition"
#' config$sample_factors=list(condition=c("ctl", "trt"))
#' out <- h0testr::f.preprocess_covariates(state, config)
#' print(out)

f.preprocess_covariates <- function(state, config, initialized=F) {
  
  f.log("preprocessing covariates", config=config)
  
  f.check_covariates(state, config, initialized=initialized)
  config$obs_col <- config$obs_id_col   ## as soon as confirm obs_id_col exists
  
  state <- f.subset_covariates(state, config)
  state <- f.set_covariate_factor_levels(state, config)
  
  return(list(state=state, config=config))
}

###############################################################################

#' Add filter statistics
#' @description
#'   Adds filtering-related statistics to \code{state$features}, 
#'     and \code{state$samples}.
#' @details 
#'   Wrapper for \code{f.samples_per_feature()}, \code{f.feature_median_expression()}, 
#'     \code{f.features_per_sample()}. Also reports quantiles of distributions. 
#' @param state List with elements formatted like the list returned by \code{f.read_data()}:
#'   \tabular{ll}{
#'     \code{expression} \cr \tab Numeric matrix with non-negative expression values. \cr
#'     \code{features}   \cr \tab A data.frame with feature meta-data for rows of expression. \cr
#'     \code{samples}    \cr \tab A data.frame with observation meta-data for columns of expression. \cr
#'   } 
#' @param config List with configuration values. Requires the following keys:
#'   \tabular{ll}{
#'     \code{feat_id_col}  \cr \tab Name of column (character) in \code{feature_file_in} that corresponds to rows of \code{data_file_in}. \cr
#'     \code{obs_id_col}   \cr \tab Name of column (character) in \code{sample_file_in} that corresponds to columns of \code{expression}. \cr
#'   }
#' @return A list (the processed state) with the following elements:
#'   \tabular{ll}{
#'     \code{expression} \cr \tab Numeric matrix with non-negative expression values. \cr
#'     \code{features}   \cr \tab A data.frame with feature meta-data for rows of expression. \cr
#'     \code{samples}    \cr \tab A data.frame with observation meta-data for columns of expression. \cr
#'   } 
#' @examples
#' exprs <- h0testr::f.sim1(n_obs=6, n_feats=8)$mat
#' feats <- data.frame(feature_id=rownames(exprs))
#' samps <- data.frame(observation_id=colnames(exprs))
#' state <- list(expression=exprs, features=feats, samples=samps)
#' config <- h0testr::f.new_config()
#' state <- h0testr::f.add_filter_stats(state, config)
#' print(state)

f.add_filter_stats <- function(state, config) {

  f.check_config(config)
  
  n <- f.samples_per_feature(state, config)
  if(!all(names(n) == state$features[, config$feat_id_col, drop=T])) {
    f.err("f.add_filter_stats: !all(names(n) == state$features[, config$feat_id_col])", config=config)
  }
  state$features[, config$n_samples_expr_col] <- n
  
  m <- f.feature_median_expression(state, config)
  if(!all(names(m) == state$features[, config$feat_id_col, drop=T])) {
    f.err("f.add_filter_stats: !all(names(m) == state$features[, config$feat_id_col, drop=T])", config=config)
  }
  state$features[, config$median_raw_col] <- m
  
  n <- f.features_per_sample(state, config)
  if(is.null(config$obs_col) || config$obs_col %in% "") {
    config$obs_col <- config$obs_id_col
  }
  if(!all(names(n) == state$samples[, config$obs_col, drop=T])) {
    f.err("f.add_filter_stats: !all(names(n) == state$samples[, config$obs_col])", config=config)
  } 
  state$samples[, config$n_features_expr_col] <- n
  
  n <- apply(state$expression, 1, function(v) sum(v > 0, na.rm=T))
  f.msg("samples per feature:", config=config)
  f.quantile(n, config, digits=0)
  
  n <- apply(state$expression, 2, function(v) sum(v > 0, na.rm=T))
  f.msg("features per sample", config=config)
  f.quantile(n, config, digits=0)
  
  return(state)
}

#' Prefilter data
#' @description
#'   Adds filtering-related statistics to \code{state$features}, 
#'     and \code{state$samples}.
#' @details 
#'   Wrapper for \code{f.samples_per_feature()}, \code{f.feature_median_expression()}, 
#'     \code{f.features_per_sample()}. Also reports quantiles of distributions. 
#' @param state List with elements formatted like the list returned by \code{f.read_data()}:
#'   \tabular{ll}{
#'     \code{expression} \cr \tab Numeric matrix with non-negative expression values. \cr
#'     \code{features}   \cr \tab A data.frame with feature meta-data for rows of expression. \cr
#'     \code{samples}    \cr \tab A data.frame with observation meta-data for columns of expression. \cr
#'   } 
#' @param config List with configuration values. Requires the following keys:
#'   \tabular{ll}{
#'     \code{feat_id_col}  \cr \tab Name of column (character) in \code{feature_file_in} that corresponds to rows of \code{data_file_in}. \cr
#'     \code{obs_id_col}   \cr \tab Name of column (character) in \code{sample_file_in} that corresponds to columns of \code{expression}. \cr
#'   }
#' @param n_samples_min minimum number of samples per feature; numeric >= 2.
#' @param n_features_min minimum number of features per sample; numeric >= 2.
#' @return A list (the filtered state) with the following elements:
#'   \tabular{ll}{
#'     \code{expression} \cr \tab Numeric matrix with non-negative expression values. \cr
#'     \code{features}   \cr \tab A data.frame with feature meta-data for rows of expression. \cr
#'     \code{samples}    \cr \tab A data.frame with observation meta-data for columns of expression. \cr
#'   } 
#' @examples
#' set.seed(101)
#' exprs <- h0testr::f.sim1(n_obs=6, n_feats=12, mcar_p=0.75)$mat
#' feats <- data.frame(feature_id=rownames(exprs))
#' samps <- data.frame(observation_id=colnames(exprs))
#' state <- list(expression=exprs, features=feats, samples=samps)
#' config <- list(feat_id_col="feature_id", obs_id_col="observation_id")
#' state2 <- h0testr::f.prefilter(state, config)
#' print(state2)

f.prefilter <- function(state, config, n_samples_min=2, n_features_min=2) {

  f.check_config(config)
  
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
#' @description Permute observation covariate 
#' @details If variable is \code{NULL}, uses \code{config$permute_var} instead. If variable is 
#'   \code{NULL} and \code{config$permute_var == ""}, skips permutation (normal execution).
#' @param state List with elements formatted like the list returned by \code{f.read_data()}:
#'   \tabular{ll}{
#'     \code{expression} \cr \tab Numeric matrix with non-negative expression values. \cr
#'     \code{features}   \cr \tab A data.frame with feature meta-data for rows of expression. \cr
#'     \code{samples}    \cr \tab A data.frame with observation meta-data for columns of expression. \cr
#'   } 
#' @param config List with configuration values. Requires the following keys:
#'   \tabular{ll}{
#'     \code{feat_id_col}    \cr \tab Name of column (character) in \code{feature_file_in} that corresponds to rows of \code{data_file_in}. \cr
#'     \code{sample_id_col}  \cr \tab Name of column (character) in \code{sample_file_in} with unique sample identifiers. \cr
#'   }
#' @param variable Character name of variable (column in samples) to permute.
#' @return A list (the permuted state) with the following elements:
#'   \tabular{ll}{
#'     \code{expression} \cr \tab Numeric matrix with non-negative expression values. \cr
#'     \code{features}   \cr \tab A data.frame with feature meta-data for rows of expression. \cr
#'     \code{samples}    \cr \tab A data.frame with observation meta-data for columns of expression. \cr
#'   } 
#' @examples
#' set.seed(101)
#' exprs <- h0testr::f.sim1(n_obs=6, n_feats=12)$mat
#' feats <- data.frame(feature_id=rownames(exprs))
#' samps <- data.frame(observation_id=colnames(exprs), age=c(rep("young", 3), rep("old", 3)))
#' state <- list(expression=exprs, features=feats, samples=samps)
#' config <- list(feat_id_col="feature_id", obs_id_col="observation_id", 
#'   sample_id_col="observation_id", permute_var="age")
#' state2 <- h0testr::f.permute(state, config)
#' print(state)
#' print(state2)

f.permute <- function(state, config, variable=NULL) {

  f.check_config(config)
  
  if(is.null(variable)) variable <- config$permute_var
  
  if(is.null(variable) || variable %in% "") {
    f.msg("skipping permutation", config=config)
  } else {
    if(!(variable %in% colnames(state$samples))) {
      f.err("f.permute: !(variable %in% colnames(state$samples))", config=config)
    }
    f.msg("permuting", variable, config=config)
    tmp <- state$samples[
      !duplicated(state$samples[, config$sample_id_col]), 
      c(config$sample_id_col, variable), 
      drop=T
    ]
    rownames(tmp) <- tmp[, config$sample_id_col]
    tmp[, variable] <- sample(tmp[, variable, drop=T], nrow(tmp), replace=F)
    state$samples[, variable] <- tmp[state$samples[, config$sample_id_col], variable, drop=T]
  } 
  
  f.check_state(state, config)
  return(state)
}

#' Load data
#' @description 
#'   Load data and metadata from files, format, and save initial copies 
#' @details Loads data from files specified in \code{config}. Prefilter uninformative rows
#'   and columns. Permute variable if requested. Save final copies.
#' @param config List with configuration values. Requires the following keys:
#'   \tabular{ll}{
#'     \code{dir_in}          \cr \tab Path to directory where input files located (character). \cr
#'     \code{feature_file_in} \cr \tab Name of file with feature metadata; assumed in \code{dir_in}. \cr
#'     \code{sample_file_in}  \cr \tab Name of file with observation metadata; assumed in \code{dir_in}. \cr
#'     \code{data_file_in}    \cr \tab Name of file with signal data; assumed in \code{dir_in}. \cr
#'     \code{feat_id_col}     \cr \tab Name of column (character) in \code{feature_file_in} that corresponds to rows of \code{data_file_in}. \cr
#'     \code{obs_id_col}      \cr \tab Name of column (character) in \code{sample_file_in} that corresponds to columns of \code{data_file_in}. \cr
#'     \code{sample_id_col}   \cr \tab Name of column (character) in \code{sample_file_in} with unique sample labels. \cr
#'   }
#' @return A list (the initial state) with the following elements:
#'   \tabular{ll}{
#'     \code{expression} \cr \tab Numeric matrix with expression values. \cr
#'     \code{features}   \cr \tab A data.frame with feature meta-data for rows of expression. \cr
#'     \code{samples}    \cr \tab A data.frame with observation meta-data for columns of expression. \cr
#'   } 
#' @examples
#' config <- h0testr::f.new_config()
#' config$dir_in <- system.file("extdata", package="h0testr")  ## where example data 
#' config$feature_file_in <- "features.tsv"
#' config$sample_file_in <- "samples.tsv"
#' config$data_file_in <- "expression.tsv" 
#' config$feat_id_col <- "feature_id"
#' config$obs_id_col <- "observation_id"
#' config$sample_id_col <- "observation_id"
#' config$frm <- ~condition
#' config$test_term <- "condition"
#' config$test_method <- "trend"
#' config$sample_factors <- list(condition=c("placebo", "drug"))
#' config$save_state <- FALSE
#' 
#' output <- h0testr::f.load_data(config)
#' 
#' names(output)
#' names(output$state)
#' output$state$expression
#' output$state$samples
#' output$state$features

f.load_data <- function(config) {
  
  f.report_config(config)
  
  state <- f.read_data(config)
  out <- f.preprocess_covariates(state, config)
  state <- out$state
  config <- out$config
  
  state <- f.add_filter_stats(state, config)
  f.check_state(state, config)
  f.report_state(state, config)
  f.save_state(state, config, prefix="1.initial")
  
  state <- f.prefilter(state, config)
  state <- f.permute(state, config)
  
  f.check_state(state, config)
  f.report_state(state, config)
  f.save_state(state, config, prefix="2.prepped")
  
  return(list(state=state, config=config))
}

