#' Read data
#' @description
#' Reads expression data, feature metadata, and observation metadata files.
#' @details 
#'   Expression data from \code{config$data_file_in} in \code{config$dir_in}.
#'   Feature metadata from \code{config$feature_file_in} in \code{config$dir_in}.
#'   Observation metadata from \code{config$sample_file_in} in \code{config$dir_in}.
#'   See documentation for \code{h0testr::new_config()} 
#'     for more detailed description of configuration parameters. 
#' @param config List with configuration values. Requires the following keys:
#'   \tabular{ll}{
#'     \code{dir_in}    \cr \tab Path to directory where input files located (character). \cr
#'     \code{feature_file_in} \cr \tab Name of file with feature metadata; assumed in \code{dir_in}. \cr
#'     \code{sample_file_in}  \cr \tab Name of file with observation metadata; assumed in \code{dir_in}. \cr
#'     \code{data_file_in}    \cr \tab Name of file with signal data; assumed in \code{dir_in}. \cr
#'     \code{feat_id_col}     \cr \tab Name of column in \code{feature_file_in} that corresponds to rows of \code{data_file_in}. \cr
#'     \code{gene_id_col}     \cr \tab Name of column in \code{feature_file_in} with unique genegroup or proteingroup ids. \cr
#'     \code{obs_id_col}      \cr \tab Name of column in \code{sample_file_in} that corresponds to columns of \code{data_file_in}. \cr
#'     \code{sample_id_col}   \cr \tab Name of column in \code{sample_file_in} with unique sample labels. \cr
#'   }
#' @return A list (the initial state) with the following elements:
#'   \tabular{ll}{
#'     \code{expression} \cr \tab Numeric matrix with non-negative expression values. \cr
#'     \code{features}   \cr \tab A data.frame with feature meta-data for rows of expression. \cr
#'     \code{samples}    \cr \tab A data.frame with observation meta-data for columns of expression. \cr
#'   }
#' @examples
#' config <- h0testr::new_config()
#' config$save_state <- FALSE
#' config$dir_in <- system.file("extdata", package="h0testr")  ## where example data 
#' config$feature_file_in <- "features.tsv"
#' config$sample_file_in <- "samples.tsv"
#' config$data_file_in <- "expression.tsv" 
#' config$feat_id_col <- "feature_id"
#' config$obs_id_col <- "observation_id"
#' 
#' state <- h0testr::read_data(config)
#' 
#' names(state)
#' print(state$features)
#' print(state$samples)
#' print(state$expression[1:6, 1:6])

read_data <- function(config) {

  check_config(config)
  
  f.log("reading data", config=config)
  
  file_path <- paste(config$dir_in, config$feature_file_in, sep="/")
  feats <- utils::read.table(file_path, header=T, sep="\t", quote="", as.is=T)
  
  file_path <- paste(config$dir_in, config$sample_file_in, sep="/")
  samps <- utils::read.table(file_path, header=T, sep="\t", quote="", as.is=T)
  
  file_path <- paste(config$dir_in, config$data_file_in, sep="/")
  exprs <- utils::read.table(file_path, header=T, sep="\t", quote="", as.is=T)
  exprs <- as.matrix(exprs)
  
  if(!(typeof(exprs) %in% c("double", "integer"))) {
    f.err("read_data: !(typeof(exprs) %in% c('double', 'integer'))", 
      "typof(exprs):", typeof(exprs), config=config)
  }
  
  state <- list(expression=exprs, features=feats, samples=samps)
  
  return(state)
}

###############################################################################

## set up types and levels of covariates; uses config$feat_id_col:

f.check_parameters <- function(state, config, initialized=F, minimal=F) {
  
  if(minimal) {
    reqd_params <- c("obs_id_col", "sample_id_col", 
      "feat_id_col", "gene_id_col", 
      "frm", "test_term", "sample_factors")
  } else {
    reqd_params <- c("n_samples_expr_col", "median_raw_col", 
      "n_features_expr_col", "obs_id_col", "sample_id_col", "feat_id_col", 
      "gene_id_col", "frm", "test_term", "sample_factors")
  }
  
  for(param in reqd_params) {
    if(!(param %in% names(config))) {
      f.err("f.check_parameters: !(param %in% names(config)), for param:", 
        param, "; names(config):", names(config), config=config)
    }
  }
  
  if(!is.null(config$save_state) && config$save_state) {
    reqd_params <- c("dir_out", "data_mid_out", "feature_mid_out", 
      "sample_mid_out", "suffix_out")
    for(param in reqd_params) {
      if(!(param %in% names(config))) {
        f.err("f.check_parameters: !(param %in% names(config)), for param:", 
          param, "; names(config):", names(config), config=config)
      }
    }
  }
  
  if(!initialized) {
    noms <- c(config$n_samples_expr_col, config$median_raw_col, config$n_features_expr_col)
    for(nom in noms) {
      if(nom %in% names(state$samples)) {
        f.err("f.check_parameters: nom %in% names(state$samples); nom:", 
          nom, "; names(state$samples):", names(state$samples), config=config)
      }
    }
  }
  
  for(nom in c(config$obs_id_col, config$sample_id_col)) {
    if(!(nom %in% names(state$samples))) {
      f.err("f.check_parameters: !(nom %in% names(state$samples)); nom:", 
        nom, "; names(state$samples):", names(state$samples), config=config)
    }
  }
  
  for(nom in c(config$feat_id_col, config$gene_id_col)) {
    if(!(nom %in% names(state$features))) {
      f.err("f.check_parameters: !(nom %in% names(state$features)); nom:",
        nom, "; names(state$features):", names(state$features), config=config)
    }
  }
  
  if(any(duplicated(state$features[[config$feat_id_col]]))) {
    f.err("f.check_parameters: any(duplicated(state$features[[config$feat_id_col]]));",
      "duplicated:", 
      state$features[[config$feat_id_col]][duplicated(state$features[[config$feat_id_col]])], 
      config=config)
  }
  
  if(any(duplicated(state$samples[[config$obs_id_col]]))) {
    f.err("f.check_parameters: any(duplicated(state$samples[[config$obs_id_col]]))",
      "duplicated:", 
      state$samples[[config$obs_id_col]][duplicated(state$samples[[config$obs_id_col]])],
      config=config)
  }
      
  return(TRUE)
}

## uses config$frm, config$sample_factors, config$obs_id_col, and config$sample_id_col:

f.subset_covariates <- function(state, config) {
  
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
    f.err("f.subset_covariates: !all(vars %in% names(state$samples)); vars:", 
      vars, config=config)
  }
  
  if(!all(names(config$sample_factors) %in% vars)) {
    f.err("f.subset_covariates: !all(names(config$sample_factors) %in% config$frm); vars:", 
      vars, config=config)
  }
  
  f.msg("subsetting sample metadata", config=config)
  state$samples <- state$samples[, unique(c(config$obs_id_col, config$sample_id_col, vars))]
  
  return(state)
}

## uses config$sample_factors

f.set_covariate_factor_levels <- function(state, config) {
  
  f.msg("setting factor levels", config=config)
  
  for(nom in names(config$sample_factors)) {
    
    ## check for potential misconfiguration first:
    if(!(nom %in% names(state$samples))) {
      f.err("f.set_covariate_factor_levels: !(nom %in% names(state$samples)); nom:", 
        nom, config=config)
    }
    
    lvls1 <- config$sample_factors[[nom]]
    lvls2 <- sort(unique(as.character(state$samples[[nom]])))
    if(!all(lvls1 %in% lvls2)) {
      f.err("f.set_covariate_factor_levels: !all(lvls1 %in% lvls2) for nom:", 
        nom, config=config)
    }
    if(!all(lvls2 %in% lvls1)) {
      f.err("f.set_covariate_factor_levels: !all(lvls2 %in% lvls1) for nom:", 
        nom, config=config)
    }
    state$samples[[nom]] <- as.character(state$samples[[nom]])
    state$samples[[nom]] <- factor(
      state$samples[[nom]], 
      levels=config$sample_factors[[nom]]
    )
  }
  
  return(state)
}

#' Preprocess covariates and parameters in config
#' @description
#'   Check \code{config} vs. \code{state}, check covariates in \code{config$frm},
#'     and set factor levels.
#' @details 
#'   Checks to make sure columns specified in \code{config} parameters are 
#'     found in \code{state$features} and \code{state$samples}.
#'   Check covariates referred to in \code{config$frm} are found in 
#'     \code{state$samples}, subsets only the needed variables into 
#'     \code{state$samples}, and sets covariate factor levels 
#'     according to \code{config$sample_factors}.
#'   Flow is:
#'     \tabular{l}{
#'       1. \code{check_config()}. 
#'       2. Subset covariates in \code{config$frm}. 
#'       3. Check feat_col and obs_col.
#'       4. Subset covariates of interest.
#'       5. Set covariate factor levels.
#'       6. Return sub-table of ANOVA results corresponding to 
#'            \code{config$test_term}.
#'     }
#'   If \code{initialized=FALSE}, then checks if \code{state$features} has 
#'     columns with names in 
#'     \code{c(config$n_samples_expr_col, config$median_raw_col)}; and if 
#'     \code{state$samples} has a column named \code{config$n_features_expr_col}. 
#'     If either is \code{TRUE}, results in error.
#'   See documentation for \code{h0testr::new_config()} 
#'     for more detailed description of configuration parameters. 
#' @param state List with elements formatted like the list returned by \code{f.read_data()}:
#'   \tabular{ll}{
#'     \code{expression} \cr \tab Numeric matrix with non-negative expression values. \cr
#'     \code{features}   \cr \tab A data.frame with feature meta-data for rows of expression. \cr
#'     \code{samples}    \cr \tab A data.frame with observation meta-data for columns of expression. \cr
#'   } 
#' @param config List with configuration values. Keys 
#'   \code{c("n_samples_expr_col", "median_raw_col", "n_features_expr_col")} 
#'     not needed if \code{minimal=TRUE}:
#'   \tabular{ll}{
#'     \code{obs_id_col}           \cr \tab Column in \code{state$samples} corresponding to \code{colnames(state$expression)}. \cr
#'     \code{sample_id_col}        \cr \tab Column in \code{state$samples} with unique sample labels. \cr
#'     \code{feat_id_col}          \cr \tab Column in \code{state$features} corresponding to \code{rownames(state$expression)}. \cr
#'     \code{gene_id_col}          \cr \tab Column in \code{state$features} with unique gene/protein group ids. \cr
#'     \code{frm}                  \cr \tab Formula object specifying formula to be fit. \cr
#'     \code{test_term}            \cr \tab Term (character) in \code{config$frm} to test for significance. \cr
#'     \code{sample_factors}       \cr \tab List with levels of factor variables in \code{config$frm}. \cr
#'     \code{n_samples_expr_col}   \cr \tab Column in \code{state$features} that corresponds to columns of \code{data_file_in}. \cr
#'     \code{median_raw_col}       \cr \tab Column in \code{state$features} that corresponds to columns of \code{data_file_in}. \cr
#'     \code{n_features_expr_col}  \cr \tab Column in \code{state$samples} that corresponds to columns of \code{data_file_in}. \cr
#'   }
#' @param initialized Logical scalar indicating if \code{state} has already 
#'   has filter statistics initialized. 
#' @param minimal Logical scalar indicating whether only minimal set of 
#'   parameters needed for formula processing should be required.
#' @return A list with the following two elements:
#'   \tabular{ll}{
#'     \code{state}   \cr \tab List with elements \code{c("expression", "features", "samples")}.
#'     \code{config}  \cr \tab List with configuration settings.
#'   }
#'   The element \code{state} is a list with the following three elements:
#'   \tabular{ll}{
#'     \code{expression} \cr \tab Numeric matrix with non-negative expression values. \cr
#'     \code{features}   \cr \tab A data.frame with feature meta-data for rows of expression. \cr
#'     \code{samples}    \cr \tab A data.frame with observation meta-data for columns of expression. \cr
#'   } 
#' @examples
#' exprs <- h0testr::sim1(n_obs=6, n_feats=12)$mat
#' feats <- data.frame(feature_id=rownames(exprs), gene_id=rownames(exprs))
#' samps <- data.frame(observation_id=colnames(exprs), 
#'   condition=c(rep("ctl", 3), rep("trt", 3)))
#' state <- list(expression=exprs, features=feats, samples=samps)
#'
#' ## minimal:
#' config <- list(
#'   obs_id_col="observation_id",
#'   sample_id_col="observation_id",
#'   feat_id_col="feature_id",
#'   gene_id_col="gene_id",
#'   frm=~condition,
#'   test_term="condition",
#'   sample_factors=list(condition=c("ctl", "trt"))
#' )
#' out <- h0testr::initialize(state, config, minimal=TRUE)
#' print(out$state)
#' str(out$config)
#'
#' config <- list(
#'   obs_id_col="observation_id",
#'   sample_id_col="observation_id",
#'   feat_id_col="feature_id",
#'   gene_id_col="gene_id",
#'   frm=~condition,
#'   test_term="condition", 
#'   sample_factors=list(condition=c("ctl", "trt")),
#'   n_samples_expr_col="n_samples_exprs",
#'   median_raw_col="median_raw",
#'   n_features_expr_col="n_features_exprs"
#' )
#' out <- h0testr::initialize(state, config)
#' print(out$state)
#' str(out$config)

initialize <- function(state, config, initialized=F, minimal=F) {
  
  f.log("initializing", config=config)
  check_config(config)  
  f.check_parameters(state, config, initialized=initialized, minimal=minimal)
  
  if(is.null(config$feat_col) || config$feat_col %in% "") {
    config$feat_col <- config$feat_id_col  ## as soon as confirm feat_id_col exists
  }
  
  if(is.null(config$obs_col) || config$obs_col %in% "") {
    config$obs_col <- config$obs_id_col    ## as soon as confirm obs_id_col exists
  }
  
  state <- f.subset_covariates(state, config)
  state <- f.set_covariate_factor_levels(state, config)
  
  return(list(state=state, config=config))
}

###############################################################################

#' Permute data
#' @description Permute observation covariate.
#' @details If \code{variable} is \code{NULL}, uses \code{config$permute_var} 
#'   instead. If \code{variable} is \code{NULL} and 
#'   \code{config$permute_var == ""}, skips permutation (normal execution).
#'   See documentation for \code{h0testr::new_config()} 
#'     for more detailed description of configuration parameters. 
#' @param state List with elements formatted like the list returned by 
#'   \code{read_data()}:
#'   \tabular{ll}{
#'     \code{expression} \cr \tab Numeric matrix with non-negative expression values. \cr
#'     \code{features}   \cr \tab A data.frame with feature meta-data for rows of expression. \cr
#'     \code{samples}    \cr \tab A data.frame with observation meta-data for columns of expression. \cr
#'   } 
#' @param config List with configuration values. Requires the following keys:
#'   \tabular{ll}{
#'     \code{obs_col}        \cr \tab Name of column in \code{feature_file_in} that corresponds to columns of \code{data_file_in}. \cr
#'     \code{sample_id_col}  \cr \tab Name of column in \code{sample_file_in} with unique sample identifiers. \cr
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
#' exprs <- h0testr::sim1(n_obs=6, n_feats=12)$mat
#' feats <- data.frame(feature_id=rownames(exprs))
#' samps <- data.frame(observation_id=colnames(exprs), age=c(rep("young", 3), rep("old", 3)))
#' state <- list(expression=exprs, features=feats, samples=samps)
#' print(state)
#' config <- list(feat_col="feature_id", obs_col="observation_id", sample_id_col="observation_id", permute_var="age")
#' state2 <- h0testr::permute(state, config)
#' print(state2)

permute <- function(state, config, variable=NULL) {

  check_config(config)
  
  if(is.null(variable)) variable <- config$permute_var
  
  if(is.null(variable) || variable %in% "") {
    f.log_block("skipping permutation", config=config)
  } else {
    if(!(variable %in% colnames(state$samples))) {
      f.err("permute: !(variable %in% colnames(state$samples))", config=config)
    }
    f.log_block("permuting", variable, config=config)
    tmp <- state$samples[
      !duplicated(state$samples[, config$sample_id_col, drop=T]), 
      c(config$sample_id_col, variable), 
      drop=T
    ]
    rownames(tmp) <- tmp[, config$sample_id_col, drop=T]
    tmp[, variable] <- sample(tmp[[variable]], nrow(tmp), replace=F)
    state$samples[, variable] <- tmp[state$samples[[config$sample_id_col]], variable, drop=T]
  } 
  
  f.check_state(state, config)
  return(state)
}

#' Load data
#' @description 
#'   Load data and metadata from files, format, and save initial copies.
#' @details Loads data from files specified in \code{config}. Prefilter uninformative rows
#'   and columns. Permute variable if requested. Save final copies.
#'   See documentation for \code{h0testr::new_config()} 
#'     for more detailed description of configuration parameters. 
#' @param config List with configuration values. Requires the following keys:
#'   \tabular{ll}{
#'     \code{dir_in}          \cr \tab Path to directory where input files located (character). \cr
#'     \code{feature_file_in} \cr \tab Name of file with feature metadata; assumed in \code{dir_in}. \cr
#'     \code{sample_file_in}  \cr \tab Name of file with observation metadata; assumed in \code{dir_in}. \cr
#'     \code{data_file_in}    \cr \tab Name of file with signal data; assumed in \code{dir_in}. \cr
#'     \code{feat_id_col}     \cr \tab Name of column in \code{feature_file_in} that corresponds to rows of \code{data_file_in}. \cr
#'     \code{gene_id_col}     \cr \tab Name of column in \code{feature_file_in} with unique gene group or protein group labels. \cr
#'     \code{obs_id_col}      \cr \tab Name of column in \code{sample_file_in} that corresponds to columns of \code{data_file_in}. \cr
#'     \code{sample_id_col}   \cr \tab Name of column in \code{sample_file_in} with unique sample labels. \cr
#'   }
#' @return A list (the initial state) with the following elements:
#'   \tabular{ll}{
#'     \code{expression} \cr \tab Numeric matrix with expression values. \cr
#'     \code{features}   \cr \tab A data.frame with feature meta-data for rows of expression. \cr
#'     \code{samples}    \cr \tab A data.frame with observation meta-data for columns of expression. \cr
#'   } 
#' @examples
#' config <- h0testr::new_config()     ## defaults
#' config$save_state <- FALSE            ## default is TRUE
#' config$dir_in <- system.file("extdata", package="h0testr")  ## where example data 
#' config$feature_file_in <- "features.tsv"
#' config$sample_file_in <- "samples.tsv"
#' config$data_file_in <- "expression.tsv" 
#' config$feat_id_col <- "feature_id"
#' config$gene_id_col <- "feature_id"
#' config$obs_id_col <- "observation_id"
#' config$sample_id_col <- "observation_id"
#' config$frm <- ~condition
#' config$test_term <- "condition"
#' config$test_method <- "trend"
#' config$sample_factors <- list(condition=c("placebo", "drug"))
#' 
#' output <- h0testr::load_data(config)
#' 
#' names(output)
#' names(output$state)
#' output$state$expression
#' output$state$samples
#' output$state$features

load_data <- function(config) {
  
  report_config(config)
  
  state <- read_data(config)
  out <- initialize(state, config)
  state <- out$state
  config <- out$config
  f.check_state(state, config)
  
  state <- add_filter_stats(state, config)
  f.check_state(state, config)
  f.report_state(state, config)
  f.save_state(state, config, prefix="1.initial")
  
  state <- prefilter(state, config)
  state <- permute(state, config)
  
  f.check_state(state, config)
  f.report_state(state, config)
  f.save_state(state, config, prefix="2.prepped")
  
  return(list(state=state, config=config))
}
