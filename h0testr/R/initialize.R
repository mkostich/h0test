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
      "frm", "test_term", "reference_levels")
  } else {
    reqd_params <- c("n_samples_expr_col", "median_raw_col", 
      "n_features_expr_col", "obs_id_col", "sample_id_col", "feat_id_col", 
      "gene_id_col", "frm", "test_term", "reference_levels")
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
    for(nom in c(config$n_samples_expr_col, config$median_raw_col)) {
      if(nom %in% names(state$features)) {
        f.err("f.check_parameters: nom %in% names(state$features); nom:",
          nom, "; names(state$features):", names(state$features), config=config)
      }
    }
    ## NOTE: n_features_expr_col not a required param when minimal, so may be
    ##   unset here; NULL %in% x is logical(0), which if() cannot handle:
    if(length(config$n_features_expr_col) %in% 1 &&
      config$n_features_expr_col %in% names(state$samples)) {
      f.err("f.check_parameters: n_features_expr_col %in% names(state$samples); nom:",
        config$n_features_expr_col,
        "; names(state$samples):", names(state$samples), config=config)
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

## Establish the missingness contract on state$expression: NA is the only
##   indicator of a missing value everywhere downstream. Raw (untransformed)
##   expression cannot be negative, and a raw zero means the feature was not
##   detected, so zeros become NA. When config$is_log_transformed is TRUE the
##   input is already on a log-like scale, in which case zeros and negative
##   values are legitimate (e.g. vsn output, or log2 of an intensity below 1)
##   and no conversion is done; the caller is then responsible for having marked
##   non-detections as NA. config$normalization_method is deliberately not
##   consulted: "none" says only that no normalization is wanted, which is a
##   separate question from what scale the values are already on, and reading a
##   scale out of it left the two able to disagree. Fills in
##   config$is_log_transformed when it was unset, so that every downstream step
##   finds an answer there:

f.zeros_to_na <- function(state, config) {

  if(!is.matrix(state$expression)) {
    f.err("f.zeros_to_na: !is.matrix(state$expression);",
      "class(state$expression):", class(state$expression), config=config)
  }

  ## unset means raw, matching the new_config() default:

  raw <- !isTRUE(config$is_log_transformed)

  if(!raw) {
    f.msg("f.zeros_to_na: config$is_log_transformed is TRUE, so treating",
      "state$expression as already transformed;",
      "zeros left as they are and negative values allowed;", "\n",
      "NA is the only indicator of a missing value", config=config)
    config$is_log_transformed <- TRUE
    return(list(state=state, config=config))
  }

  i_neg <- !is.na(state$expression) & state$expression < 0

  if(any(i_neg)) {

    rnom <- rownames(state$expression)
    cnom <- colnames(state$expression)
    if(is.null(rnom)) rnom <- as.character(1:nrow(state$expression))
    if(is.null(cnom)) cnom <- as.character(1:ncol(state$expression))

    idxs <- utils::head(which(i_neg), 5)
    rr <- ((idxs - 1) %% nrow(state$expression)) + 1
    cc <- ((idxs - 1) %/% nrow(state$expression)) + 1

    f.err("f.zeros_to_na: negative values in state$expression;",
      "raw expression values cannot be negative;", "\n",
      "if the input is already log transformed, set config$is_log_transformed",
      "to TRUE;", "\n",
      "negative values:", sum(i_neg), "of", length(i_neg), ";",
      "first offenders (feature, observation, value):", "\n",
      paste(rnom[rr], cnom[cc], state$expression[idxs], sep=", "),
      config=config)
  }

  i_zero <- !is.na(state$expression) & state$expression == 0
  state$expression[i_zero] <- NA

  f.msg("f.zeros_to_na: converted", sum(i_zero), "zero values to NA;",
    "NA is now the only indicator of a missing value", config=config)
  config$is_log_transformed <- FALSE

  return(list(state=state, config=config))
}

## uses config$frm, config$reference_levels, config$obs_id_col, and config$sample_id_col:

f.subset_covariates <- function(state, config) {
  
  ## variables referred to in formula:
  vars <- sort(unique(f.parse_frm(config$frm, config)$vars))

  ## make sure all needed variables in samps:
  if(!all(vars %in% names(state$samples))) {
    f.err("f.subset_covariates: !all(vars %in% names(state$samples)); vars:", 
      vars, "; names(state$samples):", names(state$samples), config=config)
  }
  
  if(!all(names(config$reference_levels) %in% vars)) {
    f.err(
      "f.subset_covariates:",
      "!all(names(config$reference_levels) %in% config$frm)", "\n",
      "vars:", vars, "\n",
      "names(config$reference_levels):", names(config$reference_levels),
      config=config
    )
  }
  
  f.msg("subsetting sample metadata", config=config)
  state$samples <- state$samples[, unique(c(config$obs_id_col, config$sample_id_col, vars))]
  
  return(state)
}

## uses config$reference_levels and config$frm; sets levels of the factor
##   covariates in config$frm, putting the declared reference level first and
##   the remaining levels in sorted order, and reports the levels and reference
##   level of each factor covariate, as well as the distribution of each numeric
##   (continuous) covariate, so that a misclassified covariate is visible in the
##   log. Optional types from f.covariate_types():

f.set_covariate_factor_levels <- function(state, config, types=NULL) {

  f.msg("setting factor levels", config=config)

  if(is.null(types)) types <- f.covariate_types(state, config)

  for(nom in names(types)) {

    ## check for potential misconfiguration first:
    if(!(nom %in% names(state$samples))) {
      f.err("f.set_covariate_factor_levels: !(nom %in% names(state$samples)); nom:",
        nom, "; names(state$samples):", names(state$samples), config=config)
    }

    if(types[nom] %in% "numeric") {
      f.msg("covariate", nom, ": numeric (continuous); distribution:", config=config)
      f.quantile(state$samples[[nom]], config)
      next
    }

    ## NOTE: config$reference_levels is an atomic vector, so [[ throws on a
    ##   name that is not present, rather than returning NULL:

    ref1 <- NULL
    if(nom %in% names(config$reference_levels)) {
      ref1 <- config$reference_levels[[nom]]
    }

    if(is.null(ref1)) {

      ## logical covariates are ordered FALSE, TRUE, and a covariate that is
      ##   already a factor carries its own level ordering; neither of these is
      ##   locale dependent, so neither needs to be declared:

      state$samples[[nom]] <- factor(state$samples[[nom]])

    } else {

      if(!(length(ref1) %in% 1 && is.character(ref1))) {
        f.err("f.set_covariate_factor_levels: config$reference_levels[[nom]]",
          "not scalar character, for nom:", nom, "\n", "value:", ref1,
          "; class:", class(ref1), config=config)
      }

      ## only the reference level is declared, so the remaining levels are
      ##   sorted; their order does not affect the reference level, but it does
      ##   affect the order in which coefficients are reported:

      lvls2 <- sort(unique(as.character(state$samples[[nom]])))
      if(!(ref1 %in% lvls2)) {
        f.err("f.set_covariate_factor_levels: reference level", ref1,
          "declared in config$reference_levels is not among the values of",
          "covariate", nom, ";", "\n", "its values:", lvls2, config=config)
      }
      state$samples[[nom]] <- as.character(state$samples[[nom]])
      state$samples[[nom]] <- factor(state$samples[[nom]],
        levels=c(ref1, setdiff(lvls2, ref1)))
    }

    lvls <- levels(state$samples[[nom]])
    f.msg("covariate", nom, ": factor; levels:", lvls, "; reference:", lvls[1],
      config=config)
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
#'     according to \code{config$reference_levels}.
#'   Each covariate in \code{config$frm} is classified as a factor or as
#'     numeric (continuous): variables named in \code{config$reference_levels}
#'     as well as logical variables are factors; other numeric variables are
#'     continuous; a character variable that is not named in
#'     \code{config$reference_levels} is an error, since its reference level,
#'     and hence the meaning of its coefficients, would otherwise be set by
#'     locale dependent sorting of its values. The classification is returned in
#'     \code{config$covariate_types}. Missing, non-finite, and constant
#'     covariate values are errors. A continuous covariate with
#'     \code{config$n_distinct_numeric_warn} or fewer distinct values is
#'     warned about, as it may be a miscoded factor.
#'   \code{config$reference_levels} declares only the reference (first) level of
#'     each factor covariate; the remaining levels are sorted. The resulting
#'     level ordering of each factor covariate is logged, and returned in
#'     \code{config$factor_levels}.
#'   Establishes the missingness contract on \code{state$expression}: \code{NA}
#'     is the only indicator of a missing value everywhere downstream. Raw
#'     expression values cannot be negative, and a raw zero means the feature was
#'     not detected, so a negative value is an error and zeros are converted to
#'     \code{NA}. When \code{config$is_log_transformed} is \code{TRUE} the input
#'     is already on a log-like scale, in which case zeros and negative values
#'     are legitimate and no conversion is done; the caller is then responsible
#'     for having marked non-detections as \code{NA}. Either way,
#'     \code{config$is_log_transformed} is filled in if it was unset, so that
#'     every downstream step finds the scale recorded there. Skipped when
#'     \code{initialized=TRUE}, since an already initialized state may
#'     legitimately hold negative (transformed) values.
#'   Flow is:
#'     \tabular{l}{
#'       1. \code{check_config()}. \cr
#'       2. Settle missingness: raw zeros to \code{NA}; negative raw value is an error. \cr
#'       3. Subset covariates in \code{config$frm}. \cr
#'       4. Check feat_col and obs_col. \cr
#'       5. Subset covariates of interest. \cr
#'       6. Classify covariates and check their values. \cr
#'       7. Set covariate factor levels. \cr
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
#'     \code{reference_levels}     \cr \tab Named character vector with the reference level of each factor variable in \code{config$frm}. \cr
#'     \code{n_distinct_numeric_warn} \cr \tab Warn if continuous variable in \code{config$frm} has this few distinct values. \cr
#'     \code{is_log_transformed}   \cr \tab Optional logical; whether the input is already on a log-like scale. Default \code{FALSE}, meaning raw, so zeros become \code{NA} and negative values are an error. \cr
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
#'   reference_levels=c(condition="ctl")
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
#'   reference_levels=c(condition="ctl"),
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

  ## settle missingness before anything else reads the matrix; skipped when the
  ##   state was already initialized, since by then the values may legitimately
  ##   be negative (transformed):

  if(!initialized) {
    out <- f.zeros_to_na(state, config)
    state <- out$state
    config <- out$config
  }

  ## the dependent variable is always the expression values of one feature, so
  ##   any dependent given in config$frm carries no information; drop it here,
  ##   so it is reported once, rather than by every downstream call:

  parsed <- f.parse_frm(config$frm, config)
  if(parsed$two_sided) {
    f.msg(
      "initialize:",
      "ignoring dependent variable on left-hand side of config$frm;", "\n",
      "the dependent is always the expression values of one feature.", "\n",
      "Setting config$frm to:", parsed$frm,
      config=config
    )
    config$frm <- parsed$frm
  }

  if(is.null(config$feat_col) || config$feat_col %in% "") {
    config$feat_col <- config$feat_id_col  ## as soon as confirm feat_id_col exists
  }
  
  if(is.null(config$obs_col) || config$obs_col %in% "") {
    config$obs_col <- config$obs_id_col    ## as soon as confirm obs_id_col exists
  }
  
  state$features[[config$feat_id_col]] <- as.character(state$features[[config$feat_id_col]])
  state$features[[config$gene_id_col]] <- as.character(state$features[[config$gene_id_col]])
  state$samples[[config$obs_id_col]] <- as.character(state$samples[[config$obs_id_col]])
  state$samples[[config$sample_id_col]] <- as.character(state$samples[[config$sample_id_col]])
  
  if(!any(duplicated(state$features[[config$gene_id_col]]))) {
    f.msg(
      "initialize:",
      "config$gene_id_col has no duplicates in state$features", "\n", 
      "No point in aggregating features", "\n",
      "Setting config$gene_id_col to config$feat_id_col:", 
      config$feat_id_col, "\n", 
      "Setting config$feature_aggregation to 'none'.",
      config=config
    )
    config$feat_col <- config$gene_id_col <- config$feat_id_col
    config$feature_aggregation <- "none"
  }
  
  state <- f.subset_covariates(state, config)

  ## classify the covariates in config$frm as factor or numeric (continuous)
  ##   once, here, and check their values before anything is fit; any cached
  ##   classification is discarded, so that initialize() is authoritative:

  config$covariate_types <- NULL
  types <- f.covariate_types(state, config)
  f.check_covariate_values(state, config, types=types)
  config$covariate_types <- types

  state <- f.set_covariate_factor_levels(state, config, types=types)

  ## record the resolved levels of each factor covariate, in the order used to
  ##   build the design matrix, since config$reference_levels only declares the
  ##   first of them:

  config$factor_levels <- list()
  for(nom in names(types)[types %in% "factor"]) {
    config$factor_levels[[nom]] <- levels(state$samples[[nom]])
  }

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
      f.err(
        "permute: !(variable %in% colnames(state$samples))", "\n",
        "variable:", variable, "\n", 
        "colnames(state$samples):", colnames(state$samples), 
        config=config
      )
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
#' config$reference_levels <- c(condition="placebo")
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
