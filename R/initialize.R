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
#'     \code{dir_in}    \tab Path to directory where input files located (character). \cr
#'     \code{feature_file_in} \tab Name of file with feature metadata; assumed in \code{dir_in}. \cr
#'     \code{sample_file_in}  \tab Name of file with observation metadata; assumed in \code{dir_in}. \cr
#'     \code{data_file_in}    \tab Name of file with signal data; assumed in \code{dir_in}. \cr
#'     \code{feat_id_col}     \tab Name of column in \code{feature_file_in} that corresponds to rows of \code{data_file_in}. \cr
#'     \code{gene_id_col}     \tab Name of column in \code{feature_file_in} with unique genegroup or proteingroup ids. \cr
#'     \code{obs_id_col}      \tab Name of column in \code{sample_file_in} that corresponds to columns of \code{data_file_in}. \cr
#'     \code{sample_id_col}   \tab Name of column in \code{sample_file_in} with unique sample labels. \cr
#'   }
#' @return A list (the initial state) with the following elements:
#'   \tabular{ll}{
#'     \code{expression} \tab Numeric matrix with non-negative expression values. \cr
#'     \code{features}   \tab A data.frame with feature meta-data for rows of expression. \cr
#'     \code{samples}    \tab A data.frame with observation meta-data for columns of expression. \cr
#'   }
#' @examples
#' config <- h0testr::new_config()
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
#' @export

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
      "typeof(exprs):", typeof(exprs), config=config)
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
    ## n_features_expr_col not required param when minimal, so can be
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

  ## metadata and state$expression matched by position everywhere downstream:

  f.check_expr_matrix(state, config, "f.check_parameters")

  f.check_dimnames(state, config, "feat_id_col", fn_name="f.check_parameters")
  f.check_dimnames(state, config, "obs_id_col", fn_name="f.check_parameters")

  return(TRUE)
}

## for state$expression: ensure NA is only indicator of missing value:

f.zeros_to_na <- function(state, config) {

  if(!is.matrix(state$expression)) {
    f.err("f.zeros_to_na: !is.matrix(state$expression);",
      "class(state$expression):", class(state$expression), config=config)
  }

  ## unset means raw, matching new_config() default:

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
  
  noms <- unique(c(config$obs_id_col, config$sample_id_col, vars))
  state$samples <- state$samples[, noms, drop=FALSE]

  return(state)
}

## uses config$reference_levels and config$frm; sets levels of factor
##   covariates in config$frm:

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
    
    state$samples[[nom]] <- factor(as.character(state$samples[[nom]]),
      levels=f.covariate_levels(state$samples[[nom]], nom, config,
        caller="f.set_covariate_factor_levels"))

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
#'     covariate values are errors, as are blank (empty or whitespace only)
#'     values of a factor covariate: \code{utils::read.table()}, which
#'     \code{h0testr::read_data()} uses, reads an empty field in a character
#'     column as an empty string rather than as \code{NA}, so an empty cell in
#'     the samples file would otherwise become a factor level of its own,
#'     silently adding a group made of the observations whose annotation is
#'     missing. A continuous covariate with
#'     \code{config$n_distinct_numeric_warn} or fewer distinct values is
#'     warned about, as it may be a miscoded factor.
#'   \code{config$reference_levels} declares only the reference (first) level of
#'     each factor covariate; the remaining levels are sorted. The resulting
#'     level ordering of each factor covariate is logged, and returned in
#'     \code{config$factor_levels}.
#'   Establishes uniform representation of missing values: \code{NA}
#'     is the only indicator of a missing value everywhere downstream. Raw
#'     expression values cannot be negative, and a raw zero means the feature was
#'     not detected, so a negative value is an error and zeros are converted to
#'     \code{NA}. When \code{config$is_log_transformed} is \code{TRUE} the input
#'     is already on a log-like scale, in which case zeros and negative values
#'     are legitimate and no conversion is done; the caller is then responsible
#'     for having marked non-detections as \code{NA}. Conversion is skipped
#'     when \code{initialized=TRUE}, since an already initialized state may
#'     legitimately hold negative (transformed) values, but
#'     \code{config$is_log_transformed} is filled in if it was unset,
#'     so that every downstream step finds scale recorded there. Unset means
#'     raw.
#'   Checks that the metadata lines up with \code{state$expression}: 
#'     \code{state$features} must have one row per row of the matrix and 
#'     \code{state$samples} one row per column, and \code{config$feat_id_col} and 
#'     \code{config$obs_id_col} must agree with \code{rownames(state$expression)} and
#'     \code{colnames(state$expression)}. The two are matched by position
#'     everywhere downstream.
#'   Flow is:
#'     \tabular{l}{
#'       1. \code{check_config()}. \cr
#'       2. Check parameters and ids, and that metadata lines up with the matrix. \cr
#'       3. Settle missingness: raw zeros to \code{NA}; negative raw value is an error. \cr
#'       4. Drop any dependent variable given on left-hand side of \code{config$frm}. \cr
#'       5. Set \code{config$feat_col} and \code{config$obs_col}, if unset. \cr
#'       6. Subset covariates of interest. \cr
#'       7. Classify covariates and check their values. \cr
#'       8. Set covariate factor levels. \cr
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
#'     \code{expression} \tab Numeric matrix with non-negative expression values. \cr
#'     \code{features}   \tab A data.frame with feature meta-data for rows of expression. \cr
#'     \code{samples}    \tab A data.frame with observation meta-data for columns of expression. \cr
#'   } 
#' @param config List with configuration values. Keys 
#'   \code{c("n_samples_expr_col", "median_raw_col", "n_features_expr_col")} 
#'     not needed if \code{minimal=TRUE}:
#'   \tabular{ll}{
#'     \code{obs_id_col}           \tab Column in \code{state$samples} corresponding to \code{colnames(state$expression)}. \cr
#'     \code{sample_id_col}        \tab Column in \code{state$samples} with unique sample labels. \cr
#'     \code{feat_id_col}          \tab Column in \code{state$features} corresponding to \code{rownames(state$expression)}. \cr
#'     \code{gene_id_col}          \tab Column in \code{state$features} with unique gene/protein group ids. \cr
#'     \code{frm}                  \tab Formula object specifying formula to be fit. \cr
#'     \code{test_term}            \tab Term (character) in \code{config$frm} to test for significance. \cr
#'     \code{contrast}      \tab Weighted sum (character scalar) of coefficients of \code{config$frm} to test instead of \code{config$test_term}; "" for none. \cr
#'     \code{reference_levels}     \tab Named character vector with the reference level of each factor variable in \code{config$frm}. \cr
#'     \code{n_distinct_numeric_warn} \tab Warn if continuous variable in \code{config$frm} has this few distinct values. \cr
#'     \code{is_log_transformed}   \tab Optional logical; whether the input is already on a log-like scale. Default \code{FALSE}, meaning raw, so zeros become \code{NA} and negative values are an error. \cr
#'     \code{n_samples_expr_col}   \tab Column in \code{state$features} that corresponds to columns of \code{data_file_in}. \cr
#'     \code{median_raw_col}       \tab Column in \code{state$features} that corresponds to columns of \code{data_file_in}. \cr
#'     \code{n_features_expr_col}  \tab Column in \code{state$samples} that corresponds to columns of \code{data_file_in}. \cr
#'   }
#' @param initialized Logical scalar indicating if \code{state} has already
#'   had filter statistics initialized. \code{TRUE} relaxes checking only: the
#'   three statistic columns are no longer required to be absent, and 
#'   missingness conversion is skipped, since values may by then legitimately
#'   be negative. \code{config$is_log_transformed} is still filled in if unset,
#'   and alignment of the metadata with \code{state$expression} is checked
#'   either way.
#' @param minimal Logical scalar indicating whether only minimal set of 
#'   parameters needed for formula processing should be required.
#' @return A list with the following two elements:
#'   \tabular{ll}{
#'     \code{state}   \tab List with elements \code{c("expression", "features", "samples")}. \cr
#'     \code{config}  \tab List with configuration settings. \cr
#'   }
#'   The element \code{state} is a list with the following three elements:
#'   \tabular{ll}{
#'     \code{expression} \tab Numeric matrix with non-negative expression values. \cr
#'     \code{features}   \tab A data.frame with feature meta-data for rows of expression. \cr
#'     \code{samples}    \tab A data.frame with observation meta-data for columns of expression. \cr
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
#' out <- h0testr::init_state(state, config, minimal=TRUE)
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
#' out <- h0testr::init_state(state, config)
#' print(out$state)
#' str(out$config)
#' @export

init_state <- function(state, config, initialized=F, minimal=F) {
  
  f.log("initializing", config=config)
  check_config(config)
  f.check_parameters(state, config, initialized=initialized, minimal=minimal)

  if(!initialized) {

    out <- f.zeros_to_na(state, config)
    state <- out$state
    config <- out$config

  } else if(!isTRUE(config$is_log_transformed)) {

    if(is.null(config$is_log_transformed)) {
      f.msg("init_state: initialized=TRUE, so state$expression is left as it",
        "is: no zeros converted to NA and no check for negative values;", "\n",
        "config$is_log_transformed was unset, and is recorded as FALSE (raw),",
        "matching the h0testr::new_config() default;", "\n",
        "set it to TRUE before calling if the values are already on a log-like",
        "scale, as they are once h0testr::normalize() has run",
        config=config)
    }

    config$is_log_transformed <- FALSE
  }

  ## dependent variable is always expression values of one feature:

  parsed <- f.parse_frm(config$frm, config)
  if(parsed$two_sided) {
    f.msg(
      "init_state:",
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
      "init_state:",
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
  config$covariate_types <- NULL
  types <- f.covariate_types(state, config)
  f.check_covariate_values(state, config, types=types)
  config$covariate_types <- types

  ## init_state() authoritative:
  config$factor_levels <- NULL

  state <- f.set_covariate_factor_levels(state, config, types=types)

  ## record resolved levels of each factor covariate, in same order used to
  ##   build design matrix:

  config$factor_levels <- list()
  for(nom in names(types)[types %in% "factor"]) {
    config$factor_levels[[nom]] <- levels(state$samples[[nom]])
  }

  return(list(state=state, config=config))
}

###############################################################################

## shuffle x within each level of g, by index; sample() on a length-1 numeric would
##   return permutation of seq_len(x) rather than x itself:

f.shuffle_within <- function(x, g) {
  i <- unsplit(lapply(split(seq_along(x), g), function(k) k[sample.int(length(k))]), g)
  return(x[i])
}

## shuffle `variable` across unique sample ids, optionally within strata, and map 
##   result back onto observations, so tech reps of one sample keep one label:

f.permute_col <- function(samples, config, variable, within="") {

  id <- config$sample_id_col
  cols <- unique(c(id, variable, if(nzchar(within)) within else NULL))
  tmp <- samples[!duplicated(samples[[id]]), cols, drop=F]
  rownames(tmp) <- tmp[[id]]

  tmp[[variable]] <- if(nzchar(within)) {
    f.shuffle_within(tmp[[variable]], tmp[[within]])
  } else tmp[[variable]][sample.int(nrow(tmp))]

  samples[[variable]] <- tmp[samples[[id]], variable]
  return(samples)
}

## variance of test-term coefficient up to sigma^2, from design alone; NA if the
##   term is absent from the formula or shuffle makes design singular:

f.perm_var <- function(samples, frm, term) {

  mm <- try(stats::model.matrix(frm, data=samples), silent=T)
  if(inherits(mm, "try-error")) return(NA_real_)

  k <- match(term, attr(stats::terms(frm), "term.labels"))
  if(is.na(k)) return(NA_real_)
  j <- which(attr(mm, "assign") %in% k)
  if(!length(j)) return(NA_real_)

  v <- try(solve(crossprod(mm)), silent=T)
  if(inherits(v, "try-error")) return(NA_real_)
  return(mean(diag(v)[j]))
}

## how many strata hold more than one level of the test term; a stratum holding one level
##   contributes an identity permutation and so no null at all:

f.n_informative <- function(samples, term, within) {
  n <- tapply(samples[[term]], samples[[within]], function(v) length(unique(v)))
  return(sum(n > 1, na.rm=T))
}

#' Check that a permutation scheme preserves design
#' @description
#'   Reports whether shuffling \code{config$test_term} degrades precision with which
#'     that term is estimated, and if so which column of \code{samples} to permute within
#'     to avoid it.
#' @details
#'   Unresricted permutation is only valid where the observations are exchangeable. In a blocked or
#'     paired design they are not: shuffling the test term across blocks breaks the
#'     pairing, so the term is estimated from fewer informative units, its standard error
#'     inflates, and permuted results carry systematically fewer hits than they should.
#'     Where those permuted hits are used to estimate an FDR, as in
#'     \code{h0testr::tune_check()}, the FDR comes out too low.
#'   The statistic is the variance of the test-term coefficient up to \code{sigma^2},
#'     \code{(X'X)^-1}, which needs the design matrix only and no expression data. It is
#'     computed for the design as given and for \code{n_shuffles} shuffles of it; the
#'     inflation is the ratio of the median shuffled value to the unshuffled one. A ratio
#'     near 1 means the shuffle costs nothing.
#'   Suggestions are drawn from the variables of \code{config$frm} only, since a blocking
#'     factor that is not in the model is not one this check can price. A candidate is
#'     suggested when permuting within it brings the inflation back under
#'     \code{config$permute_inflation_max} and leaves at least two strata holding more
#'     than one level of the test term.
#'   See documentation for \code{h0testr::new_config()}
#'     for more detailed description of configuration parameters.
#' @param config List with configuration values. Requires the following keys:
#'   \tabular{ll}{
#'     \code{frm}                   \tab Formula with the variable of interest and covariates. \cr
#'     \code{test_term}             \tab Term in \code{frm} that is tested and permuted. \cr
#'     \code{sample_id_col}         \tab Name of column in \code{samples} with sample identifiers. \cr
#'     \code{permute_within}        \tab Strata column, or \code{""} for a free shuffle. \cr
#'     \code{permute_inflation_max} \tab Largest inflation ratio treated as acceptable. \cr
#'   }
#' @param samples A data.frame of observation meta-data, like \code{state$samples}.
#' @param variable Character name of the variable that will be shuffled. \code{NULL} uses
#'   \code{config$permute_var}, falling back to \code{config$test_term} where that is
#'   \code{""}. The variance priced is always that of \code{config$test_term}.
#' @param n_shuffles Number of shuffles to median over. Default 20.
#' @return A list with the following elements:
#'   \tabular{ll}{
#'     \code{inflation}     \tab Ratio of median shuffled to unshuffled coefficient variance (numeric). \cr
#'     \code{threshold}     \tab \code{config$permute_inflation_max} (numeric). \cr
#'     \code{within}        \tab Strata column the ratio was computed for (character). \cr
#'     \code{ok}            \tab TRUE if \code{inflation} is at or under \code{threshold} (logical). \cr
#'     \code{suggest}       \tab Columns that would bring the ratio under it (character). \cr
#'     \code{n_informative} \tab Strata holding more than one level of the test term (integer). \cr
#'   }
#' @examples
#' ## a paired design: one case, two conditions, is not freely exchangeable
#' samps <- data.frame(
#'   sample_id=paste0("s", 1:12),
#'   case=rep(paste0("c", 1:6), each=2),
#'   grp=rep(c("ctl", "trt"), 6))
#' config <- h0testr::new_config()
#' config$frm <- ~grp+case
#' config$test_term <- "grp"
#' config$sample_id_col <- "sample_id"
#' set.seed(101)
#' out <- h0testr::check_permutation(config, samps)
#' out$inflation                        ## well above 1: a free shuffle costs precision
#' out$suggest                          ## "case"
#' @export

check_permutation <- function(config, samples, variable=NULL, n_shuffles=20) {

  check_config(config)

  term <- config$test_term
  within <- if(is.null(config$permute_within)) "" else config$permute_within
  thr <- if(is.null(config$permute_inflation_max)) 1.2 else config$permute_inflation_max

  out <- list(inflation=NA_real_, threshold=thr, within=within, ok=T,
    suggest=character(0), n_informative=NA_integer_)

  ## with $contrast set in place of $test_term there is no single term to price:
  if(is.null(term) || !nzchar(term)) return(out)

  if(is.null(variable)) variable <- config$permute_var
  if(is.null(variable) || !nzchar(variable)) variable <- term

  vars <- all.vars(config$frm)
  if(!all(c(vars, variable, config$sample_id_col) %in% colnames(samples))) return(out)

  v0 <- f.perm_var(samples, config$frm, term)
  if(is.na(v0) || v0 <= 0) return(out)

  f.ratio <- function(g) {
    v <- replicate(n_shuffles,
      f.perm_var(f.permute_col(samples, config, variable, g), config$frm, term))
    return(stats::median(v, na.rm=T) / v0)
  }

  out$inflation <- f.ratio(within)
  if(nzchar(within)) out$n_informative <- f.n_informative(samples, variable, within)
  out$ok <- !is.na(out$inflation) && out$inflation <= thr
  if(out$ok) return(out)

  for(nom in setdiff(vars, variable)) {
    if(!(is.character(samples[[nom]]) || is.factor(samples[[nom]]))) next
    if(f.n_informative(samples, variable, nom) < 2) next
    r <- f.ratio(nom)
    if(!is.na(r) && r <= thr) out$suggest <- c(out$suggest, nom)
  }

  return(out)
}

###############################################################################

#' Permute data
#' @description Permute observation covariate.
#' @details If \code{variable} is \code{NULL}, uses \code{config$permute_var}
#'   instead. If \code{variable} is \code{NULL} and
#'   \code{config$permute_var == ""}, skips permutation (normal execution).
#'   The shuffle is over unique \code{config$sample_id_col}, so technical replicates of one
#'     sample keep a single label.
#'   \code{config$permute_within} names a column of strata to shuffle within, leaving the
#'     blocking structure of a paired or blocked design intact; \code{""}, the default,
#'     shuffles freely over all samples. Before shuffling, \code{h0testr::check_permutation()}
#'     prices the scheme and the call fails if it inflates the variance of the
#'     \code{config$test_term} coefficient by more than \code{config$permute_inflation_max},
#'     naming the column that would fix it. Set \code{config$permute_force} to permute
#'     regardless.
#'   See documentation for \code{h0testr::new_config()}
#'     for more detailed description of configuration parameters.
#' @param state List with elements formatted like the list returned by
#'   \code{read_data()}:
#'   \tabular{ll}{
#'     \code{expression} \tab Numeric matrix with non-negative expression values. \cr
#'     \code{features}   \tab A data.frame with feature meta-data for rows of expression. \cr
#'     \code{samples}    \tab A data.frame with observation meta-data for columns of expression. \cr
#'   }
#' @param config List with configuration values. Requires the following keys:
#'   \tabular{ll}{
#'     \code{obs_col}               \tab Name of column in \code{feature_file_in} that corresponds to columns of \code{data_file_in}. \cr
#'     \code{sample_id_col}         \tab Name of column in \code{sample_file_in} with unique sample identifiers. \cr
#'     \code{permute_within}        \tab Strata column to permute within, or \code{""} for a free shuffle. \cr
#'     \code{permute_inflation_max} \tab Largest acceptable inflation of the test-term coefficient variance. \cr
#'     \code{permute_force}         \tab TRUE permutes even where that is exceeded. \cr
#'   }
#' @param variable Character name of variable (column in samples) to permute.
#' @return A list (the permuted state) with the following elements:
#'   \tabular{ll}{
#'     \code{expression} \tab Numeric matrix with non-negative expression values. \cr
#'     \code{features}   \tab A data.frame with feature meta-data for rows of expression. \cr
#'     \code{samples}    \tab A data.frame with observation meta-data for columns of expression. \cr
#'   } 
#' @examples
#' set.seed(101)
#' exprs <- h0testr::sim1(n_obs=6, n_feats=12)$mat
#' feats <- data.frame(feature_id=rownames(exprs))
#' samps <- data.frame(observation_id=colnames(exprs), age=c(rep("young", 3), rep("old", 3)))
#' state <- list(expression=exprs, features=feats, samples=samps)
#' print(state)
#' config <- list(feat_col="feature_id", obs_col="observation_id",
#'   sample_id_col="observation_id", permute_var="age")
#' state2 <- h0testr::permute(state, config)
#' print(state2)
#' @export

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
    within <- if(is.null(config$permute_within)) "" else config$permute_within

    if(nzchar(within)) {
      if(!(within %in% colnames(state$samples))) {
        f.err(
          "permute: !(config$permute_within %in% colnames(state$samples))", "\n",
          "config$permute_within:", within, "\n",
          "colnames(state$samples):", colnames(state$samples),
          config=config
        )
      }
      if(any(is.na(state$samples[[within]]))) {
        f.err("permute: config$permute_within has", sum(is.na(state$samples[[within]])),
          "NA value(s), so those observations belong to no stratum; column:", within,
          config=config)
      }
      ## shuffle is over unique sample ids, so a stratum has to be a property of 
      ## sample and not of observation:
      n_str <- tapply(state$samples[[within]], state$samples[[config$sample_id_col]],
        function(v) length(unique(v)))
      if(any(n_str > 1)) {
        f.err("permute: config$permute_within is not constant within",
          paste0("config$sample_id_col; ", sum(n_str > 1)), "sample(s) span more than one",
          "stratum; permute_within:", within, "; sample_id_col:", config$sample_id_col,
          config=config)
      }
      if(f.n_informative(state$samples, variable, within) < 1) {
        f.err("permute: no stratum holds more than one level of", variable, "so permuting",
          "within", within, "returns the data unchanged and estimates no null;",
          variable, "is nested in", within, config=config)
      }
    }

    chk <- check_permutation(config, state$samples, variable=variable)
    if(!chk$ok) {
      fix <- if(length(chk$suggest)) {
        paste0("set config$permute_within to one of: ", paste(chk$suggest, collapse=", "))
      } else "no column of config$frm restores it; consider a different test"
      msg <- c("permuting", variable,
        if(nzchar(within)) paste("within", within) else "over all samples",
        "inflates the variance of the", config$test_term, "coefficient",
        paste0(signif(chk$inflation, 3), "x,"), "above config$permute_inflation_max",
        paste0(chk$threshold, ";"), "permuted results will carry too few hits and any fdr",
        "estimated from them will be too low;", fix,
        "\n  to permute anyway, set config$permute_force=TRUE")
      if(isTRUE(config$permute_force)) {
        f.msg("WARNING: permute:", f.cat_args(msg), config=config)
      } else f.err(msg, config=config)
    }

    f.log_block("permuting", variable,
      if(nzchar(within)) paste("within", within) else "over all samples", config=config)
    state$samples <- f.permute_col(state$samples, config, variable, within)
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
#'     \code{dir_in}          \tab Path to directory where input files located (character). \cr
#'     \code{feature_file_in} \tab Name of file with feature metadata; assumed in \code{dir_in}. \cr
#'     \code{sample_file_in}  \tab Name of file with observation metadata; assumed in \code{dir_in}. \cr
#'     \code{data_file_in}    \tab Name of file with signal data; assumed in \code{dir_in}. \cr
#'     \code{feat_id_col}     \tab Name of column in \code{feature_file_in} that corresponds to rows of \code{data_file_in}. \cr
#'     \code{gene_id_col}     \tab Name of column in \code{feature_file_in} with unique gene group or protein group labels. \cr
#'     \code{obs_id_col}      \tab Name of column in \code{sample_file_in} that corresponds to columns of \code{data_file_in}. \cr
#'     \code{sample_id_col}   \tab Name of column in \code{sample_file_in} with unique sample labels. \cr
#'   }
#' @return A list (the initial state) with the following elements:
#'   \tabular{ll}{
#'     \code{expression} \tab Numeric matrix with expression values. \cr
#'     \code{features}   \tab A data.frame with feature meta-data for rows of expression. \cr
#'     \code{samples}    \tab A data.frame with observation meta-data for columns of expression. \cr
#'   } 
#' @examples
#' config <- h0testr::new_config()     ## defaults
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
#' @export

load_data <- function(config) {
  
  report_config(config)
  
  state <- read_data(config)
  out <- init_state(state, config)
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
