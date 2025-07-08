#' Get configuration template with mix of defaults and example values.
#' @description
#'   Returns a configuration list filled with defaults and
#'     example values. Configuration meant to be customized then passed to
#'     other functions.
#' @details 
#'   You must customize \code{frm}, \code{test_term}, and \code{sample_factors}.
#'     You often need to customize \code{feat_id_col}, \code{gene_id_col}, 
#'     \code{sample_id_col}, and \code{obs_id_col}.  
#' @return list of configuration values
#' @examples
#' config <- f.new_config()
#' 
#' ## you must customize frm, test_term, and sample_factors:
#' config$frm <- ~ age + sex + age:sex
#' config$test_term <- "age:sex"
#' config$sample_factors <- list(age=c("young", "old"), sex=c("female", "male"))
#'
#' ## these may or may not need customization, depending on your file formats:
#' config$feat_id_col <- "peptide_id"
#' config$gene_id_col <- "gene_id"
#' config$sample_id_col <- "sample"
#' config$obs_id_col <- "observation"
#' 
#' str(config)
#' f.report_config(config)

f.new_config <- function() {

  config <- list(

    ## input/output paths:
    feature_file_in="features.tsv",      ## feature annotation .tsv; row features
    sample_file_in="samples.tsv",        ## sample annotation .tsv; row observations
    data_file_in="expression.tsv",       ## quantification matrix .tsv; rowname features, colname observations
    dir_in=".",                          ## data_file_in, feature_file_in, and sample_file_in found here
    dir_out=".",                         ## output directory
    
    ## formula for testing: actual formula can have '+' and ':'; not tested w/ e.g. '*' yet.
    frm=~age+gender+age:gender,          ## formula with variable of interest and covariates
    test_term="age:gender",              ## term (scalar character) in $frm on which test is to be performed
    permute_var="",                      ## name (scalar character) of variable to permute; "" for no permutation (normal execution)
    sample_factors=list(                 ## set levels of factor variables in $frm
      age=c("young", "old"),             ## by default, numeric treated as numeric; if levels set here, treated as factor
      gender=c("Male", "Female")         ## by default, character treated as factor with alphabetically ordered levels
    ),
    
    ## samps and feats column names (new cols are introduced by the code):
    feat_id_col="feature_id",            ## column (scalar character) in feature_file_in matching rownames of data_file_in
    gene_id_col="gene_id",               ## column (scalar character) in feature_file_in with gene group or protein group ids 
    obs_id_col="observation_id",         ## column (scalar character) in sample_file_in matching colnames of data_file_in
    sample_id_col="sample_id",           ## column (scalar character) in sample_file_in with sample ids; not unique if tech reps; same as obs_id_col if no tech reps
    n_samples_expr_col="n_samps_expr",   ## new col (scalar character) for feature metadata; n samples expressing feature
    median_raw_col="median_raw",         ## new col (scalar character) for feature metadata; median feature expression in expressing samples
    n_features_expr_col="n_feats_expr",  ## new col (scalar character) for sample metadata; n features expressed
    obs_col="",                          ## for internal use; leave ""; samps[, obs_col] == colnames(exprs) throughout script
    feat_col="",                         ## for internal use; leave ""; feats[, feat_col] == rownames(exprs) throughout script
    
    ## output file naming:
    log_file="",                         ## log file path (character); or "" for log to console                 
    feature_mid_out=".features",         ## midfix for output feature files
    sample_mid_out=".samples",           ## midfix for output samples file
    data_mid_out=".expression",          ## midfix for output expression files
    result_mid_out=".results",           ## midfix for output results file
    suffix_out=".tsv",                   ## suffix for output files 
    
    ## tunable options: defaults are usually ok, except:
    ##   for dia: usually works ok: RLE:unif_sample_lod:0.05 for norm_method:impute_method:impute_quantile
    ##   for dda: usually works ok: quantile:0.75:unif_sample_lod:0 for norm_method:norm_quantile:impute_method:impute_quantile
    norm_method="quantile",              ## in c("TMM", "TMMwsp", "RLE", "upperquartile", "quantile", "cpm", "vsn", "qquantile", "log2", "none")
    norm_quantile=0.75,                  ## for quantile normalization; 0.5 is median; 0.75 is upper quartile;
    n_samples_min=2,                     ## min samples/feature w/ feature expression > 0 to keep feature
    n_features_min=1000,                 ## min features/sample w/ expression > 0 to keep sample
    ## in: c("sample_lod", "unif_sample_lod", "unif_global_lod", "rnorm_feature", "glm_binom", "loess_logit", "glmnet", "rf", "none")
    impute_method="sample_lod",
    impute_quantile=0.01,                ## quantile for unif_ imputation methods
    impute_scale=1,                      ## for rnorm_feature, adjustment on sd of distribution [1: no change];
    impute_span=0.5,                     ## loess span for f.impute_loess_logit
    impute_n_pts=1e7,                    ## granularity of imputed values for f.impute_glm_binom and f.impute_loess_logit
    test_method="trend",                 ## in c("voom", "trend")
    ## run_order character vector with elements from {"normalize", "combine_reps", "filter", "impute"}:
    run_order=c("normalize", "combine_reps", "filter", "impute"),   ## determines order of workflow operations
    
    ## misc; 
    save_state=TRUE,                     ## whether to save output files; might set to FALSE for tuning/testing
    probs=c(0, 0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99, 1.0),
    width=110,
    verbose=T
  )
  
  return(config)
}

#' Check configuration
#' @description
#'   Check configuration list for recognizable names and proper value types.
#' @details
#'   Only checks parameters in the configuration. Does not complain about 
#'     missing settings.
#' @param config List with configuration values like those returned by \code{f.new_config()}.
#' @return Logical scalar \code{TRUE} if configuration ok. Otherwise throws error.
#' @examples
#' config <- list()
#' h0testr::f.check_config(config)
#'
#' config <- h0testr::f.new_config()
#' h0testr::f.check_config(config)
#'
#' config$impute_quantile <- "0.01"
#' h0testr::f.check_config(config)

f.check_config <- function(config) {

  if(length(config) %in% 0) {
    f.msg("empty config ok", config=config)
    return(TRUE)
  }
  
  scalar_character <- c("feature_file_in", "sample_file_in", "data_file_in", 
    "dir_in", "dir_out", "test_term", "permute_var", 
    "feat_id_col", "gene_id_col", "feat_col",
    "obs_id_col", "sample_id_col", "obs_col", "n_samples_expr_col", 
    "median_raw_col", "n_features_expr_col", "log_file", "feature_mid_out", 
    "sample_mid_out", "data_mid_out", "result_mid_out", "suffix_out", 
    "norm_method", "impute_method", "test_method")
    
  scalar_counts <- c("n_samples_min", "n_features_min", "impute_n_pts", "width")
  scalar_props <- c("norm_quantile", "impute_quantile", 
    "impute_scale", "impute_span")
  scalar_positive <- c("impute_span")
  scalar_logical <- c("save_state", "verbose")
  scalar_formula <- c("frm")
  vector_character <- c("run_order")
  vector_props <- c("probs")
  list_character <- c("sample_factors")
  
  ## check all param names in config are non-emtpy and recognized:
  
  noms <- names(config)
  all_noms <- c(
    scalar_character, scalar_counts, scalar_props, 
    scalar_positive, scalar_logical, scalar_formula, vector_props, 
    list_character
  )
  
  if(is.null(noms)) f.err("f.check_config: is.null(names(config))", config=config)
  for(idx in 1:length(config)) {
    if(length(noms[idx]) %in% 0) {
      f.err("f.check_config: empty parameter name in config at idx:", 
        idx, config=config)
    }
    if(!(noms[idx] %in% all_noms)) {
      f.err("f.check_config: unexpected parameter name:", noms[idx], 
        config=config)
    }
  }  
  
  ## check param values:
  
  for(nom in scalar_character) {
    if(nom %in% config) {
      if(!(is.character(config[[nom]]) && length(config[[nom]] == 1))) {
        f.err("f.check_config: param not scalar character; param:",  nom, 
          "; value:", config[[nom]], config=config)
      }
    }
  }
  
  for(nom in scalar_counts) {
    if(nom %in% config) {
      if(!(is.numeric(config[[nom]]) && length(config[[nom]] == 1))) {
        f.err("f.check_config: param not scalar count; param:",  nom, 
          "; value:", config[[nom]], config=config)
      }
      if(config[[nom]] < 0) {
        f.err("f.check_config: param not non-negative integer:", nom,
          "; value:", config[[nom]], config=config)
      }
      if(config[[nom]] != round(config[[nom]])) {
        f.err("f.check_config: param not an integer:", nom,
          "; value:", config[[nom]], config=config)
      }
    }
  }
  
  for(nom in scalar_props) {
    if(nom %in% config) {
      if(!is.numeric(config[[nom]]) && length(config[[nom]] == 1)) {
        f.err("f.check_config: param not scalar proportion; param:",  nom, 
          "; value:", config[[nom]], config=config)
      }
      if(config[[nom]] < 0 || config[[nom]] > 1) {
        f.err("f.check_config: param not between 0 and 1:", nom,
          "; value:", config[[nom]], config=config)
      }
    }
  }
  
  for(nom in scalar_positive) {
    if(nom %in% config) {
      if(!is.numeric(config[[nom]]) && length(config[[nom]] == 1)) {
        f.err("f.check_config: param not scalar positive numeric; param:",  nom, 
          "; value:", config[[nom]], config=config)
      }
      if(config[[nom]] < 0) {
        f.err("f.check_config: param not non-negative:", nom,
          "; value:", config[[nom]], config=config)
      }
    }
  }
  
  for(nom in scalar_logical) {
    if(nom %in% config) {
      if(!is.logical(config[[nom]]) && length(config[[nom]] == 1)) {
        f.err("f.check_config: param not scalar logical; param:",  nom, 
          "; value:", config[[nom]], config=config)
      }
    }
  }
  
  for(nom in scalar_formula) {
    if(nom %in% config) {
      if(!methods::is(config[[nom]], "formula") && length(config[[nom]] == 1)) {
        f.err("f.check_config: param not scalar formula; param:",  nom, 
          "; value:", config[[nom]], config=config)
      }
    }
  }
  
  for(nom in vector_character) {
    if(nom %in% config) {
      if(!is.character(config[[nom]])) {
        f.err("f.check_config: param not vector of character; param:",  nom, 
          "; value:", config[[nom]], config=config)
      }
    }
  }
  
  for(nom in vector_props) {
    if(nom %in% config) {
      if(!is.numeric(config[[nom]])) {
        f.err("f.check_config: param not vector of proportions; param:",  nom, 
          "; value:", config[[nom]], config=config)
      }
      if(any(config[[nom]] < 0) || any(config[[nom]] > 1)) {
        f.err("f.check_config: proportions out of range [0, 1]; param:", nom,
          "; value:", config[[nom]], config=config)
      }
    }
  }
  
  for(nom in list_character) {
    if(nom %in% config) {
      if(!is.numeric(config[[nom]])) {
        f.err("f.check_config: param not list of character; param:",  nom, 
          "; value:", config[[nom]], config=config)
      }
    }
  }

  return(TRUE)
}

#' Report configuration
#' @description
#' Reports configuration settings used for run, and checks if
#'   \code{config$test_term} is compatible with \code{config$frm}.
#' @details Report written to \code{config$log_file}; if \code{config$log_file == ""},
#'   written to standard out (console or terminal). Only supports non-lists
#'     and lists of non-lists (not lists of lists) as \code{config} values.
#'   Throws error if \code{config$test_term} is not compatible with \code{config$frm}.
#' @param config List with configuration values.
#' @return NULL
#' @examples
#' config <- f.new_config()
#' 
#' ## you must customize frm, test_term, and sample_factors:
#' config$frm <- ~ age + sex + age:sex
#' config$test_term <- "age:sex"
#' config$sample_factors <- list(age=c("young", "old"), sex=c("female", "male"))
#'
#' ## these may or may not need customization, depending on your file formats:
#' config$feat_id_col <- "peptide_id"
#' config$gene_id_col <- "gene_id"
#' config$sample_id_col <- "sample"
#' config$obs_id_col <- "observation"
#' 
#' f.report_config(config)

f.report_config <- function(config) {

  f.check_config(config)

  for(k1 in names(config)) {
    v1 <- config[[k1]]
    if(is.list(v1)) {
      for(k2 in names(v1)) {
        f.msg(k1, ":", k2, ":", paste(as.character(v1[[k2]]), sep=", "), config=config)
      }
    } else f.msg(k1, ":", paste(as.character(v1), sep=", "), config=config)
  }
  
  ## check config$test_term compatible with config$frm; throws error if not, 
  ##   else returns NULL:

  trms <- as.character(config$frm)[2]
  trms <- gsub("[[:space:]]+", "", trms)
  trms <- unlist(strsplit(trms, split="\\+"))
  
  if(length(trms) %in% 0) {
    f.err("f.report_config: length(trms) %in% 0", config=config)
  }  
  
  if(!(config$test_term %in% trms)) {
    f.err("f.report_config: config$test_term", config$test_term, "not in config$frm", 
      as.character(config$frm), config=config)
  }
}
