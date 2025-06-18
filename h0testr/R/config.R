#' Get configuration template with mix of defaults and example values.
#' @description
#'   Returns a configuration list filled with defaults and
#'     example values. Configuration meant to be customized then passed to
#'     other functions.
#' @details Report written to \code{config$log_file}; if \code{config$log_file == ""},
#'   written to standard out (console or terminal). Only supports non-lists
#'     and lists of non-lists (not lists of lists) as \code{config} values.
#'   Throws error if \code{config$test_term} is not compatible with \code{config$frm}.
#' @return list of configuration values
#' @examples
#'   \dontrun{
#'     config <- f.new_config()
#'     config$frm <- ~ age + sex + age:sex
#'     config$test_term <- "age:sex"
#'     config$sample_factors <- list(age=c("young", "old"), sex=c("female", "male"))
#'     config$feat_id_col <- "gene"
#'     config$sample_id_col <- "sample"
#'     config$obs_id_col <- "observation"
#'     config$log_file <- "/path/to/log.txt"
#'     f.report_config(config)
#'   }

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
    test_term="age:gender",              ## term (character) in $frm on which test is to be performed
    permute_var="",                      ## name (character) of variable to permute; "" for no permutation (normal execution)
    sample_factors=list(                 ## set levels of factor variables in $frm
      age=c("young", "old"),             ## by default, numeric treated as numeric; if levels set here, treated as factor
      gender=c("Male", "Female")         ## by default, character treated as factor with alphabetically ordered levels
    ),
    
    ## samps and feats column names (new cols are introduced by the code):
    feat_id_col="gene_id",               ## column (character) in feature_file_in matching rownames of data_file_in
    obs_id_col="observation_id",         ## column in sample_file_in matching colnames of data_file_in
    sample_id_col="sample_id",           ## column in sample_file_in with sample ids; not unique if tech reps; same as obs_id_col if no tech reps
    obs_col="",                          ## for internal use; leave ""; samps[, obs_col] == colnames(exprs) throughout script
    n_samples_expr_col="n_samps_expr",   ## new col for feature metadata; n samples expressing feature
    median_raw_col="median_raw",         ## new col for feature metadata; median feature expression in expressing samples
    n_features_expr_col="n_feats_expr",  ## new col for sample metadata; n features expressed
    
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
    impute_granularity=0.0001,           ## granularity of imputed values for f.impute_glm_binom and f.impute_loess_logit
    impute_span=0.25,                    ## loess span for f.impute_loess_logit 
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
