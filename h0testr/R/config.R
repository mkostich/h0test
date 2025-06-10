#' Get configuration template with mix of defaults and example values.
#' @description
#' `f.new_config` returns a configuration list filled with defaults and
#'   example values. Configuration meant to be customized then passed to
#'   other functions.
#' @details Report written to config$log_file; if `config$log_file == ""`,
#'   written to standard out (console or terminal). Only supports non-lists
#'     and lists of non-lists (not lists of lists) as `config` values.
#'   Throws error if `config$test_term` is not compatible with `config$frm`.
#' @return list of configuration values
#' @examples
#'   \dontrun{
#'     config <- f.new_config()
#'     config$frm <- ~ age + sex + age:sex
#'     config$test_term <- "age:sex"
#'     config$sample_factors <- list(age=c("young", "old"), sex=c("female", "male"))
#'     config$feat_id_col <- "gene"
#'     config$sample_col <- "sample"
#'     config$obs_col <- "observation"
#'     f.report_config(config)
#'   }

f.new_config <- function() {

  config <- list(
    ## input/output paths:
    feature_file_in="features.tsv",      ## feature annotation .tsv; row features
    sample_file_in="samples.tsv",        ## sample annotation .tsv; row observations
    data_file_in="expression.tsv",       ## lfq normalized quant matrix .tsv; row features, column observations
    dir_in=".",                          ## where DATA_FILE_IN, FETURE_FILE_IN, and SAMPLE_FILE_IN live
    dir_out=".",                         ## output directory
    ## formula for testing: actual formula can have '+' and ':'; not tested w/ e.g. '*' yet.
    frm=~age+strain+gender+age:strain,   ## formula with variable of interest and covariates
    test_term="age:strain",              ## term in FRM on which test is to be performed
    permute_var="",                      ## variable to permute; "" for no permutation (normal execution)
    sample_factors=list(                 ## set levels of factor variables in sample metadata
      age=c("4", "12", "24"),            ## numeric treated as numeric unless levels set here, then as factor
      strain=c("B6", "A_J", "Balbc_J"),
      gender=c("Male", "Female")
    ),
    ## samps and feats column names (new cols are introduced by the code):
    feat_id_col="gene_id",               ## feats[, FEAT_ID_COL] == rownames(exprs)
    obs_col="assay_id",                  ## unique id for observations; samps[, OBS_COL] == colnames(exprs)
    sample_col="sample_id",              ## id for samples in samps; not unique if tech reps;
    n_samples_expr_col="n_samps_expr",   ## new col; n samples expressing feature
    median_raw_col="median_raw",         ## new col; median feature expression in expressing samples
    n_features_expr_col="n_feats_expr",  ## new col; n features express
    ## output file naming:
    log_file="log.txt",                  ## log file path; or "" for log to console                 
    feature_mid_out=".features",         ## midfix for output feature files
    sample_mid_out=".samples",           ## midfix for output samples file
    data_mid_out=".expression",          ## midfix for output expression files
    result_mid_out=".results",           ## prefix for output results file
    suffix_out=".tsv",                   ## suffix for output files  ## tunable options: defaults are usually ok, except:
    ##   for dia: usually works ok: RLE:unif_sample_lod:0.05 for NORM_METHOD:IMPUTE_METHOD:IMPUTE_QUANTILE
    ##   for dda: usually works ok: quantile:0.75:unif_sample_lod:0 for NORM_METHOD:NORM_QUANTILE:IMPUTE_METHOD:IMPUTE_QUANTILE
    norm_method="quantile",              ## in c("vsn","cpm","quantile","qquantile","TMM","TMMwsp","RLE","upperquartile")
    norm_quantile=0.75,                  ## for quantile normalization; 0.5 is median; 0.75 is upper quartile;
    transform_method="log2",             ## in c('log2', 'log10', 'none')
    n_samples_min=2,                     ## min samples/feature w/ feature expression > 0 to keep feature
    n_features_min=1000,                 ## min features/sample w/ expression > 0 to keep sample
    impute_method="unif_sample_lod",     ## in c("unif_global_lod", "unif_sample_lod", "sample_lod", "rnorm_feature")
    impute_quantile=0.01,                ## quantile for unif_ imputation methods
    impute_scale=1,                      ## for rnorm_feature, adjustment on sd of distribution [1: no change];
    impute_granularity=0.0001,           ## granularity of imputed values for f.impute_glm_binom and f.impute_loess_logit
    impute_span=0.25,                    ## loess span for f.impute_loess_logit 
    test_method="trend",                 ## in c("voom", "trend")
    ## misc:
    probs=c(0, 0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99, 1.0),
    width=110,
    verbose=T
  )
  
  return(config)
}


