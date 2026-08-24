#' Get configuration template with default and example values
#' @description
#'   Returns a configuration list filled with defaults and
#'     example values. Configuration meant to be customized then passed to
#'     other functions.
#' @details 
#'   This function is meant to simplify generation of configurations for 
#'     higher level functions, like \code{h0testr::run()}, or 
#'     \code{h0testr::tune()}. For most other functions, you can pass a 
#'     simpler config as a list containing only the needed parameters. See 
#'     documentation and examples for the function of interest for the minimal 
#'     configuration needed.
#'   For hypothesis testing or calls to \code{h0testr::initialize()},
#'     customize \code{frm}, \code{test_term}, and \code{reference_levels}.
#'   \code{contrast} is the alternative to \code{test_term}: it names a weighted
#'     sum of the coefficients of \code{frm} to test within the full model, written
#'     as arithmetic over coefficient names, e.g. \code{"grpb - grpc"} or
#'     \code{"(grpb + grpc)/2 - grpd"}. A coefficient whose name is not a syntactic
#'     name, such as an interaction, has to be backquoted, e.g.
#'     \code{"`sexM:batchb2`"}. Constants may scale a coefficient; two coefficient
#'     names may not be multiplied together, that not being a weighted sum.
#'   The two are mutually exclusive and one run tests one hypothesis, so set
#'     \code{test_term} to \code{""} to test a contrast. The difference is what
#'     each answers: \code{test_term} obeys marginality, testing the named term
#'     together with every term containing it, which for a factor of more than two
#'     levels is a joint test over all of its coefficients; a contrast is one
#'     degree of freedom, so it can compare two particular levels, and reaches
#'     \code{test_method="deqms"}, which cannot run a joint test at all.
#'   When using the config to load data from files (e.g. by calling 
#'     \code{h0testr::load_data(config)}), calling 
#'     \code{h0testr::initialize()}, or for aggregating multiple 
#'     observations per sample (e.g. by calling \code{h0testr::combine_replicates()}),
#'     or for aggregating precursors into gene/protein groups (e.g. by calling 
#'     \code{h0testr::combine_precursors()}), you should customize 
#'     \code{feat_id_col}, \code{gene_id_col}, \code{sample_id_col}, 
#'     and \code{obs_id_col}. In these cases, leave \code{feat_col} and
#'     \code{obs_col} as \code{""} (they will be automatically set and 
#'     changed after observation or precursor aggregation.
#'   Otherwise, you may only need to set \code{feat_col} and \code{obs_col}.
#'     see examples for the function of interest to see what is needed.
#'   \code{is_log_transformed} declares the scale of \code{state$expression} and
#'     is the single place the rest of the package looks to find it out. The
#'     default \code{FALSE} says the input is raw, so \code{h0testr::initialize()}
#'     converts zeros to \code{NA} and rejects negative values, and
#'     \code{h0testr::normalize()} may transform it. Set it to \code{TRUE} for
#'     input that has already been log transformed or otherwise put on a log-like
#'     scale: zeros and negative values are then left alone,
#'     \code{h0testr::normalize()} accepts no method other than \code{"none"},
#'     and the imputation and testing functions know not to transform it again.
#'     \code{h0testr::normalize()} sets it to \code{TRUE} once it has transformed
#'     the data, so a workflow does not have to track the scale itself.
#'   \code{impute_floor_offset} is an offset, in log2 units, from the global
#'     minimum observed value, and it is used only when
#'     \code{is_log_transformed} is \code{TRUE}. The \code{unif_} imputation
#'     methods draw from an interval whose lower bound is
#'     \code{min(state$expression, na.rm=TRUE) + impute_floor_offset}. The
#'     default \code{-1} puts that bound one log2 unit below the dimmest value
#'     actually measured, so the interval has width even when
#'     \code{impute_quantile} is \code{0}. A less negative value narrows the
#'     interval, which reduces the variability of the imputed values and makes
#'     the test anti-conservative; \code{0} leaves no width at all when
#'     \code{impute_quantile} is \code{0}. A more negative value imputes dimmer
#'     values, spread more widely. It is ignored for raw input, where the
#'     imputation floor is zero abundance.
#' @return list of configuration values
#' @examples
#' config <- new_config()    ## all possible settings with defaults
#' str(config)                 ## check out the default settings
#' 
#' ## you MUST customize frm, test_term, and reference_levels:
#' config$frm <- ~ age + sex + age:sex
#' config$test_term <- "age:sex"
#' config$reference_levels <- c(age="young", sex="female")
#'
#' ## these may or may not need customization, depending on your file formats:
#' config$feat_id_col <- "peptide_id"
#' config$gene_id_col <- "gene_id"
#' config$sample_id_col <- "sample"
#' config$obs_id_col <- "observation"
#' 
#' str(config)
#' report_config(config)

new_config <- function() {

  config <- list(

    ## input/output paths:
    feature_file_in="features.tsv",      ## feature annotation .tsv; row features
    sample_file_in="samples.tsv",        ## sample annotation .tsv; row observations
    data_file_in="expression.tsv",       ## quantification matrix .tsv; rowname features, colname observations
    dir_in=".",                          ## data_file_in, feature_file_in, and sample_file_in found here
    dir_out=".",                         ## output directory
    
    ## match up expression dimnames with features metadata and samples metadata:
    feat_id_col="feature_id",            ## column (scalar character) in feature_file_in matching rownames of data_file_in
    gene_id_col="gene_id",               ## column (scalar character) in feature_file_in with gene group or protein group ids 
    obs_id_col="observation_id",         ## column (scalar character) in sample_file_in matching colnames of data_file_in
    sample_id_col="sample_id",           ## column (scalar character) in sample_file_in with sample ids; not unique if tech reps; same as obs_id_col if no tech reps
    obs_col="",                          ## for internal use; leave ""; samps[, obs_col] == colnames(exprs) throughout script
    feat_col="",                         ## for internal use; leave ""; feats[, feat_col] == rownames(exprs) throughout script
    
    ## formula for testing: actual formula can have '+' and ':'; not tested w/ e.g. '*' yet.
    frm=~age+gender+age:gender,          ## formula with variable of interest and covariates
    test_term="age:gender",              ## term (scalar character) in $frm on which test is to be performed; "" iff $contrast is set
    contrast="",                         ## weighted sum (scalar character) of coefficients of $frm to test, e.g. "grpb - grpc"; "" for none (test $test_term instead)
    permute_var="",                      ## name (scalar character) of variable to permute; "" for no permutation (normal execution)
    reference_levels=c(                  ## reference level of each factor variable in $frm
      age="young",                       ## numeric variable treated as continuous unless named here
      gender="Male"                      ## character variable must be named here; else error
    ),
    n_distinct_numeric_warn=5,           ## warn if continuous variable in $frm has this few distinct values

    ## new cols introduced into metadata data.frames by the code:
    n_samples_expr_col="n_samps_expr",   ## new col (scalar character) for feature metadata; n samples expressing feature
    median_raw_col="median_raw",         ## new col (scalar character) for feature metadata; median feature expression in expressing samples
    n_features_expr_col="n_feats_expr",  ## new col (scalar character) for sample metadata; n features expressed
    n_feats_col="n_feats",               ## new col (scalar character) for gene-level feature metadata; n features aggregated into each gene; preserved if already present
    combine_method_col="combine_method", ## new col (scalar character) for gene-level feature metadata; which aggregator summarized each gene; preserved if already present
    df_test_col="df_test",               ## new col (scalar character) for feature metadata; estimable df for config$test_term
    df_resid_col="df_resid",             ## new col (scalar character) for feature metadata; residual df of model fitted to feature

    ## output file naming:
    log_file="",                         ## log file path (character); or "" for log to console                 
    feature_mid_out=".features",         ## midfix for output feature files
    sample_mid_out=".samples",           ## midfix for output samples file
    data_mid_out=".expression",          ## midfix for output expression files
    result_mid_out=".results",           ## midfix for output results file
    suffix_out=".tsv",                   ## suffix for output files 
    
    ## tunable options: defaults are usually ok, except:
    ##   for dia: usually works ok: RLE:unif_sample_lod:0.05 for normalization_method:impute_method:impute_quantile
    ##   for dda: usually works ok: quantile:0.75:unif_sample_lod:0 for normalization_method:normalization_quantile:impute_method:impute_quantile
    is_log_transformed=FALSE,            ## whether state$expression is already on a log-like scale; FALSE means raw, so initialize() converts zeros to NA and rejects negative values; set TRUE by normalize()
    normalization_method="RLE",          ## normalization method; h0testr::normalize_methods() retuns options.
    normalization_quantile=0.75,         ## for quantile normalization; 0.5 is median; 0.75 is upper quartile;
    normalization_span=0.7,              ## span for normalize_loess()
    n_samples_min=2,                     ## min samples/feature w/ non-NA feature expression to keep feature
    n_features_min=1000,                 ## min features/sample w/ non-NA expression to keep sample
    estimability="test",                 ## estimability required of config$test_term; in c("test", "term", "full")
    df_resid_min=2,                      ## min residual degrees of freedom per feature to keep feature
    feature_aggregation="medianPolish",  ## in c("medianPolish", "robustSummary", "none")
    impute_method="sample_lod",          ## method for imputing missing values; h0testr::impute_methods() returns options.
    impute_quantile=0.01,                ## quantile for unif_ imputation methods
    impute_floor_offset=-1,              ## offset (non-positive; log2 units) from the global observed minimum giving the lower bound of the unif_ imputation interval, when config$is_log_transformed is TRUE
    impute_scale=1,                      ## for rnorm_feature, adjustment on sd of distribution [1: no change];
    impute_span=0.5,                     ## loess span for impute_loess_logit()
    impute_k=7,                          ## k for impute_knn() or impute_lls()
    impute_npcs=5,                       ## number of PCs for PCA-based imputations "bpca", "ppca", and "svdImpute"
    impute_alpha=1,                      ## alpha mixing parameter for impute_glmnet()
    impute_n_pts=1e7,                    ## granularity of imputed values for impute_glm_binom() and impute_loess_logit()
    impute_aug_steps=3,                  ## data augmentation iterations for impute_rf() and impute_glmnet()
    test_method="trend",                 ## hypothesis test method; h0testr::test_methods() returns options.
    test_prior_df=3,                     ## prior df for test_proda()
    test_moderate=TRUE,                  ## whether to shrink the per-feature error variance across features before testing; used by test_prolfqua(), which is the only method here that can be told not to; the limma-based methods always moderate and proda always shrinks
    test_trend=FALSE,                    ## whether the prior variance of that shrinkage is fitted against mean feature intensity instead of being flat; honored by test_method "prolfqua", where it makes the prior the one test_method="trend" uses, and by "deqms", where it sets the limma prior beneath DEqMS's count-based one; "trend" always trends and the remaining methods cannot, which test() warns about rather than refusing; overridden by test(trend=); unrelated to test_method
    test_random_obs=TRUE,                ## whether the feature level mixed model paths, test_method %in% c("prolfqua_lmer", "msqrob_agg"), add a random observation effect to the random feature effect; TRUE is the calibrated model; FALSE gives the feature-only structure both packages document, which is anti-conservative; see test_prolfqua() and test_msqrob()
    test_ridge=FALSE,                    ## whether test_method="msqrob_agg" penalizes the fixed effects, msqrob2::msqrobAggregate(ridge=TRUE); FALSE is msqrob2's own default throughout; TRUE shrinks the coefficients toward zero, so the reported logFC is not comparable to what the other methods report, renames the fitted parameters, and refuses a mean model with fewer than two non-intercept columns; see test_msqrob()

    ## run_order character vector with elements from f.run_order_steps(), which is {"normalize", "combine_replicates", "combine_features", "filter", "impute"}; check_config() refuses any other name, and run() warns about a repeat:
    run_order=c("normalize", "combine_replicates", "combine_features", "filter", "impute"),   ## order of workflow operations
    
    ## misc; 
    save_state=TRUE,                     ## whether to save output files; might set to FALSE for tuning/testing
    probs=c(0, 0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99, 1.0),
    width=110,
    verbose=T
  )
  
  return(config)
}

## helper for check_config(): a lower bound on the number of non-intercept columns
##   stats::model.matrix() will build from config$frm, or NA when the config does not
##   settle it. Every term contributes at least one column, so a formula with two or
##   more terms already has at least two and the count need not be exact. A single term
##   is counted only when it names one variable whose type initialize() has resolved: a
##   numeric contributes one column, and a factor one per level, less one for a design
##   that keeps its intercept. Anything else, an interaction or a variable not yet
##   classified, is NA, since how many columns it becomes depends on the data:

f.frm_min_noint_cols <- function(config) {

  if(!f.is_formula(config$frm)) return(NA_integer_)

  trms <- try(stats::terms(config$frm), silent=T)
  if(inherits(trms, "try-error")) return(NA_integer_)

  labs <- attr(trms, "term.labels")           ## excludes any dependent variable
  if(length(labs) %in% 0) return(0L)
  if(length(labs) > 1) return(length(labs))

  nom <- labs[1]
  if(!(nom %in% names(config$covariate_types))) return(NA_integer_)

  typ <- config$covariate_types[[nom]]
  if(is.na(typ)) return(NA_integer_)
  if(typ %in% "numeric") return(1L)
  if(!(nom %in% names(config$factor_levels))) return(NA_integer_)

  n_lev <- length(config$factor_levels[[nom]])
  if(n_lev %in% 0) return(NA_integer_)

  return(as.integer(n_lev - (attr(trms, "intercept") %in% 1)))
}

#' Check configuration
#' @description
#'   Check configuration list for recognizable names and proper value types.
#' @details
#'   Only checks parameters in the configuration. Does not complain about
#'     missing settings.
#'   A numeric parameter has to be a single finite value. \code{NA} and
#'     \code{NaN} are refused by name, having reached a comparison as R's own
#'     "missing value where TRUE/FALSE needed" before, and so is an infinite
#'     count, which used to pass every check since \code{Inf == round(Inf)}. An
#'     infinite proportion is refused as out of range rather than as not finite,
#'     which says more about a proportion.
#'   Beyond the type of each value, the combinations that no other setting can rescue
#'     are refused here rather than at the step that would meet them:
#'     \code{config$test_method} outside \code{h0testr::test_methods()},
#'     \code{config$normalization_method} outside \code{h0testr::normalize_methods()},
#'     \code{config$impute_method} outside \code{h0testr::impute_methods()}, and
#'     \code{config$feature_aggregation} outside \code{"medianPolish"},
#'     \code{"robustSummary"} and \code{"none"}, each of which would otherwise be found
#'     only when that step ran, which for \code{config$impute_method} is after the whole
#'     rest of the workflow; \code{""} is allowed for all four and means unset, the step's
#'     own \code{method} argument naming the method instead;
#'     \code{config$contrast} and \code{config$test_term} both naming something to
#'     test, one run testing one hypothesis; \code{config$test_ridge=TRUE} on a
#'     mean model that cannot carry a penalty, which needs at least two non-intercept
#'     columns; and a \code{config$run_order} naming a step that is not one of the
#'     workflow steps \code{h0testr::run()} can walk, which would otherwise fail as an
#'     "object not found" partway through the workflow, with the steps before it already
#'     done. A step named twice is not refused, since running one again may be meant;
#'     \code{h0testr::run()} warns about it once, its output files being named after the
#'     first occurrence. The limits that depend on the data rather than on the configuration,
#'     such as how many features a gene has, stay with the engine that meets them; see
#'     \code{h0testr::test_msqrob()} and \code{h0testr::test_deqms()}.
#'   A setting that the chosen method does not consult is not refused and is not
#'     reported here: it is a \code{NOTE} from \code{h0testr::test()}, which every
#'     workflow calls once, this function being called once per step.
#'   See documentation for \code{h0testr::new_config()}
#'     for more detailed description of configuration parameters.
#' @param config List with configuration values like those returned by \code{new_config()}.
#' @return Logical scalar \code{TRUE} if configuration ok. Otherwise throws error.
#' @examples
#' config <- list()
#' h0testr::check_config(config)
#'
#' config <- h0testr::new_config()
#' h0testr::check_config(config)
#'
#' ## invalid: value must be numeric, not character; throws an error:
#' config$impute_quantile <- "0.01"
#' try(h0testr::check_config(config))
#'
#' ## invalid: not one of h0testr::test_methods(); throws an error:
#' config <- h0testr::new_config()
#' config$test_method <- "trrend"
#' try(h0testr::check_config(config))
#'
#' ## invalid: not one of h0testr::impute_methods(); throws an error here rather than
#' ##   after the rest of the workflow has already run:
#' config <- h0testr::new_config()
#' config$impute_method <- "unif_sample_lodd"
#' try(h0testr::check_config(config))
#'
#' ## invalid: a config written against a version that had the key; says what happened
#' ##   to it rather than only naming it:
#' config <- h0testr::new_config()
#' config$feature_aggregation_scaled <- TRUE
#' try(h0testr::check_config(config))
#'
#' ## invalid: the penalty needs at least two non-intercept columns, and a two level
#' ##   factor with an intercept gives one; throws an error once initialize() has
#' ##   resolved the levels:
#' config <- h0testr::new_config()
#' config$test_method <- "msqrob_agg"
#' config$test_ridge <- TRUE
#' config$frm <- ~grp
#' config$covariate_types <- c(grp="factor")
#' config$factor_levels <- list(grp=c("ctl", "trt"))
#' try(h0testr::check_config(config))
#'
#' ## ok: the same penalty on a design without an intercept, two columns:
#' config$frm <- ~0 + grp
#' h0testr::check_config(config)
#'
#' ## ok: a setting the chosen method does not consult is not this function's business;
#' ##   h0testr::test() notes it instead:
#' config <- h0testr::new_config()
#' config$test_method <- "trend"
#' config$test_ridge <- TRUE
#' h0testr::check_config(config)
#'
#' ## invalid: a step h0testr::run() has no function for; throws an error before the
#' ##   workflow starts rather than after the steps before it have run:
#' config <- h0testr::new_config()
#' config$run_order <- c("normalize", "combine_featurez", "filter")
#' try(h0testr::check_config(config))
#'
#' ## ok: no steps at all, meaning load_data() and then test(), for data that arrives
#' ##   already normalized, aggregated and complete:
#' config$run_order <- character(0)
#' h0testr::check_config(config)

check_config <- function(config) {

  if(length(config) %in% 0) {
    f.msg("empty config ok", config=config)
    return(TRUE)
  }
  
  scalar_character <- c("feature_file_in", "sample_file_in", "data_file_in", 
    "dir_in", "dir_out", "test_term", "contrast", "permute_var",
    "feat_id_col", "gene_id_col", "feat_col",
    "obs_id_col", "sample_id_col", "obs_col", "n_samples_expr_col",
    "median_raw_col", "n_features_expr_col", "n_feats_col", "combine_method_col",
    "df_test_col", "df_resid_col",
    "log_file", "feature_mid_out",
    "sample_mid_out", "data_mid_out", "result_mid_out", "suffix_out",
    "normalization_method", "feature_aggregation", "impute_method", "test_method",
    "estimability")

  scalar_counts <- c("n_samples_min", "n_features_min", "impute_n_pts",
    "impute_k", "impute_npcs", "impute_aug_steps", "test_prior_df",
    "n_distinct_numeric_warn", "width", "df_resid_min")
  
  scalar_props <- c("normalization_quantile", "impute_quantile", "impute_span",
    "impute_alpha", "normalization_span")
  scalar_positive <- c("impute_scale")
  scalar_nonpositive <- c("impute_floor_offset")
  ## log_from_raw is set by normalize(), not by the user; see f.check_state():
  scalar_logical <- c("save_state", "verbose",
    "is_log_transformed", "log_from_raw", "test_moderate", "test_trend",
    "test_random_obs", "test_ridge")
  scalar_formula <- c("frm")
  ## covariate_types is set by initialize(), not by the user; see
  ##   f.covariate_types():
  vector_character <- c("run_order", "covariate_types", "reference_levels")
  vector_props <- c("probs")
  ## factor_levels is set by initialize(), not by the user; see
  ##   f.set_covariate_factor_levels():
  list_character <- c("factor_levels")
  
  ## check all param names in config are non-emtpy and recognized:
  
  noms <- names(config)
  all_noms <- c(
    scalar_character, scalar_counts, scalar_props,
    scalar_positive, scalar_nonpositive, scalar_logical, scalar_formula,
    vector_character, vector_props,
    list_character
  )

  ## retired parameters, checked before the unrecognized name loop below so that
  ##   a config written against an older version says what to do instead of just
  ##   naming the offending parameter:

  retired <- c(
    zeros_to_na=paste("use is_log_transformed instead, with the opposite sense:",
      "is_log_transformed=FALSE means the input is raw, so zeros become NA and",
      "negative values are an error"),
    feature_aggregation_scaled=paste("rescaling before aggregation is no longer",
      "supported, having only ever been refused: it divided each feature by its",
      "own mean, which is a raw scale operation, but aggregation requires log",
      "scale data; drop the key. See h0testr::combine_features()")
  )
  for(nom in names(retired)) {
    if(nom %in% noms) {
      f.err("check_config: retired parameter:", nom, ";", retired[[nom]],
        config=config)
    }
  }
  
  if(is.null(noms)) f.err("check_config: is.null(names(config))", config=config)
  for(idx in 1:length(config)) {
    if(nchar(noms[idx]) %in% 0) {
      f.err("check_config: empty parameter name in config at idx:", 
        idx, config=config)
    }
    if(!(noms[idx] %in% all_noms)) {
      f.err("check_config: unexpected parameter name:", noms[idx], 
        config=config)
    }
  }  
  
  ## check param values:
  
  for(nom in scalar_character) {
    if(nom %in% names(config)) {
      if(!(is.character(config[[nom]]) && length(config[[nom]]) == 1)) {
        f.err("check_config: param not scalar character; param:",  nom, 
          "; value:", config[[nom]], config=config)
      }
    }
  }
  
  for(nom in scalar_counts) {
    if(nom %in% names(config)) {
      if(!(is.numeric(config[[nom]]) && length(config[[nom]]) == 1)) {
        f.err("check_config: param not scalar count; param:",  nom, 
          "; value:", config[[nom]], config=config)
      }

      ## is.numeric() is TRUE for NA_real_ and for NaN, either of which then
      ##   reached the comparison below as R's own "missing value where
      ##   TRUE/FALSE needed", naming neither the parameter nor this function
      ##   (measured). A logical NA is caught by the is.numeric() test above:

      if(is.na(config[[nom]])) {
        f.err("check_config: param is NA or NaN; param:", nom,
          "; value:", config[[nom]], config=config)
      }

      ## Inf passed every check below it: Inf == round(Inf), and Inf is not
      ##   < 0, so an infinite count was accepted outright (measured on
      ##   n_features_min and impute_k). Same test and wording as the
      ##   scalar_nonpositive loop further down. -Inf is reported as not
      ##   finite now rather than as not non-negative, which it also is:

      if(!is.finite(config[[nom]])) {
        f.err("check_config: param not finite:", nom,
          "; value:", config[[nom]], config=config)
      }
      if(config[[nom]] < 0) {
        f.err("check_config: param not non-negative integer:", nom,
          "; value:", config[[nom]], config=config)
      }
      if(config[[nom]] != round(config[[nom]])) {
        f.err("check_config: param not an integer:", nom,
          "; value:", config[[nom]], config=config)
      }
    }
  }
  
  for(nom in scalar_props) {
    if(nom %in% names(config)) {
      if(!(is.numeric(config[[nom]]) && length(config[[nom]]) == 1)) {
        f.err("check_config: param not scalar proportion; param:",  nom, 
          "; value:", config[[nom]], config=config)
      }

      ## as in the scalar_counts loop above: NA_real_ and NaN are numeric, and
      ##   the range check below is where R, rather than this function,
      ##   reported them:

      if(is.na(config[[nom]])) {
        f.err("check_config: param is NA or NaN; param:", nom,
          "; value:", config[[nom]], config=config)
      }
      if(config[[nom]] < 0 || config[[nom]] > 1) {
        f.err("check_config: param not between 0 and 1:", nom,
          "; value:", config[[nom]], config=config)
      }
    }
  }
  
  for(nom in scalar_positive) {
    if(nom %in% names(config)) {
      if(!(is.numeric(config[[nom]]) && length(config[[nom]]) == 1)) {
        f.err("check_config: param not scalar positive numeric; param:",  nom, 
          "; value:", config[[nom]], config=config)
      }

      ## as in the two loops above: is.numeric() is TRUE for NA_real_ and for
      ##   NaN, which then reached the comparison below as R's own "missing
      ##   value where TRUE/FALSE needed", and Inf was accepted outright
      ##   (both measured on impute_scale):

      if(is.na(config[[nom]])) {
        f.err("check_config: param is NA or NaN; param:", nom,
          "; value:", config[[nom]], config=config)
      }

      if(!is.finite(config[[nom]])) {
        f.err("check_config: param not finite:", nom,
          "; value:", config[[nom]], config=config)
      }

      if(config[[nom]] < 0) {
        f.err("check_config: param not non-negative:", nom,
          "; value:", config[[nom]], config=config)
      }
    }
  }

  ## an offset from an observed value, so it has to be finite; a positive offset
  ##   would put the bound it derives above the value it is measured from:

  for(nom in scalar_nonpositive) {
    if(nom %in% names(config)) {
      if(!(is.numeric(config[[nom]]) && length(config[[nom]]) == 1)) {
        f.err("check_config: param not scalar non-positive numeric; param:", nom,
          "; value:", config[[nom]], config=config)
      }
      if(!is.finite(config[[nom]])) {
        f.err("check_config: param not finite:", nom,
          "; value:", config[[nom]], config=config)
      }
      if(config[[nom]] > 0) {
        f.err("check_config: param not non-positive:", nom,
          "; value:", config[[nom]], config=config)
      }
    }
  }

  for(nom in scalar_logical) {
    if(nom %in% names(config)) {
      if(!(is.logical(config[[nom]]) && length(config[[nom]]) == 1)) {
        f.err("check_config: param not scalar logical; param:",  nom,
          "; value:", config[[nom]], config=config)
      }
      if(is.na(config[[nom]])) {
        f.err("check_config: param is NA; param:", nom, config=config)
      }
    }
  }
  
  for(nom in scalar_formula) {
    if(nom %in% names(config)) {
      if(!f.is_formula(config[[nom]])) {
        f.err("check_config: param not a formula; param:",  nom,
          "; value:", config[[nom]], config=config)
      }
      ## checks operators and variables; throws error if unsupported:
      f.parse_frm(config[[nom]], config)
    }
  }
  
  for(nom in vector_character) {
    if(nom %in% names(config)) {
      if(!is.character(config[[nom]])) {
        f.err("check_config: param not vector of character; param:",  nom, 
          "; value:", config[[nom]], config=config)
      }
    }
  }
  
  for(nom in vector_props) {
    if(nom %in% names(config)) {
      if(!is.numeric(config[[nom]])) {
        f.err("check_config: param not vector of proportions; param:",  nom, 
          "; value:", config[[nom]], config=config)
      }

      ## one NA_real_ or NaN anywhere in the vector made any(x < 0) itself NA,
      ##   which is R's "missing value where TRUE/FALSE needed" (measured).
      ##   An infinite element is left to the range check below, which names
      ##   it as the out-of-range proportion it is:

      if(anyNA(config[[nom]])) {
        f.err("check_config: param has an NA or NaN; param:", nom,
          "; value:", config[[nom]], config=config)
      }

      if(any(config[[nom]] < 0) || any(config[[nom]] > 1)) {
        f.err("check_config: proportions out of range [0, 1]; param:", nom,
          "; value:", config[[nom]], config=config)
      }
    }
  }
  
  for(nom in list_character) {
    if(nom %in% names(config)) {
      if(!is.list(config[[nom]])) {
        f.err("check_config: param not list of character; param:",  nom,
          "; value:", config[[nom]], config=config)
      }
    }
  }

  ## estimability is an ordered requirement, each level strictly stronger than
  ##   the one before it; see filter_features_by_estimability():

  if("estimability" %in% names(config)) {
    allowed <- c("test", "term", "full")
    if(!(config$estimability %in% allowed)) {
      f.err("check_config: unexpected estimability:", config$estimability, "\n",
        "allowed:", allowed, config=config)
    }
  }

  ## config$test_method names the engine test() will dispatch to, and a name that is
  ##   not one of them is a configuration error that no other setting can fix, so it is
  ##   caught here rather than at the end of the dispatch chain in test(). Empty is
  ##   allowed and means unset: test(method=) overrides config$test_method, and test()
  ##   refuses only when both are unset. NA is not, a missing setting being expressed
  ##   by its absence:

  if("test_method" %in% names(config)) {
    if(is.na(config$test_method) ||
        !(config$test_method %in% c(test_methods(), ""))) {
      f.err("check_config: unexpected test_method:", config$test_method, "\n",
        "allowed:", test_methods(), "\n",
        "or \"\", meaning unset, with test(method=) naming the engine instead",
        config=config)
    }
  }

  ## config$normalization_method, config$impute_method and config$feature_aggregation name
  ##   the methods normalize(), impute() and combine_features() dispatch to, and a name that
  ##   is not one of them is settled by the configuration alone, so it is refused here
  ##   rather than at the step that would meet it. That matters most for a misspelled
  ##   impute_method: impute() is the last step of the default config$run_order, so the whole
  ##   pipeline runs before the typo is found. Empty is allowed and means unset, the step's
  ##   own method= argument naming the method instead; NA is not, a missing setting being
  ##   expressed by its absence. Only the name is checked: whether the method suits the data
  ##   it is given stays with the step, which is the one that has the data:

  method_params <- list(
    normalization_method=normalize_methods(),
    impute_method=impute_methods(),
    feature_aggregation=c("medianPolish", "robustSummary", "none")
  )

  for(nom in names(method_params)) {
    if(!(nom %in% names(config))) next
    if(is.na(config[[nom]]) || !(config[[nom]] %in% c(method_params[[nom]], ""))) {
      f.err("check_config: unexpected", paste0(nom, ":"), config[[nom]], "\n",
        "allowed:", method_params[[nom]], "\n",
        "or \"\", meaning unset, with the step's method= argument naming the method",
        "instead", config=config)
    }
  }

  ## config$run_order names the pipeline steps run() walks, fetching each with get(), so
  ##   a name that is not one of them fails as "object not found" at the point that step
  ##   would have run, with every step before it already done and its output already
  ##   written. The set is small and fixed, so a misspelling is settled by the
  ##   configuration alone and belongs here. A zero length run_order is allowed and means
  ##   load_data() then test(), which is a real workflow for data that arrives prepared.
  ##   A repeated step is not refused, normalizing again after aggregation being something
  ##   someone may mean; run() warns about it once, f.save_state() naming its files after
  ##   the first occurrence:

  if("run_order" %in% names(config) && length(config$run_order)) {

    steps <- f.run_order_steps()
    bad <- config$run_order[is.na(config$run_order) |
      !(config$run_order %in% steps)]

    if(length(bad)) {
      f.err("check_config: unexpected step in config$run_order:", bad, "\n",
        "allowed:", steps, "\n",
        "  run() fetches each step by name, so this would otherwise fail partway",
        "through the workflow rather than before it starts;", "\n",
        "  config$run_order:", config$run_order, config=config)
    }
  }

  ## config$test_term and config$contrast are two different hypotheses about the
  ##   same model, and one run tests one hypothesis: the results file carries one
  ##   row per feature, and run.R's permutation aggregation reads it back that way.
  ##   Which of the two was tested cannot be recovered from the config if both are
  ##   set, and neither is a safe default to prefer: new_config() ships a non-empty
  ##   test_term, so preferring the contrast would silently ignore a populated key,
  ##   and preferring test_term would silently ignore the key that was deliberately
  ##   added. So set config$test_term to "" to test a contrast:

  if(f.contrast_set(config) && length(config$test_term) %in% 1 &&
      !is.na(config$test_term) && nzchar(trimws(config$test_term))) {
    f.err("check_config: config$contrast and config$test_term both name",
      "something to test, and one run tests one hypothesis;", "\n",
      "config$contrast:", config$contrast, "; config$test_term:",
      config$test_term, "\n",
      "set config$test_term to \"\" to test the contrast, or config$contrast to",
      "\"\" to test the term", config=config)
  }

  ## config$test_ridge=TRUE is msqrob2::msqrobAggregate(ridge=TRUE), which refits every
  ##   non-intercept column of the mean model as a level of a random effect and refuses
  ##   a mean model with fewer than two of them, so ~grp for a two level factor has to
  ##   be written ~0+grp; see test_msqrob(). How many columns a factor contributes
  ##   depends on how many levels it has, which is data and not configuration, so only
  ##   what the config settles is refused here: a formula with no non-intercept term at
  ##   all, and a single term whose type and levels initialize() has already resolved
  ##   into config$covariate_types and config$factor_levels. Anything less definite is
  ##   left to msqrob2. Only for the one method that reads the setting; on any other it
  ##   is an ignored setting, which the note below covers:

  if(isTRUE(config$test_ridge) && "test_method" %in% names(config) &&
      config$test_method %in% "msqrob_agg") {

    n_cols <- f.frm_min_noint_cols(config)

    if(!is.na(n_cols) && n_cols < 2) {
      f.err("check_config: config$test_ridge is TRUE, which penalizes the fixed",
        "effects and so needs a mean model with at least two non-intercept columns,",
        "and config$frm", config$frm, "has", n_cols, ";", "\n",
        if(n_cols %in% 0) {
          paste("  to fix, give config$frm a covariate for the penalty to apply to,",
            "or set config$test_ridge FALSE")
        } else {
          paste("  to fix, suppress the intercept, config$frm = ~0 + ..., which makes",
            "a two level factor two columns rather than one, or add a covariate, or",
            "set config$test_ridge FALSE")
        },
        config=config)
    }
  }

  ## reference_levels holds one reference level per factor variable in
  ##   config$frm, so every element needs a variable name, and every value has
  ##   to be a usable level:

  for(nom in c("reference_levels", "covariate_types", "factor_levels")) {

    if(!(nom %in% names(config)) || length(config[[nom]]) %in% 0) next

    noms1 <- names(config[[nom]])

    if(is.null(noms1) || any(is.na(noms1)) || any(nchar(noms1) %in% 0)) {
      f.err("check_config:", nom, "needs a variable name for every element;",
        "names:", noms1, config=config)
    }

    if(any(duplicated(noms1))) {
      f.err("check_config: duplicated variable names in", nom, ":",
        noms1[duplicated(noms1)], config=config)
    }

    if(nom %in% "factor_levels") next

    if(any(is.na(config[[nom]])) || any(nchar(config[[nom]]) %in% 0)) {
      f.err("check_config:", nom, "values must be non-empty and non-NA;",
        "values:", paste(noms1, config[[nom]], sep="="), config=config)
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
#'   See documentation for \code{h0testr::new_config()} 
#'     for more detailed description of configuration parameters. 
#' @param config List with configuration values.
#' @return NULL
#' @examples
#' config <- new_config()
#' 
#' ## you must customize frm, test_term, and reference_levels:
#' config$frm <- ~ age + sex + age:sex
#' config$test_term <- "age:sex"
#' config$reference_levels <- c(age="young", sex="female")
#'
#' ## these may or may not need customization, depending on your file formats:
#' config$feat_id_col <- "peptide_id"
#' config$gene_id_col <- "gene_id"
#' config$sample_id_col <- "sample"
#' config$obs_id_col <- "observation"
#'
#' report_config(config)

report_config <- function(config) {

  check_config(config)

  for(k1 in names(config)) {
    v1 <- config[[k1]]
    if(is.list(v1)) {
      for(k2 in names(v1)) {
        f.msg(k1, ":", k2, ":", paste(as.character(v1[[k2]]), collapse=", "), config=config)
      }
    } else if(!is.null(names(v1))) {
      ## named vectors (e.g. config$covariate_types) reported as name=value:
      f.msg(k1, ":", paste(names(v1), as.character(v1), sep="="), config=config)
    } else if(f.is_formula(v1)) {
      ## as.character() of a formula splits it at the tilde, which comma-joining
      ##   would then report as "~, age + sex"; deparse() keeps it one expression:
      f.msg(k1, ":", paste(deparse(v1), collapse=" "), config=config)
    } else f.msg(k1, ":", paste(as.character(v1), collapse=", "), config=config)
  }
  
  ## check config$test_term compatible with config$frm; throws error if not,
  ##   else returns NULL. A run testing config$contrast has no config$test_term to
  ##   check, check_config() having just refused a config that sets both; the
  ##   contrast is checked against the coefficients of the design by
  ##   f.design_contrast(), which needs the observations to build it:

  if(!f.contrast_set(config)) f.normalize_terms(config)
}
