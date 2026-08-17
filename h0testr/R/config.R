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
    test_trend=FALSE,                    ## whether the prior variance of that shrinkage is fitted against mean feature intensity instead of being flat; used by test_prolfqua(), where it makes the prior the one test_method="trend" uses; unrelated to test_method
    test_random_obs=TRUE,                ## whether the feature level mixed model paths, test_method %in% c("prolfqua_lmer", "msqrob_agg"), add a random observation effect to the random feature effect; TRUE is the calibrated model; FALSE gives the feature-only structure both packages document, which is anti-conservative; see test_prolfqua() and test_msqrob()
    test_ridge=FALSE,                    ## whether test_method="msqrob_agg" penalizes the fixed effects, msqrob2::msqrobAggregate(ridge=TRUE); FALSE is msqrob2's own default throughout; TRUE shrinks the coefficients toward zero, so the reported logFC is not comparable to what the other methods report, renames the fitted parameters, and refuses a mean model with fewer than two non-intercept columns; see test_msqrob()

    ## run_order character vector with elements from {"normalize", "combine_replicates", "combine_features", "filter", "impute"}:
    run_order=c("normalize", "combine_replicates", "combine_features", "filter", "impute"),   ## order of workflow operations
    
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
      "negative values are an error")
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
        f.msg(k1, ":", k2, ":", paste(as.character(v1[[k2]]), sep=", "), config=config)
      }
    } else if(!is.null(names(v1))) {
      ## named vectors (e.g. config$covariate_types) reported as name=value:
      f.msg(k1, ":", paste(names(v1), as.character(v1), sep="="), config=config)
    } else f.msg(k1, ":", paste(as.character(v1), sep=", "), config=config)
  }
  
  ## check config$test_term compatible with config$frm; throws error if not,
  ##   else returns NULL. A run testing config$contrast has no config$test_term to
  ##   check, check_config() having just refused a config that sets both; the
  ##   contrast is checked against the coefficients of the design by
  ##   f.design_contrast(), which needs the observations to build it:

  if(!f.contrast_set(config)) f.normalize_terms(config)
}
