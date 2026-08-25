## The pipeline steps run() will walk, which is what config$run_order may name. run()
##   fetches each by get(), so a name that is not one of these is an "object not found"
##   at the point the step would have run, with the steps before it already done; kept in
##   one place here so that check_config(), which refuses such a name up front, and the
##   comment on config$run_order in new_config() cannot drift apart. Not the order they
##   have to appear in: normalize() before combine_features() is required, aggregation
##   needing log scale data, but that is combine_features()'s own refusal to make, and
##   filter() and impute() are usable either way round:

f.run_order_steps <- function() {
  return(
    c("normalize", "combine_replicates", "combine_features", "filter", "impute")
  )
}

#' Run a basic workflow
#' @description
#'   Run a basic workflow according to: \code{config$run_order}.
#' @details
#'   Run a basic workflow: \code{load_data() -> config$run_order -> test()},
#'     where \code{config$run_order} is vector of functions which are run in
#'     the specified order.
#'   Each name in \code{config$run_order} must be one of \code{"normalize"},
#'     \code{"combine_replicates"}, \code{"combine_features"}, \code{"filter"} and
#'     \code{"impute"}; \code{h0testr::check_config()} refuses anything else before the
#'     workflow starts. A zero length \code{config$run_order} runs
#'     \code{h0testr::load_data()} and then \code{h0testr::test()}. Naming a step twice
#'     runs it twice, which is allowed but warned about, since \code{f.save_state()} names
#'     its output files after the first occurrence of a step and the later one overwrites
#'     them.
#'   Unlike \code{h0testr::tune()}, which skips \code{combine_features()} for the test
#'     methods that take feature level input, this function walks
#'     \code{config$run_order} as given: for \code{config$test_method} of
#'     \code{"prolfqua_lmer"} or \code{"msqrob_agg"}, leave \code{"combine_features"} out
#'     of it, or those methods will refuse the aggregated state at the last step.
#'   \code{config$test_method="none"} runs the whole workflow and skips the test step, so
#'     that this function can be used to normalize, combine, filter and impute and stop
#'     there; \code{$original}, \code{$standard} and \code{$fit} of the result are then
#'     \code{NULL} and \code{$state} is the processed state. \code{"none"} is not one of
#'     \code{h0testr::test_methods()}, which names the engines;
#'     \code{h0testr::check_config()} is where the values the key may hold are listed.
#'   See documentation for \code{h0testr::new_config()}
#'     for more detailed description of configuration parameters.
#' @param config List with configuration values like those returned by 
#'   \code{new_config()}.
#' @return A list with the following elements: 
#'   \tabular{ll}{
#'     \code{state}    \cr \tab List with elements \code{$expression}, \code{$features}, and \code{$samples}. \cr
#'     \code{config}   \cr \tab List with configuration settings. \cr
#'     \code{original} \cr \tab A \code{data.frame} with native results of test. \cr
#'     \code{standard} \cr \tab A \code{data.frame} with results in standardized format. \cr
#'     \code{fit}      \cr \tab Fitted model from selected testing procedure.
#'   }
#'   The last three are \code{NULL} when \code{config$test_method} is \code{"none"}.
#' @examples
#' config <- h0testr::new_config()          ## defaults
#' config$save_state <- FALSE               ## default is TRUE
#' config$dir_in <- system.file("extdata", package="h0testr")  ## where example data 
#' config$feature_file_in <- "features.tsv"
#' config$sample_file_in <- "samples.tsv"
#' config$data_file_in <- "expression.tsv" 
#' config$feat_id_col <- config$gene_id_col <- "feature_id"
#' config$obs_id_col <- config$sample_id_col <- "observation_id"
#' config$frm <- ~condition
#' config$test_term <- "condition"
#' config$test_method <- "trend"
#' config$reference_levels <- c(condition="placebo")
#' config$n_features_min <- 10     ## default 1000 too big for small demo dataset
#' config$run_order <- c("normalize", "combine_replicates", "filter", "impute")
#'
#' print(config$run_order)
#'
#' result <- h0testr::run(config)     ## run workflow
#' head(result$original)              ## hit table as returned by underlying software
#' head(result$standard)              ## hit table in a standardized format
#' print(result$fit)                ## model fit by selected testing procedure

run <- function(config) {

  ## no report_config(config) here: load_data() below calls it as its first statement,
  ##   so a call here only prints the same configuration to the log twice:

  f.log_block("starting load_data", config=config)
  out <- load_data(config)

  ## a step named twice runs twice, which is not obviously wrong: normalizing again after
  ##   aggregation, or filtering again after imputation, are things someone may mean. But
  ##   f.save_state() names its output files after the first occurrence of a step in
  ##   config$run_order, so a second occurrence overwrites the first one's files. Said
  ##   here, once, rather than in check_config(), which every step calls; the names
  ##   themselves are refused there, where a config error belongs:

  dups <- unique(config$run_order[duplicated(config$run_order)])

  if(length(dups)) {
    f.msg("WARNING: run: config$run_order names", dups, "more than once, so",
      if(length(dups) > 1) "those steps run" else "that step runs", "again;", "\n",
      " f.save_state() names its files after the first occurrence of a step, so the",
      "later one overwrites the files of the earlier;", "\n",
      " config$run_order:", config$run_order, config=config)
  }

  ## the scale of the data lives in out$config$is_log_transformed, put there by
  ##   initialize() and updated by normalize(), so each step below reads it from
  ##   the config it is handed rather than from a second copy kept here:

  for(f_name in config$run_order) {

    f.log_block("starting", f_name, config=config)
    fn <- get(f_name)
    out <- fn(out$state, out$config)
  }

  f.log_block("starting test", config=out$config)
  result <- test(out$state, out$config)
  
  return(list(state=out$state, config=out$config, 
    original=result$original, standard=result$standard, fit=result$fit))
}

## helper for tune(); normalize and combine reps:

f.tune1 <- function(state, config, normalization_method) {
  
  if(normalization_method %in% "q50") {
    config$normalization_method <- "quantile"
    config$normalization_quantile <- 0.5
  } else if(normalization_method %in% c("q75")) {
    config$normalization_method <- "quantile"
    config$normalization_quantile <- 0.75
  } else if(normalization_method %in% "upperquartile") {
    config$normalization_quantile <- 0.75
  }
  
  ## normalize and combine reps:
  f.log_block("f.tune:1: normalize", config=config)
  out <- normalize(state, config)
  
  f.log_block("f.tune:1: combine_replicates", config=config)
  out <- combine_replicates(out$state, out$config)
  
  f.log_block("f.tune:1: return", config=config)
  return(out)
}

## helper for tune(); filter, impute, and test:

## helper for f.tune2(); the result row for a parameter combination that was not
##   tested, so that one unusable combination does not abort a whole sweep. Same
##   columns as a real result row, with nhits and ntests unset:

f.tune2_na_row <- function(config) {
  return(
    data.frame(norm=config$normalization_method, nquant=config$normalization_quantile,
      impute=config$impute_method, iquant=config$impute_quantile,
      scale=config$impute_scale, span=config$impute_span,
      npcs=config$impute_npcs, k=config$impute_k, test=config$test_method,
      perm=config$permute_var, nhits=NA, ntests=NA,
      time=format(Sys.time(), "%H:%M:%S"), stringsAsFactors=F
    )
  )
}

f.tune2 <- function(state, config) {

  ## no is_log_transformed argument: the scale of the data lives in
  ##   config$is_log_transformed, put there by initialize() and updated by normalize(),
  ##   and impute() and test() below both resolve an unset argument from the config they
  ##   are handed. tune() was reading that same field out of the config and passing it back
  ##   in beside the config, so the two could only ever agree, and a second copy of a
  ##   setting is a second thing to keep in step for no gain:

  ## there were two guards here for engines that could not express the hypothesis
  ##   config$test_term implied, which skipped the combination rather than letting it
  ##   stop the sweep. Neither is left. The first was for an engine bounded by terms
  ##   rather than by coefficients: prolfqua, while it read the rows of a per-term anova
  ##   table; test_prolfqua() now compares an explicit full and reduced design, and the
  ##   guard went with f.test_max_terms(), which by then returned Inf for every method.
  ##   The second was for an engine that could report only one coefficient at a time.
  ##   deqms was the last of those, and f.deqms_moderated_f() forms the joint test from
  ##   the variance prior DEqMS::spectraCounteBayes() fits, so no engine is bounded that
  ##   way either and f.test_max_cols() and f.design_test_cols_max() went with it. Every
  ##   method now runs every formula and test_term that the filters accept, so a sweep
  ##   has nothing to skip on these grounds:

  f.log_block("f.tune:2: filter", config=config)
  out <- filter(state, config)

  if(length(unique(out$state$samples[[out$config$sample_id_col]])) < 4) {
    f.msg("WARNING: f.tune2: post-filter <4 samples left; return NA",
      config=config)
    return(f.tune2_na_row(config))
  }

  if(length(unique(out$state$features[[out$config$gene_id_col]])) < 20) {
    f.msg("WARNING: f.tune2: post-filter <20 genes left; return NA",
      config=config)
    return(f.tune2_na_row(config))
  }

  f.log_block("f.tune:2: impute", config=config)
  out <- impute(out$state, out$config)

  ## test_deqms() refuses a run in which every gene has the same number of features,
  ##   there being no spread for its variance prior to be fitted against. Unlike the
  ##   guards above this is a property of the data reaching the test rather than of
  ##   config$frm, so it is checked here, after filtering, and only for the one method
  ##   it applies to. Checked rather than caught so that the reason is specific:

  if(config$test_method %in% "deqms") {
    counts <- f.gene_counts(out$state, out$config, "f.tune2")
    if(length(unique(counts)) < 2) {
      f.msg("WARNING: f.tune2: test_method deqms needs the number of features per",
        "gene to vary, and every one of", length(counts), "genes has",
        unique(counts), "; skipping and returning NA", config=config)
      return(f.tune2_na_row(config))
    }
  }

  ## the mixed model paths of test_prolfqua() and test_msqrob() need feature level
  ##   input, which is what tune() hands them by not aggregating for a gene level
  ##   method. A config naming one column as both the feature and the gene id has no
  ##   feature level to model, and no other parameter in the sweep changes that:

  if(config$test_method %in% c("prolfqua_lmer", "msqrob_agg") &&
      out$config$feat_id_col %in% out$config$gene_id_col) {
    f.msg("WARNING: f.tune2: test_method", config$test_method, "needs feature level",
      "input, and config$feat_id_col and config$gene_id_col both name '",
      out$config$feat_id_col, "', so there is no feature level to model;",
      "skipping and returning NA", config=config)
    return(f.tune2_na_row(config))
  }

  ## a sweep exists to fill in a matrix of parameter combinations, so one combination
  ##   that an engine cannot fit should cost that cell and not the rest of the run. The
  ##   guards above catch the failures that can be predicted; anything else is caught
  ##   here, logged in full so that the run says where it ran into trouble, and recorded
  ##   as a combination that was never tested. tune_check() already reads such a row as
  ##   untested rather than as one that found nothing. f.err() logs its whole message
  ##   before stopping, so the detail is in the log immediately above the warning below,
  ##   which only records which combination the log entry belongs to:

  f.log_block("f.tune:2: test", config=config)
  result <- try(test(out$state, out$config), silent=T)

  if(inherits(result, "try-error")) {
    f.msg("WARNING: f.tune2: test_method", config$test_method, "failed on this",
      "combination; skipping and returning NA;", "\n",
      "  ", as.character(result), config=config)
    return(f.tune2_na_row(config))
  }

  tbl <- result$standard

  result <- data.frame(norm=config$normalization_method, nquant=config$normalization_quantile, 
    impute=config$impute_method, iquant=config$impute_quantile, 
    scale=config$impute_scale, span=config$impute_span, 
    npcs=config$impute_npcs, k=config$impute_k, test=config$test_method, 
    ## na.rm, so that a single feature with an undefined adjusted p-value, from a
    ##   singular or non-converged per-feature fit, counts as no hit rather than
    ##   turning the whole row's nhits into NA, which tune_check() then reads as a
    ##   combination that ran and found nothing:
    perm=config$permute_var, nhits=sum(tbl$adj_pval < 0.05, na.rm=T), ntests=nrow(tbl),
    time=format(Sys.time(), "%H:%M:%S"), stringsAsFactors=F)
  
  f.log_block("f.tune:2: return", config=config)
  return(result)
}

#' Basic tuning loop
#' @description
#'   Run a basic tuning loop according.
#' @details
#'   Run a basic tuning loop that iterates over parameter combinations, 
#'     performing the following:
#'     \tabular{l}{
#'       1. Read data, prefilter, and optionally permute with \code{load_data()}. \cr
#'       1. Inter-observation normalization with \code{normalize()}. \cr
#'       2. Combine replicate observations wtih \code{combine_replicates()}. \cr
#'       3. Combine peptides into gene/protein groups with \code{combine_features()}. \cr
#'       4. Filter uninformative features and observations with \code{filter()}. \cr
#'       5. Impute missing values with \code{impute()}. \cr
#'       6. Hypothesis testing with \code{test()}. \cr
#'     }
#'   Normally, one does one run with \code{config$permute_var=""}, and 
#'     \code{N} runs (we recommend \code{N >= 20}) with 
#'     \code{config$permute_var} set to the name of a variable in 
#'     \code{config$test_term}.
#'   Tunes the following \code{config} values: \code{normalization_method},
#'     \code{normalization_quantile}, \code{impute_method}, \code{impute_quantile},
#'     \code{impute_scale}, \code{impute_span}, \code{impute_k}, \code{impute_npcs}, and
#'     \code{test_method}. Every other setting is taken from \code{config} as given and held
#'     fixed for the whole sweep, so a sweep says nothing about it:
#'     \code{normalization_span}, \code{impute_floor_offset}, \code{impute_alpha},
#'     \code{impute_n_pts}, \code{impute_aug_steps}, \code{feature_aggregation},
#'     \code{n_samples_min}, \code{n_features_min}, \code{estimability},
#'     \code{df_resid_min}, \code{test_prior_df}, \code{test_moderate}, \code{test_trend},
#'     \code{test_random_obs} and \code{test_ridge}. To compare two values of one of those,
#'     run one sweep per value.
#'   \code{config$run_order} is not read at all: the step sequence above is fixed. A
#'     \code{config$run_order} that reorders or omits steps, which \code{h0testr::run()}
#'     honors, therefore sweeps a different pipeline than the one it names. The one
#'     departure from the sequence is \code{combine_features()}, which is skipped for the
#'     test methods that take feature level input; see \code{h0testr::run()}.
#'   See documentation for \code{h0testr::new_config()}
#'     for more detailed description of configuration parameters.
#' @param config List with configuration values like those returned by \code{new_config()}.
#' @param normalization_methods Character vector of methods to try. One or more element of
#'   \code{h0testr::normalize_methods()}, or \code{"q50"} or \code{"q75"}, which are
#'   \code{"quantile"} at a \code{normalization_quantile} of \code{0.5} and \code{0.75}.
#'   Defaults to every method \code{h0testr::normalize_methods()} names, with
#'   \code{"quantile"} entered as those two. \code{"loess"} is the slowest of them, so
#'   drop it from the list if the sweep takes too long. A name outside that set is refused
#'   before the sweep starts, rather than surfacing partway through from inside
#'   \code{normalize()}, which reports the argument it was handed rather than the sweep that
#'   handed it over. \code{"none"} is accepted: not normalizing is something a sweep has a
#'   use for comparing against.
#' @param impute_methods Character vector of methods to try. One or more of:
#'   \code{c("sample_lod", "unif_sample_lod", "unif_global_lod", "rnorm_feature", "glm_binom", "loess_logit", "glmnet", "rf", "knn", "min_det", "min_prob", "qrilc", "bpca", "ppca", "svdImpute", "lls", "missforest", "none")}.
#'   A name that is not one of \code{h0testr::impute_methods()} is refused before the sweep
#'   starts, rather than three loops down, after every earlier combination has already been
#'   run. \code{"none"} is accepted here too, for the same reason.
#' @param impute_quantiles Numeric vector of quantiles to try for \code{impute_unif_*} methods. 
#'   One or more values between \code{0.0} and \code{1.0}.
#' @param impute_scales Numeric vector of scales to try for \code{impute_rnorm_feature}. 
#'   See \code{impute_rnorm_feature()} \code{scale}. parameter.
#' @param impute_spans Numeric vector of spans to try for \code{impute_loess_logit}.
#' @param impute_npcs Numeric vector of N PCs to try for \code{impute_method \%in\% c("bpca", "ppca", "svdImpute")}.
#' @param impute_ks Numeric vector of \code{k} to use for \code{impute_method \%in\% c("knn", "lls")}.
#' @param test_methods Character vector with one or more element of
#'   \code{h0testr::test_methods()}. Defaults to every method except
#'   \code{"prolfqua_lmer"} and \code{"msqrob_agg"}, each of which fits a mixed model
#'   per gene and so costs orders of magnitude more time per cell than the rest of the
#'   sweep put together; name either explicitly to include it. A name that is not one of
#'   \code{h0testr::test_methods()} is refused before the sweep starts, rather than after
#'   the combination it appears in has been normalized and imputed. That includes
#'   \code{"none"}, which \code{config$test_method} accepts as "skip the test step": a
#'   sweep over not testing has nothing to compare.
#' @return A data.frame with the following columns:
#'   \tabular{ll}{
#'     \code{norm}   \cr \tab Normalization method (character). \cr
#'     \code{nquant} \cr \tab Normalization quantile (numeric). \cr
#'     \code{impute} \cr \tab Imputation method (character). \cr
#'     \code{iquant} \cr \tab Imputation quantile (numeric). \cr
#'     \code{scale}  \cr \tab Imputation scale, for \code{impute_rnorm_feature} and
#'       \code{impute_qrilc} (numeric). \cr
#'     \code{span}   \cr \tab Imputation span, for \code{impute_loess_logit} (numeric). \cr
#'     \code{npcs}   \cr \tab Number of PCs, for \code{"bpca"}, \code{"ppca"} and
#'       \code{"svdImpute"} (numeric). \cr
#'     \code{k}      \cr \tab Number of neighbors, for \code{"knn"} and \code{"lls"}
#'       (numeric). \cr
#'     \code{test}   \cr \tab Test method (character). \cr
#'     \code{perm}   \cr \tab Permuted variable (character). \cr
#'     \code{nhits}  \cr \tab Number of hits (numeric); \code{NA} if not tested. \cr
#'     \code{ntests} \cr \tab Number of tests (numeric); \code{NA} if not tested. \cr
#'     \code{time}   \cr \tab Timestamp. \cr
#'   }
#'   A column for a parameter that this combination does not use is still filled in, from
#'     \code{config}, rather than left \code{NA}: \code{iquant} carries
#'     \code{config$impute_quantile} even where \code{impute} is \code{"rf"}. The nine
#'     columns from \code{norm} through \code{test} are what identify a combination, and
#'     \code{h0testr::tune_check()} joins the permuted results to the unpermuted ones on
#'     those nine by name, so renaming or dropping one of them there is refused rather
#'     than silently changing what is being compared.
#'   A combination that could not be tested yields a row with \code{nhits} and
#'     \code{ntests} set to \code{NA} rather than aborting the sweep, and the reason
#'     is written to \code{config$log_file}. This happens when too few samples or
#'     genes survive filtering. It no longer happens because \code{test_method} cannot
#'     express the test \code{config$test_term} implies: every method now tests as many
#'     coefficients jointly as \code{config$test_term} carries. \code{"deqms"} and
#'     \code{"msqrob"} were the two that could not, and both were bounded at their APIs
#'     rather than by their fits: \code{msqrob2::hypothesisTest()} returns one table per
#'     column of the contrast, and \code{DEqMS::spectraCounteBayes()} moderates one
#'     coefficient's t-statistic, while a joint test needs only the fitted coefficients
#'     and, for \code{"deqms"}, the per-gene variance that same function fits. Both now
#'     report a joint test computed from those; see \code{h0testr::test_msqrob()} and
#'     \code{h0testr::test_deqms()}.
#'     \code{"proda"} and \code{"prolfqua"} were bounded by neither,
#'     testing several coefficients jointly: \code{"proda"} by likelihood ratio and
#'     \code{"prolfqua"} by an F-test comparing the full design against the design
#'     with the tested columns removed, with the error variance moderated across
#'     features as the \code{limma}-based methods do.
#'     \code{"prolfqua_lmer"} is unbounded in the same way, testing the same
#'     coefficients by a Satterthwaite Wald F from one mixed model per gene, and
#'     \code{"msqrob_agg"} by a joint Wald F from one \code{msqrob2} mixed model per
#'     gene.
#'     \code{"deqms"} is also skipped when every gene has the same number of
#'     features, there being no spread for its variance prior to be fitted against;
#'     see \code{h0testr::test_deqms()}. \code{"prolfqua_lmer"} and
#'     \code{"msqrob_agg"} are skipped when
#'     \code{config$feat_id_col} and \code{config$gene_id_col} name one column, there
#'     being no feature level for them to model; see
#'     \code{h0testr::test_prolfqua()} and \code{h0testr::test_msqrob()}.
#'     Any other failure of the test step is caught the same way, so that one
#'     combination an engine cannot fit costs that cell and not the rest of the
#'     sweep; the error is written to \code{config$log_file} in full, immediately
#'     above the line naming the combination it belongs to.
#'     \code{h0testr::tune_check()} counts these \code{NA} \code{nhits} as zero hits,
#'     but gives such a row no \code{fdr}, so that a combination which never ran is
#'     not ranked above every combination that did.
#' @examples
#' ## set up configuration:
#' config <- h0testr::new_config()     ## defaults
#' config$save_state <- FALSE          ## default is TRUE
#' config$dir_in <- system.file("extdata", package="h0testr")  ## where example data 
#' config$feature_file_in <- "features.tsv"
#' config$sample_file_in <- "samples.tsv"
#' config$data_file_in <- "expression.tsv" 
#' config$feat_id_col <- config$gene_id_col <- "feature_id"
#' config$obs_id_col <- config$sample_id_col <- "observation_id"
#' config$n_features_min <- 10         ## default 1000 too big for small demo dataset
#' config$frm <- ~condition
#' config$test_term <- "condition"
#' config$test_method <- "trend"
#' config$reference_levels <- c(condition="placebo")
#' 
#' ## one run with unpermuted data:
#' config$permute_var <- ""            ## no permutation
#' set.seed(101)
#' out1 <- h0testr::tune(config,
#'   normalization_methods=c("RLE", "q75", "cpm", "log2"),
#'   impute_methods=c("sample_lod", "unif_sample_lod", "none"),
#'   impute_quantiles=c(0, 0.05, 0.1),
#'   test_methods=c("trend", "msqrob", "proda", "prolfqua")
#' )
#' ## write.table(out1, "0.condition.tune.tsv", quote=F, sep="\t", row.names=F)
#' 
#' ## one run with permuted data; run 20+ such runs w/ suffices 1:20:
#' config$permute_var <- "condition"   ## permute variable in test_term
#' set.seed(101)
#' out2 <- h0testr::tune(config,
#'   normalization_methods=c("RLE", "q75", "cpm", "log2"),
#'   impute_methods=c("sample_lod", "unif_sample_lod", "none"),
#'   impute_quantiles=c(0, 0.05, 0.1),
#'   test_methods=c("trend", "msqrob", "proda", "prolfqua")
#' )
#' ## write.table(out2, "1.condition.tune.tsv", quote=F, sep="\t", row.names=F)

tune <- function(
    config,  
    ## every method normalize_methods() names, with "quantile" entered as the "q50"
    ##   and "q75" that f.tune1() turns back into it, since a quantile is a second
    ##   parameter and the sweep varies one name at a time. "loess" was left out while
    ##   normalize() handed raw intensities to limma::normalizeCyclicLoess() and then
    ##   log2(x + 1)'d the negative fitted values it returns for the smallest of them
    ##   into NaN: 5707 of the 146841 measured values on the rdtc_seer2 protein
    ##   groups, so the cell was scored on a matrix the sweep had damaged. normalize()
    ##   now transforms before that fit instead, so nothing is lost, but it is still
    ##   the slowest method here; drop it from this list if the sweep takes too long.
    ##   "quantiles.robust" is left out for a different reason: normalize_mscoreutils()
    ##   refuses it on data with any missing value: it assigns values by rank within each
    ##   observation, so a gap comes back at that rank in every observation rather than in
    ##   the one it was missing from (see normalize_mscoreutils()), and
    ##   normalization comes here before imputation, so the sweep would stop at that
    ##   cell on any real matrix. Pass it explicitly to sweep it on a complete one:
    normalization_methods=c("RLE", "upperquartile", "q50", "q75", "cpm", "max",
      "sum", "div.mean", "div.median", "TMM", "TMMwsp", "vsn", "qquantile",
      "loess", "log2", "none"),
    impute_methods=c("sample_lod", "unif_sample_lod", "unif_global_lod", 
      "rnorm_feature", "glm_binom", "loess_logit", "glmnet", "rf", 
      "knn", "min_det", "min_prob", "qrilc", "bpca", "ppca", "svdImpute", 
      "lls", "missforest", "none"),
    impute_quantiles=c(0, 0.01, 0.05, 0.1), 
    impute_scales=c(1, 0.33, 0.1),
    impute_spans=c(0.25, 0.5, 0.75),
    impute_npcs=c(3, 5, 10),
    impute_ks=c(5, 10, 20), 
    test_methods=c("lm", "trend", "deqms", "msqrob", "proda", "prolfqua", "voom")) {

  ## the sweep assigns each of these to config$test_method in turn, so a name that is not
  ##   an engine is caught here rather than inside the loop, where it would fail only after
  ##   the normalization, aggregation and imputation of that combination had been computed.
  ##   "none" is a legal config$test_method, meaning skip the test step, and is refused
  ##   here for the same reason as a typo: a sweep over not testing measures nothing. The
  ##   argument shadows the function of the same name, which R resolves anyway, a call
  ##   looking only at function bindings:

  bad <- setdiff(test_methods, test_methods())

  if(length(bad)) {
    f.err("tune: unexpected test_methods:", bad, "\n",
      "allowed:", test_methods(), "\n",
      "config$test_method also accepts \"none\", meaning skip the test step, which a",
      "sweep has no use for, there being nothing to compare", config=config)
  }

  ## the same reason applies to the other two lists the sweep assigns from, and neither
  ##   was checked. A name that is not an impute_method reaches f.err() three loops down,
  ##   after every earlier combination has been run, and one that is not a normalization
  ##   method is not checked at all and surfaces from inside normalize() with a message
  ##   about the argument it was handed rather than about the sweep that handed it over.
  ##   The allowed normalization set is normalize_methods() plus "q50" and "q75", which
  ##   normalize() does not accept: f.tune1() maps each onto normalization_method
  ##   "quantile" with the matching config$normalization_quantile. "none" is left in both
  ##   sets, unlike for test_methods above: not normalizing and not imputing are both
  ##   settings a sweep has a use for comparing against:

  bad <- setdiff(impute_methods, impute_methods())

  if(length(bad)) {
    f.err("tune: unexpected impute_methods:", bad, "\n",
      "allowed:", impute_methods(), config=config)
  }

  norm_ok <- union(normalize_methods(), c("q50", "q75"))
  bad <- setdiff(normalization_methods, norm_ok)

  if(length(bad)) {
    f.err("tune: unexpected normalization_methods:", bad, "\n",
      "allowed:", norm_ok, config=config)
  }

  ## load data:
  f.log_block("loading data", config=config)
  out <- load_data(config)              ## overwritten at each iteration
  state1 <- out$state                   ## save for subsequent iterations
  config1 <- out$config                 ## save for subsequent iterations
  rslt <- NULL

  for(normalization_method in normalization_methods) {
    
    config1$normalization_method <- normalization_method
    f.msg("normalization_method:", normalization_method, config=config1)
    
    f.log_block("normalize and combine reps", config=config1)
    out <- f.tune1(state1, config1, normalization_method=normalization_method)
    state2 <- out$state                 ## save for subsequent iterations
    config2 <- out$config               ## save for subsequent iterations

    ## normalize(), called by f.tune1(), records the scale in config2, which is what
    ##   reaches f.tune2() below, so the sweep does not keep a second copy of it:

    for(test_method in test_methods) {
      
      config2$test_method <- test_method
      f.msg("test_method:", test_method, config=config2)
      
      if(!f.gene_level_method(test_method)) {
        ## for methods that do not use peptides for gene testing:
        f.log_block("combine_features", config=config2)
        out <- combine_features(state2, config2)
        state3 <- out$state
        config3 <- out$config
      } else {
        ## for methods that do use peptides for gene testing:
        f.log_block("skipping combine_features", config=config2)
        state3 <- state2
        config3 <- config2
      }
      
      for(impute_method in impute_methods) {
        
        config3$impute_method <- impute_method
        f.msg("impute_method:", impute_method, config=config3)
        
        if(impute_method %in% c("unif_global_lod", "unif_sample_lod", "min_det")) {
          for(impute_quantile in impute_quantiles) {
            f.log_block("normalization_method:", normalization_method, 
              "; impute_method:", impute_method, "; test_method:", 
              test_method, config=config3)
            f.msg("impute_quantile:", impute_quantile, config=config3)
            config3$impute_quantile <- impute_quantile
            
            f.log_block("filter, impute, and test", config=config3)
            rslt_i <- f.tune2(state3, config3)
            rslt <- rbind(rslt, rslt_i)
            f.log_obj(rslt, config=config3)
          }
        } else if(impute_method %in% c("qrilc", "rnorm_feature")) {
          for(impute_scale in impute_scales) {
            f.log_block("normalization_method:", normalization_method, 
              "; impute_method:", impute_method, "; test_method:", 
              test_method, config=config3)
            f.msg("impute_scale:", impute_scale, config=config3)
            config3$impute_scale <- impute_scale
            
            f.log_block("filter, impute, and test", config=config3)
            rslt_i <- f.tune2(state3, config3)
            rslt <- rbind(rslt, rslt_i)
            f.log_obj(rslt, config=config3)
          }
        } else if(impute_method %in% c("min_prob")) {
          for(impute_quantile in impute_quantiles) {
            for(impute_scale in impute_scales) {
              f.log_block("normalization_method:", normalization_method, 
                "; impute_method:", impute_method, 
                "; test_method:", test_method, config=config3)
              f.msg("impute_quantile:", impute_quantile, 
                "; impute_scale:", impute_scale, config=config3)
              config3$impute_quantile <- impute_quantile
              config3$impute_scale <- impute_scale
              
              f.log_block("filter, impute, and test", config=config3)
              rslt_i <- f.tune2(state3, config3)
              rslt <- rbind(rslt, rslt_i)
              f.log_obj(rslt, config=config3)
            }
          }
        } else if(impute_method %in% c("loess_logit")) {
          for(impute_span in impute_spans) {
            f.log_block("normalization_method:", normalization_method, 
              "; impute_method:", impute_method, "; test_method:", 
              test_method, config=config3)
            f.msg("impute_span:", impute_span, config=config3)
            config3$impute_span <- impute_span
            
            f.log_block("filter, impute, and test", config=config3)
            rslt_i <- f.tune2(state3, config3)
            rslt <- rbind(rslt, rslt_i)
            f.log_obj(rslt, config=config3)
          }
        } else if(impute_method %in% c("bpca", "ppca", "svdImpute")) {
          for(npcs in impute_npcs) {
            f.log_block("normalization_method:", normalization_method, 
              "; impute_method:", impute_method, "; test_method:", 
              test_method, config=config3)
            f.msg("npcs:", npcs, config=config3)
            config3$impute_npcs <- npcs
            
            f.log_block("filter, impute, and test", config=config3)
            rslt_i <- f.tune2(state3, config3)
            rslt <- rbind(rslt, rslt_i)
            f.log_obj(rslt, config=config3)
          }
        } else if(impute_method %in% c("knn", "lls")) {
          for(impute_k in impute_ks) {
            f.log_block("normalization_method:", normalization_method, 
              "; impute_method:", impute_method, "; test_method:", 
              test_method, config=config3)
            f.msg("impute_k:", impute_k, config=config3)
            config3$impute_k <- impute_k
            
            f.log_block("filter, impute, and test", config=config3)
            rslt_i <- f.tune2(state3, config3)
            rslt <- rbind(rslt, rslt_i)
            f.log_obj(rslt, config=config3)
          }
        } else if(impute_method %in% c("sample_lod", "glm_binom", "glmnet", 
            "rf", "missforest", "none")) {
          
          f.log_block("normalization_method:", normalization_method, 
            "; impute_method:", impute_method, "; test_method:", test_method, 
            config=config3)
          ## a guard here, labelled OIL_WATER, skipped test_method "msqrob" and "voom"
          ##   whenever state3$expression carried any NA. It is gone, for three reasons.
          ##   It read the pre-imputation matrix, so it fired on every combination in
          ##   this branch and not just the one where the NAs survive: all of sample_lod,
          ##   glm_binom, glmnet, rf and missforest fill them, and were dropped anyway.
          ##   Neither engine needs it even at impute_method "none": test_msqrob() counts
          ##   non-missing values per feature and hands them to msqrob2's weighting, so
          ##   missingness is what it is built for, and test_voom() holds out the features
          ##   carrying an NA, saying how many, and refuses only when fewer than two
          ##   complete features are left, from the post-filter matrix that actually
          ##   reaches it rather than from this one. And it advanced the loop with next
          ##   instead of recording anything, so the combination was absent from the
          ##   result rather than present with nhits NA, which is what every other skip
          ##   here produces and what tune_check() joins on. The one case that genuinely
          ##   fails, voom on a state imputation left too sparse, is caught by the try()
          ##   in f.tune2() and recorded with test_voom()'s own message:

          f.log_block("filter, impute, and test", config=config3)
          rslt_i <- f.tune2(state3, config3)
          rslt <- rbind(rslt, rslt_i)
          f.log_obj(rslt, config=config3)
        } else {
          f.err("tune: unexpected impute_method:", 
            impute_method, config=config3)
        } ## if impute_method %in% ...
      }   ## for impute_method in impute_methods
    }     ## for test_method in test_methods
  }       ## for normalization_method in normalization_methods
  
  ## rslt starts NULL and grows by rbind(), so a sweep that never reached f.tune2()
  ##   returns NULL rather than an empty table, and the write.table() in the usage above
  ##   errors on it. An empty method list does that, and so does any future skip that
  ##   advances the loop without recording a row. Columns come from f.tune2_na_row() so
  ##   that the empty table cannot drift out of step with the rows a real sweep returns:

  if(is.null(rslt)) {
    f.msg("WARNING: tune: no parameter combination was tested; returning an empty",
      "result table", config=config)
    rslt <- f.tune2_na_row(config)[0, , drop=F]
  }

  f.log_block("returning result", config=config)
  return(rslt)
}

#' Check tuning results
#' @description
#'   Check results of tuning.
#' @details
#'   Imports data from basic tuning loop, comparing results from unpermuted
#'     data with those from permuted data. FDR is estimated from the permuted
#'     data results. Recommend that tuning use at least 20 iterations with
#'     permuted data.
#'   Combinations that \code{h0testr::tune()} skipped, because the test method
#'     cannot run the test \code{config$test_term} names, and combinations that
#'     lost every feature to filtering, have \code{ntests} of \code{0} and are
#'     given an \code{NA} \code{fdr}, which sorts them below every combination
#'     that ran. They are not combinations that found nothing; they are
#'     combinations that were never tested.
#'   A combination present in the unpermuted results but in none of the
#'     permuted ones has nothing to estimate its FDR from, so it too is given
#'     an \code{NA} \code{fdr} and sorted below the combinations that have one.
#'     A warning names how many combinations each set has that the other does
#'     not, in both directions.
#'   Permuted result rows with \code{ntests} of \code{0} were skipped as well,
#'     so they are dropped before \code{max1}, \code{mid1}, \code{avg1} and
#'     \code{sd1} are taken, rather than counted as permutations that found no
#'     false positives. A combination left with no permuted run that tested
#'     anything is then treated as one with no permuted counterpart, and gets
#'     an \code{NA} \code{fdr} rather than an \code{fdr} of \code{0}.
#'   A \code{dir_in} holding the unpermuted file and no permuted one has nothing
#'     to estimate an FDR against at all, and is refused with a message naming the
#'     prefix and suffix that were looked for.
#' @param dir_in Character scalar with path to directory containing tuning results.
#' @param prefix Character scalar with prefix (if any) of tuning result filenames.
#'   Only files beginning with it are read, so that two sweeps sharing
#'   \code{suffix} and \code{dir_in} do not pool their results; a warning names
#'   how many files were passed over. \code{""} reads every file ending in
#'   \code{suffix}.
#' @param suffix Character scalar with distinctive suffix (required) of tuning results filenames.
#' @param config List with at least \code{log_file} defined (can be \code{""}).
#' @param fdr_cutoff Numeric scalar between \code{0} and \code{1.0} specifying 
#'   cutoff for false discovery rate. Trials not meeting cutoff are moved to the 
#'   bottom of the output \code{data.frame}.
#' @return A \code{data.frame} with the following columns:
#'   \tabular{ll}{
#'     \code{nhits}      \cr \tab Number of significant hits. \cr
#'     \code{ntests}     \cr \tab Number of features tested; \code{0} if the combination did not run. \cr
#'     \code{fdr}        \cr \tab False discovery rate; \code{NA} where \code{ntests}
#'       is \code{0}, or where there is no permuted counterpart. \cr
#'     \code{max1}       \cr \tab Maximum number of hits in any permutation. \cr
#'     \code{mid1}       \cr \tab Median number of hits across permutations. \cr
#'     \code{avg1}       \cr \tab Average number of hits across permutations. \cr
#'     \code{sd1}        \cr \tab Standard deviation of number of hits across permutations. \cr
#'     \code{norm}       \cr \tab Normalization method. \cr
#'     \code{nquant}     \cr \tab Normalization quantile. \cr
#'     \code{impute}     \cr \tab Imputation method. \cr
#'     \code{iquant}     \cr \tab Imputation quantile. \cr
#'     \code{scale}      \cr \tab Scale for imputation. \cr
#'     \code{span}       \cr \tab Span for loess-based imputation. \cr
#'     \code{npcs}       \cr \tab Number of principle components for imputation. \cr
#'     \code{k}          \cr \tab Number of nearest neighbors or groups for imputation. \cr
#'     \code{test}       \cr \tab Test method. \cr
#'   }
#' @examples
#' dir_in <- system.file("extdata/tune", package="h0testr")
#' prefix <- ""
#' suffix <- ".condition.tune.tsv"
#' config <- list()
#' tbl <- h0testr::tune_check(dir_in, prefix, suffix, config)
#' print(tbl)

tune_check <- function(dir_in, prefix, suffix, config, fdr_cutoff=0.05) {

  ## prefix was taken and then ignored: the pattern was built from suffix alone, so every
  ##   file in dir_in ending in suffix was treated as part of this sweep, and the permuted
  ##   files of a second sweep sharing the suffix were pooled into this one's null without
  ##   a word about it. Anchor prefix at the front of the name too. A prefix of "" still
  ##   matches every name ending in suffix, as before. The suffix-only pattern is kept
  ##   because it is also what strips the suffix off a filename below:

  sfx_re <- gsub("(\\W)", "\\\\\\1", suffix)
  prfx_re <- gsub("(\\W)", "\\\\\\1", prefix)
  pat_sfx <- paste0(sfx_re, "$")
  pat <- paste0("^", prfx_re, ".*", pat_sfx)

  files_sfx <- sort(list.files(path=dir_in, pattern=pat_sfx))
  files <- files_sfx[grepl(pat, files_sfx)]
  other <- setdiff(files_sfx, files)

  if(length(other)) {
    f.msg("WARNING: tune_check:", length(other), "file(s) in", dir_in, "end in", suffix,
      "but do not begin with the prefix", paste0("'", prefix, "',"),
      "so they are another sweep's results and are ignored; first few:",
      utils::head(other, 3), config=config)
  }

  unperm_file <- paste0(prefix, "0", suffix)
  i0 <- files %in% unperm_file
  if(sum(i0) != 1) f.err("tune_check: no unperm file found; looking for:", 
    unperm_file, config=config)
  perm_files <- files[!i0]
  f.msg("found 1 unpermuted file and", length(perm_files), "permuted files\n",
    config=config)

  ## with no permuted file there is no null to estimate an fdr against, and nothing
  ##   downstream said so in terms anyone could act on. do.call(rbind, list()) below is
  ##   NULL, the ntests fill turns that into a one-element list, and the drop of untested
  ##   rows dies on it with "incorrect number of dimensions", which names neither this
  ##   function nor the files that are missing; before that drop existed the same case
  ##   reached the join column check instead and was reported as a missing column, which
  ##   was no better. Say it here, where the count is already in hand, and name what was
  ##   looked for, since a mistyped prefix or suffix is the likely cause:

  if(length(perm_files) < 1) {
    f.err("tune_check: found the unpermuted file", unperm_file, "but no permuted results",
      "file to estimate an fdr from; looked in", dir_in, "for other files beginning with",
      paste0("'", prefix, "'"), "and ending in", paste0("'", suffix, "';"),
      "h0testr::tune() writes those with config$permute_var set", config=config)
  }
  dat0 <- utils::read.table(paste(dir_in, unperm_file, sep="/"), header=T, 
    sep="\t", quote="", as.is=T)

  dat0$nhits[is.na(dat0$nhits)] <- 0
  dat0$ntests[is.na(dat0$ntests)] <- 0

  obj <- list()
  for(perm_file in perm_files) {
    f.msg("reading", perm_file, config=config)
    prfx <- sub(pat_sfx, "", perm_file)   ## pat now matches the whole name, so use pat_sfx
    dat_i <- utils::read.table(paste(dir_in, perm_file, sep="/"), header=T, 
      sep="\t", quote="", as.is=T)
    dat_i$perm_prfx <- prfx
    obj[[perm_file]] <- dat_i
  }

  dat1 <- do.call(rbind, obj)
  rownames(dat1) <- NULL

  dat1$ntests[is.na(dat1$ntests)] <- 0

  ## a permuted run that h0testr::tune() skipped, or that lost every feature to
  ##   filtering, has ntests 0 and nhits NA. Setting that nhits to 0 and summarizing it
  ##   beside the runs that ran says the permutation looked for false positives and found
  ##   none, when it never looked: a combination whose permuted runs were all skipped came
  ##   out of here with max1 0 and therefore an fdr of 0, the best score in the table,
  ##   sorted near the top on no evidence at all; one with only some of its permuted runs
  ##   skipped had max1, mid1, avg1 and sd1 pulled toward 0 by runs that never happened.
  ##   The unpermuted side carries the same rows and is already refused an fdr for them,
  ##   below. Drop them here instead, which leaves a combination with no permuted run that
  ##   tested anything looking like one with no permuted row at all; for estimating an fdr
  ##   the two are the same thing, and the mismatch warning below reports it:

  i_run <- !(dat1$ntests %in% 0)

  if(any(!i_run)) {
    f.msg("WARNING: tune_check:", sum(!i_run), "of", nrow(dat1),
      "permuted result row(s) have ntests 0, so nothing was tested in them; they are",
      "dropped rather than counted as permutations that found no hits", config=config)
  }

  dat1 <- dat1[i_run, , drop=F]

  if(nrow(dat1) < 1) {
    f.err("tune_check: no permuted result row has ntests above 0, so there is nothing",
      "to estimate an fdr from; permuted files read:", length(perm_files), config=config)
  }

  dat1$nhits[is.na(dat1$nhits)] <- 0

  ## the join key was dat0[, 1:9] and dat1[, 1:9], the first nine columns being the
  ##   parameter combination as f.tune2() happens to emit it. Add a column to that
  ##   data.frame, or reorder it, and the key silently changes meaning: it still pastes
  ##   nine values together and the two sides still match each other, just on the wrong
  ##   nine, so every fdr below is computed against the wrong permuted rows and nothing
  ##   errors. Naming the columns makes that a stop instead. Both tables are checked, and
  ##   separately: the permuted files are separate runs of the sweep and can have been
  ##   written by a different version of it than the unpermuted one:

  key_cols <- c("norm", "nquant", "impute", "iquant", "scale", "span", "npcs", "k",
    "test")

  bad <- setdiff(key_cols, names(dat0))

  if(length(bad)) {
    f.err("tune_check:", unperm_file, "is missing join column(s):", bad, "\n",
      "  columns found:", names(dat0), config=config)
  }

  bad <- setdiff(key_cols, names(dat1))

  if(length(bad)) {
    f.err("tune_check: permuted results are missing join column(s):", bad, "\n",
      "  columns found:", names(dat1), config=config)
  }

  k0 <- apply(dat0[, key_cols, drop=F], 1, paste, collapse=":")
  k1 <- apply(dat1[, key_cols, drop=F], 1, paste, collapse=":")

  ## perm_max[k0] below is NA for any unpermuted combination with no permuted
  ##   counterpart, so that combination gets an NA fdr and the sort files it with the
  ##   combinations that never ran, which is not what it is: it ran, and there is just
  ##   nothing to measure it against. Combinations whose permuted runs were all dropped
  ##   just above, having tested nothing, arrive here the same way and are counted with
  ##   them. A permuted combination absent from the unpermuted table goes the other way
  ##   and is dropped from the output without trace. Neither is an error, the two sides
  ##   being separate runs, but both mean the sweeps do not line up and both were silent:

  miss0 <- setdiff(unique(k0), unique(k1))
  miss1 <- setdiff(unique(k1), unique(k0))

  if(length(miss0)) {
    f.msg("WARNING: tune_check:", length(miss0), "of", length(unique(k0)),
      "unpermuted combination(s) have no permuted counterpart that tested anything, so",
      "they get no fdr and sort below every combination that has one; first few:",
      utils::head(miss0, 3), config=config)
  }

  if(length(miss1)) {
    f.msg("WARNING: tune_check:", length(miss1), "of", length(unique(k1)),
      "permuted combination(s) are absent from", unperm_file,
      "and are ignored; first few:", utils::head(miss1, 3), config=config)
  }

  ## permuted results in dat1; take max, median, mean, and sd of 10 permutation results:
  perm_max <- tapply(dat1$nhits, k1, max, na.rm=T)
  perm_mid <- tapply(dat1$nhits, k1, stats::median, na.rm=T)
  perm_avg <- tapply(dat1$nhits, k1, mean, na.rm=T)
  perm_sd  <- tapply(dat1$nhits, k1, stats::sd, na.rm=T)
  
  ## get them in the same order as dat1 (k1 made from dat1):
  dat0$max1 <- perm_max[k0]   ## max number of hits in permutations
  dat0$mid1 <- perm_mid[k0]   ## median number of hits in permutations
  dat0$avg1 <- perm_avg[k0]   ## mean number of hits in permutations
  dat0$sd1  <- perm_sd[k0]    ## sd(nhits) in permutations
  dat0$perm <- NULL

  ## average number of false positives == average number of hits across 10 sets of permutation results;
  ##   false positive rate: (average number of false positives) / (total number of positives)
  
  nhits <- dat0$nhits
  nhits[nhits %in% 0] <- 1
  dat0$fdr <- dat0$max1 / nhits   ## used to be $avg
  dat0$fdr[dat0$fdr > 1] <- 1.0

  ## a combination f.tune2() skipped, or that lost every feature to filtering, has
  ##   ntests 0 and nhits NA, which the substitution above turns into 0 hits out of
  ##   0 tests; its permuted runs are skipped identically, so max1 is 0 too and the
  ##   fdr computed for it is 0, which is the best score there is. Such a
  ##   combination is not a good one, it is one that never ran, so give it no fdr
  ##   and let the sort below put it below every combination that did run. ntests is
  ##   reported so the difference from 0 hits out of many tests is visible:

  dat0$fdr[dat0$ntests %in% 0] <- NA

  dat0 <- dat0[, c("nhits", "ntests", "fdr", "max1", "mid1", "avg1", "sd1", "norm",
    "nquant", "impute", "iquant", "scale", "span", "npcs", "k", "test")]

  i <- dat0$fdr < fdr_cutoff
  i[is.na(i)] <- FALSE
  dat0a <- dat0[i, ]
  dat0b <- dat0[!i, ]
  
  dat0a <- dat0a[order(dat0a$nhits, -dat0a$fdr, decreasing=T),  ]
  dat0b <- dat0b[order(dat0b$nhits, -dat0b$fdr, decreasing=T), ]

  dat0 <- rbind(dat0a, dat0b)
  rownames(dat0) <- NULL

  return(dat0)
}

