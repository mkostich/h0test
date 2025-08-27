#' Run a basic workflow
#' @description
#'   Run a basic workflow according to: \code{config$run_order}.
#' @details
#'   Run a basic workflow: \code{load_data() -> config$run_order -> test()}, 
#'     where \code{config$run_order} is vector of functions which are run in 
#'     the specified order.
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
#' config$sample_factors <- list(condition=c("placebo", "drug"))
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
  
  report_config(config)
  
  f.log_block("starting load_data", config=config)
  out <- load_data(config)
  is_log_transformed <- FALSE
  
  for(f_name in config$run_order) {
    
    f.log_block("starting", f_name, config=config)
    fn <- get(f_name)
    
    if(f_name %in% c("impute", "test")) {
      out <- fn(out$state, out$config, 
        is_log_transformed=is_log_transformed)
    } else {
      out <- fn(out$state, out$config)
    }
    
    if(f_name %in% "normalize") {
      if(!config$normalization_method %in% c("none")) {
        is_log_transformed <- TRUE
      }
    }
  }
  
  f.log_block("starting test", config=out$config)
  result <- test(out$state, out$config, is_log_transformed=is_log_transformed)
  
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

f.tune2 <- function(state, config, is_log_transformed=is_log_transformed) {
  
  f.log_block("f.tune:2: filter", config=config)
  out <- filter(state, config)
  
  if(length(unique(out$state$samples[[out$config$sample_id_col]])) < 4) {
    f.msg("WARNING: f.tune2: post-filter <4 samples left; return NA", 
      config=config) 
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
  
  if(length(unique(out$state$features[[out$config$gene_id_col]])) < 20) {
    f.msg("WARNING: f.tune2: post-filter <20 genes left; return NA", 
      config=config) 
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
  
  f.log_block("f.tune:2: impute", config=config)
  out <- impute(out$state, out$config, 
    is_log_transformed=is_log_transformed)
  
  f.log_block("f.tune:2: test", config=config)
  result <- test(out$state, out$config, 
    is_log_transformed=is_log_transformed)
  
  tbl <- result$standard
  
  result <- data.frame(norm=config$normalization_method, nquant=config$normalization_quantile, 
    impute=config$impute_method, iquant=config$impute_quantile, 
    scale=config$impute_scale, span=config$impute_span, 
    npcs=config$impute_npcs, k=config$impute_k, test=config$test_method, 
    perm=config$permute_var, nhits=sum(tbl$adj_pval < 0.05), ntests=nrow(tbl), 
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
#'   Currently tunes the following \code{config} values: \code{normalization_method}, 
#'     \code{normalization_quantile}, \code{impute_method}, \code{impute_quantile}, 
#'     \code{impute_scale}, \code{impute_span}, \code{impute_k}, 
#'     \code{impute_npcs}, and \code{test_method}. Notably, does not currently 
#'     tune \code{impute_alpha}, \code{impute_aug_steps}, or \code{run_order}. 
#'   See documentation for \code{h0testr::new_config()} 
#'     for more detailed description of configuration parameters. 
#' @param config List with configuration values like those returned by \code{new_config()}.
#' @param normalization_methods Character vector of methods to try. One or more of:
#'   \code{c("RLE", "upperquartile", "q50", "q75", "quantiles.robust", "cpm", "max", "div.mean", "TMMwsp", "vsn", "qquantile", "log2", "none")}.
#' @param impute_methods Character vector of methods to try. One or more of:
#'   \code{c("sample_lod", "unif_sample_lod", "unif_global_lod", "rnorm_feature", "glm_binom", "loess_logit", "glmnet", "rf", "knn", "min_det", "min_prob", "qrilc", "bpca", "ppca", "svdImpute", "lls", "missforest", "none")}.
#' @param impute_quantiles Numeric vector of quantiles to try for \code{impute_unif_*} methods. 
#'   One or more values between \code{0.0} and \code{1.0}.
#' @param impute_scales Numeric vector of scales to try for \code{impute_rnorm_feature}. 
#'   See \code{impute_rnorm_feature()} \code{scale}. parameter.
#' @param impute_spans Numeric vector of spans to try for \code{impute_loess_logit}.
#' @param impute_npcs Numeric vector of N PCs to try for \code{impute_method \%in\% c("bpca", "ppca", "svdImpute")}.
#' @param impute_ks Numeric vector of \code{k} to use for \code{impute_method \%in\% c("knn", "lls")}.
#' @param test_methods Character vector with one or more of: 
#'   \code{c("voom", "trend", "deqms", "msqrob", "proda")}.
#' @return A data.frame with the following columns:
#'   \tabular{ll}{
#'     \code{norm}       \cr \tab Normalization method (character). \cr
#'     \code{norm_quant} \cr \tab Normalization quantile (numeric). \cr
#'     \code{impute}     \cr \tab Imputation method (character). \cr
#'     \code{imp_quant}  \cr \tab Imputation quantile (numeric). \cr
#'     \code{scale}      \cr \tab Imputation scale for \code{impute_rnorm_feature}. \cr
#'     \code{test}       \cr \tab Test method (character). \cr
#'     \code{perm}       \cr \tab Permuted variable (character). \cr
#'     \code{nhits}      \cr \tab Number of hits (numeric). \cr
#'     \code{ntests}     \cr \tab Number of tests (numeric). \cr
#'     \code{time}       \cr \tab Timestamp. \cr
#'   }
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
#' config$sample_factors <- list(condition=c("placebo", "drug"))
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
    normalization_methods=c("RLE", "upperquartile", "q50", "q75", 
      "quantiles.robust", "cpm", "max", "div.mean", "TMMwsp", "vsn", 
      "qquantile", "log2", "none"),
    impute_methods=c("sample_lod", "unif_sample_lod", "unif_global_lod", 
      "rnorm_feature", "glm_binom", "loess_logit", "glmnet", "rf", 
      "knn", "min_det", "min_prob", "qrilc", "bpca", "ppca", "svdImpute", 
      "lls", "missforest", "none"),
    impute_quantiles=c(0, 0.01, 0.05, 0.1), 
    impute_scales=c(1, 0.33, 0.1),
    impute_spans=c(0.25, 0.5, 0.75),
    impute_npcs=c(3, 5, 10),
    impute_ks=c(5, 10, 20), 
    test_methods=c("trend", "deqms", "msqrob", "proda", "prolfqua", "voom")) {
  
  ## load data:
  f.log_block("loading data", config=config)
  out <- load_data(config)              ## overwritten at each iteration
  state1 <- out$state                   ## save for subsequent iterations
  config1 <- out$config                 ## save for subsequent iterations
  rslt <- NULL
  is_log_transformed <- FALSE
  
  for(normalization_method in normalization_methods) {
    
    config1$normalization_method <- normalization_method
    f.msg("normalization_method:", normalization_method, config=config1)
    
    f.log_block("normalize and combine reps", config=config1)
    out <- f.tune1(state1, config1, normalization_method=normalization_method)
    state2 <- out$state                 ## save for subsequent iterations
    config2 <- out$config               ## save for subsequent iterations
    if(normalization_method %in% c("none")) {
      is_log_transformed <- FALSE
    } else {
      is_log_transformed <- TRUE
    }
    
    for(test_method in test_methods) {
      
      config2$test_method <- test_method
      f.msg("test_method:", test_method, config=config2)
      
      if(!(test_method %in% c("deqms", "msqrob"))) {
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
            rslt_i <- f.tune2(state3, config3, 
              is_log_transformed=is_log_transformed)
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
            rslt_i <- f.tune2(state3, config3, 
              is_log_transformed=is_log_transformed)
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
              rslt_i <- f.tune2(state3, config3, 
                is_log_transformed=is_log_transformed)
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
            rslt_i <- f.tune2(state3, config3, 
              is_log_transformed=is_log_transformed)
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
            rslt_i <- f.tune2(state3, config3, 
              is_log_transformed=is_log_transformed)
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
            rslt_i <- f.tune2(state3, config3, 
              is_log_transformed=is_log_transformed)
            rslt <- rbind(rslt, rslt_i)
            f.log_obj(rslt, config=config3)
          }
        } else if(impute_method %in% c("sample_lod", "glm_binom", "glmnet", 
            "rf", "missforest", "none")) {
          
          f.log_block("normalization_method:", normalization_method, 
            "; impute_method:", impute_method, "; test_method:", test_method, 
            config=config3)
          if(any(is.na(c(state3$expression))) && test_method %in% c("msqrob", "voom")) {
            f.msg("WARNING: OIL_WATER: skipping test_method", test_method, 
              "because of NAs in expression", config=config3)
            next  ## next impute_method in impute_methods
          }          
          f.log_block("filter, impute, and test", config=config3)
          rslt_i <- f.tune2(state3, config3, 
            is_log_transformed=is_log_transformed)
          rslt <- rbind(rslt, rslt_i)
          f.log_obj(rslt, config=config3)
        } else {
          f.err("tune: unexpected impute_method:", 
            impute_method, config=config3)
        } ## if impute_method %in% ...
      }   ## for impute_method in impute_methods
    }     ## for test_method in test_methods
  }       ## for normalization_method in normalization_methods
  
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
#' @param dir_in Character scalar with path to directory containing tuning results.
#' @param prefix Character scalar with prefix (if any) of tuning result filenames.
#' @param suffix Character scalar with distinctive suffix (required) of tuning results filenames.
#' @param config List with at least \code{log_file} defined (can be \code{""}).
#' @param fdr_cutoff Numeric scalar between \code{0} and \code{1.0} specifying 
#'   cutoff for false discovery rate. Trials not meeting cutoff are moved to the 
#'   bottom of the output \code{data.frame}.
#' @return A \code{data.frame} with the following columns:
#'   \tabular{ll}{
#'     \code{nhits}      \cr \tab Number of significant hits. \cr
#'     \code{fdr}        \cr \tab False discovery rate. \cr
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

  pat <- paste0(gsub("(\\W)", "\\\\\\1", suffix), "$")
  files <- sort(list.files(path=dir_in, pattern=pat))
  
  unperm_file <- paste0(prefix, "0", suffix)
  i0 <- files %in% unperm_file
  if(sum(i0) != 1) f.err("tune_check: no unperm file found; looking for:", 
    unperm_file, config=config)
  perm_files <- files[!i0]
  f.msg("found 1 unpermuted file and", length(perm_files), "permuted files\n", 
    config=config)
  dat0 <- utils::read.table(paste(dir_in, unperm_file, sep="/"), header=T, 
    sep="\t", quote="", as.is=T)

  dat0$nhits[is.na(dat0$nhits)] <- 0
  dat0$ntests[is.na(dat0$ntests)] <- 0

  obj <- list()
  for(perm_file in perm_files) {
    f.msg("reading", perm_file, config=config)
    prfx <- sub(pat, "", perm_file)
    dat_i <- utils::read.table(paste(dir_in, perm_file, sep="/"), header=T, 
      sep="\t", quote="", as.is=T)
    dat_i$perm_prfx <- prfx
    obj[[perm_file]] <- dat_i
  }

  dat1 <- do.call(rbind, obj)
  rownames(dat1) <- NULL

  dat1$nhits[is.na(dat1$nhits)] <- 0
  dat1$ntests[is.na(dat1$ntests)] <- 0

  k0 <- apply(dat0[, 1:9], 1, paste, collapse=":")
  k1 <- apply(dat1[, 1:9], 1, paste, collapse=":")
  
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
  dat0$fdr <- dat0$max / nhits   ## used to be $avg
  dat0$fdr[dat0$fdr > 1] <- 1.0

  dat0 <- dat0[, c("nhits", "fdr", "max1", "mid1", "avg1", "sd1", "norm", "nquant", 
    "impute", "iquant", "scale", "span", "npcs", "k", "test")]

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

