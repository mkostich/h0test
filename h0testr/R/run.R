#' Run a basic workflow
#' @description
#'   Run a basic workflow according to: \code{config$run_order}.
#' @details
#'   Run a basic workflow: \code{f.load_data -> config$run_order -> f.test}, 
#'     where \code{config$run_order} is vector of functions which are run in 
#'     the specified order.
#'   See documentation for \code{h0testr::f.new_config()} 
#'     for more detailed description of configuration parameters. 
#' @param config List with configuration values like those returned by 
#'   \code{f.new_config()}.
#' @return A list with the following elements:
#'   \tabular{ll}{
#'     \code{state}  \cr \tab A list with elements \code{$expression}, \code{$features}, and \code{$samples}. \cr
#'     \code{config} \cr \tab A list with configuration settings. \cr
#'     \code{tbl}    \cr \tab A data.frame containing results of test.. \cr
#'   }
#' @examples
#' config <- h0testr::f.new_config()
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
#' config$run_order <- c("normalize", "combine_reps", "filter", "impute")
#' config$save_state <- FALSE
#'
#' print(config$run_order)
#'
#' out <- h0testr::f.run(config)   ## run workflow
#' head(out$original)              ## hit table as returned by underlying software
#' head(out$standard)              ## hit table in a standardized format

f.run <- function(config) {

  f.report_config(config)

  f.log_block("starting f.load_data", config=config)
  out <- f.load_data(config)
  
  for(f in config$run_order) {
    f <- paste0("f.", f)
    f.log_block("starting", f, config=config)
    f <- get(f)
    out <- f(out$state, out$config)
  }
  
  f.log_block("starting f.test", config=out$config)
  tbl_list <- f.test(out$state, out$config)
  
  return(list(state=out$state, config=out$config, 
    original=tbl_list$original, standard=tbl_list$standard))
}


## helper for f.tune():

f.tune1 <- function(state, config, norm_method) {

    if(norm_method %in% "q50") {
      config$norm_method <- "quantile"
      config$norm_quantile <- 0.5
    } else if(norm_method %in% c("q75")) {
      config$norm_method <- "quantile"
      config$norm_quantile <- 0.75
    } else if(norm_method %in% "upperquartile") {
      config$norm_quantile <- 0.75
    }
    
    ## normalize and combine reps:
    out <- f.normalize(state, config)
    out <- f.combine_reps(out$state, out$config)
    
    return(out)
}

## helper for f.tune():

f.tune2 <- function(state, config) {
  
  out <- f.filter(state, config)
  out <- f.impute(out$state, out$config)
  tbl_list <- f.test(out$state, out$config)
  tbl <- tbl_list$standard
    
  rslt <- data.frame(norm=config$norm_method, norm_quant=config$norm_quantile, 
    impute=config$impute_method, imp_quant=config$impute_quantile, 
    scale=config$impute_scale, span=config$impute_span, test=config$test_method, 
    perm=config$permute_var, nhits=sum(tbl$adj.P.Val < 0.05), ntests=nrow(tbl), 
    time=format(Sys.time(), "%H:%M:%S"), stringsAsFactors=F)
  
  return(rslt)
}

#' Basic tuning loop
#' @description
#'   Run a basic tuning loop according.
#' @details
#'   Run a basic tuning loop that iterates over parameter combinations, 
#'     performing: 
#'     \code{f.normalize() -> f.combine_reps() -> f.filter() -> f.impute()}. 
#'     Do one run with \code{config$permute_var=""}, and \code{N} runs (we 
#'     recommend \code{N >= 20} with \code{config$permute_var} set to the 
#'     name of a variable in \code{config$test_term}.
#'   See documentation for \code{h0testr::f.new_config()} 
#'     for more detailed description of configuration parameters. 
#' @param config List with configuration values like those returned by \code{f.new_config()}.
#' @param norm_methods Character vector of methods to try. One or more of:
#'   \code{c("TMM", "TMMwsp", "RLE", "upperquartile", "q50", "q75", "cpm", "vsn", "qquantile", "log2", "none")}.
#' @param impute_methods Character vector of methods to try. One or more of:
#'   \code{c("sample_lod", "unif_sample_lod", "unif_global_lod", "rnorm_feature", "glm_binom", "loess_logit", "glmnet", "rf", "none")}.
#' @param impute_quantiles Numeric vector of quantiles to try for \code{f.impute_unif_*} methods. 
#'   One or more values between \code{0.0} and \code{1.0}.
#' @param impute_scales Numeric vector of scales to try for \code{f.impute_rnorm_feature}. 
#'   See \code{f.impute_rnorm_feature()} \code{scale}. parameter.
#' @param impute_spans Numeric vector of spans to try for \code{f.impute_loess_logit}.
#' @param test_methods Character vector with one or more of: 
#'   \code{c("voom", "trend", "deqms", "msqrob", "proda")}.
#' @return A data.frame with the following columns:
#'   \tabular{ll}{
#'     \code{norm}       \cr \tab Normalization method (character). \cr
#'     \code{norm_quant} \cr \tab Normalization quantile (numeric). \cr
#'     \code{impute}     \cr \tab Imputation method (character). \cr
#'     \code{imp_quant}  \cr \tab Imputation quantile (numeric). \cr
#'     \code{scale}      \cr \tab Imputation scale for \code{f.impute_rnorm_feature}. \cr
#'     \code{test}       \cr \tab Test method (character). \cr
#'     \code{perm}       \cr \tab Permuted variable (character). \cr
#'     \code{nhits}      \cr \tab Number of hits (numeric). \cr
#'     \code{ntests}     \cr \tab Number of tests (numeric). \cr
#'     \code{time}       \cr \tab Timestamp. \cr
#'   }
#' @examples
#' config <- h0testr::f.new_config()
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
#' config$n_features_min <- 10         ## default 1000 too big for small demo dataset
#' config$save_state <- FALSE
#' 
#' ## one run with unpermuted data:
#' config$permute_var <- ""            ## no permutation
#' set.seed(101)
#' out1 <- h0testr::f.tune(config,
#'   norm_methods=c("RLE", "q75", "cpm", "log2"),
#'   impute_methods=c("sample_lod", "unif_sample_lod", "none"),
#'   impute_quantiles=c(0, 0.05, 0.1),
#'   test_methods=c("trend", "msqrob", "proda", "prolfqua")
#' )
#' ## write.table(out1, "0.condition.tune.tsv", quote=F, sep="\t", row.names=F)
#' 
#' ## one run with permuted data; run 20+ such runs w/ suffices 1:20:
#' config$permute_var <- "condition"   ## permute variable in test_term
#' set.seed(101)
#' out2 <- h0testr::f.tune(config,
#'   norm_methods=c("RLE", "q75", "cpm", "log2"),
#'   impute_methods=c("sample_lod", "unif_sample_lod", "none"),
#'   impute_quantiles=c(0, 0.05, 0.1),
#'   test_methods=c("trend", "msqrob", "proda")
#' )
#' ## write.table(out2, "1.condition.tune.tsv", quote=F, sep="\t", row.names=F)

f.tune <- function(
    config,  
    norm_methods=c("TMMwsp", "RLE", "upperquartile", "q50", "q75", 
      "max", "div.mean", "cpm", "vsn", "qquantile", "quantiles.robust",
      "log2", "none"),
    impute_methods=c("sample_lod", "unif_sample_lod", "unif_global_lod", 
      "rnorm_feature", "glm_binom", "loess_logit", "glmnet", "rf", "none"),
    impute_quantiles=c(0, 0.01, 0.05, 0.1, 0.25), 
    impute_scales=c(1, 0.5, 0.25, 0.1),
    impute_spans=c(0.25, 0.5, 0.75),
    test_methods=c("voom", "trend", "deqms", "msqrob", "proda")) {
    
  f.report_config(config)

  ## load data:
  f.log_block("loading data", config=config)
  out <- f.load_data(config)            ## overwritten at each iteration
  state1 <- out$state                   ## save for subsequent iterations
  config1 <- out$config                 ## save for subsequent iterations
  rslt <- NULL
  
  for(norm_method in norm_methods) {
    config1$norm_method <- norm_method
    f.msg("norm_method:", norm_method, config=config)
    f.log_block("normalize and combine", config=config1)
    out <- f.tune1(state1, config1, norm_method=norm_method)
    state2 <- out$state                 ## save for subsequent iterations
    config2 <- out$config               ## save for subsequent iterations
    
    for(test_method in test_methods) {
      config2$test_method <- test_method
      f.msg("test_method:", test_method, config=config2)
      for(impute_method in impute_methods) {
        config2$impute_method <- impute_method
        f.msg("impute_method:", impute_method, config=config2)
        if(impute_method %in% c("unif_global_lod", "unif_sample_lod")) {
          for(impute_quantile in impute_quantiles) {
            config2$impute_quantile <- impute_quantile
            f.msg("impute_quantile:", impute_quantile, config=config2)
            f.log_block("impute and test", config=config2)
            rslt_i <- f.tune2(state2, config2)
            rslt <- rbind(rslt, rslt_i)
            f.log_obj(rslt, config=config2)
          }
        } else if(impute_method %in% c("rnorm_feature")) {
          for(impute_scale in impute_scales) {
            config2$impute_scale <- impute_scale
            f.msg("impute_scale:", impute_scale, config=config2)
            f.log_block("impute and test", config=config2)
            rslt_i <- f.tune2(state2, config2)
            rslt <- rbind(rslt, rslt_i)
            f.log_obj(rslt, config=config2)
          }
        } else if(impute_method %in% c("loess_logit")) {
          for(impute_span in impute_spans) {
            config2$impute_span <- impute_span
            f.msg("impute_span:", impute_span, config=config2)
            f.log_block("impute and test", config=config2)
            rslt_i <- f.tune2(state2, config2)
            rslt <- rbind(rslt, rslt_i)
            f.log_obj(rslt, config=config2)
          }
        } else if(impute_method %in% c("sample_lod", "glm_binom", "glmnet", "rf")) {
          f.log_block("impute and test", config=config2)
          rslt_i <- f.tune2(state2, config2)
          rslt <- rbind(rslt, rslt_i)
          f.log_obj(rslt, config=config2)
        }
      }
    }
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
#'   See documentation for \code{h0testr::f.new_config()} 
#'     for more detailed description of configuration parameters. 
#' @param dir_in Character scalar with path to directory containing tuning results.
#' @param sfx Character scalar with distinctive suffix of tuning results files.
#' @param config List with at least \code{log_file} defined (can be \code{""}).
#' @return A data.frame with the following columns:
#'   \tabular{ll}{
#'     \code{nhits}      \cr \tab Number of significant hits (numeric). \cr
#'     \code{fdr}        \cr \tab False discovery rate (numeric). \cr
#'     \code{max0}       \cr \tab Maximum number of hits in any permutation (numeric). \cr
#'     \code{mid0}       \cr \tab Median number of hits across permutations (numeric). \cr
#'     \code{avg0}       \cr \tab Average number of hits across permutations (numeric). \cr
#'     \code{sd0}        \cr \tab Standard deviation of number of hits across permutations (numeric). \cr
#'     \code{norm}       \cr \tab Normalization method (character). \cr
#'     \code{norm_quant} \cr \tab Normalization quantile (numeric). \cr
#'     \code{impute}     \cr \tab Imputation method (character). \cr
#'     \code{imp_quant}  \cr \tab Imputation quantile (numeric). \cr
#'     \code{scale}      \cr \tab Scale (numeric). \cr
#'     \code{test}       \cr \tab Test method (character). \cr
#'   }
#' @examples
#' dir_in <- system.file("extdata/tune", package="h0testr")
#' sfx <- ".condition.tune.tsv"
#' config <- list()

#' tbl <- h0testr::f.tune_check(dir_in, sfx, config)
#' print(tbl)

f.tune_check <- function(dir_in, sfx, config) {

  pat <- paste0(gsub("(\\W)", "\\\\\\1", sfx), "$")
  files <- sort(list.files(path=dir_in, pattern=pat))
  
  unperm_file <- paste0("0", sfx)
  i0 <- files %in% unperm_file
  if(sum(i0) != 1) f.err("f.tune_check: no unperm file found; looking for:", 
    unperm_file, config=config)
  perm_files <- files[!i0]
  f.msg("found 1 unpermuted file and", length(perm_files), "permuted files\n", 
    config=config)
  dat0 <- utils::read.table(paste(dir_in, unperm_file, sep="/"), header=T, 
    sep="\t", quote="", as.is=T)
  
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
  
  k0 <- apply(dat0[, 1:6], 1, paste, collapse=":")
  k1 <- apply(dat1[, 1:6], 1, paste, collapse=":")
  
  ## permuted results in dat1; take max, median, mean, and sd of 10 permutation results:
  perm_max <- tapply(dat1$nhits, k1, max, na.rm=T)
  perm_mid <- tapply(dat1$nhits, k1, stats::median, na.rm=T)
  perm_avg <- tapply(dat1$nhits, k1, mean, na.rm=T)
  perm_sd  <- tapply(dat1$nhits, k1, stats::sd, na.rm=T)
  
  ## get them in the same order as dat1 (k1 made from dat1):
  dat0$max1 <- perm_max[k0]
  dat0$mid1 <- perm_mid[k0]
  dat0$avg1 <- perm_avg[k0]
  dat0$sd1  <- perm_sd[k0]
  dat0$perm <- NULL
  
  ## average number of false positives == average number of hits across 10 sets of permutation results;
  ##   false positive rate: (average number of false positives) / (total number of positives)
  
  dat0$nhits[dat0$nhits %in% 0] <- 1
  dat0$fdr <- dat0$avg / dat0$nhits
  
  dat0 <- dat0[
    order(dat0$nhits, -dat0$fdr, decreasing=T), 
    c("nhits", "fdr", "max1", "mid1", "avg1", "sd1", 
      "norm", "norm_quant", "impute", "imp_quant", "scale", "test")
  ]
  rownames(dat0) <- NULL
  
  return(dat0)
}
