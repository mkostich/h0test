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
  tbl <- f.test(out$state, out$config)
    
  rslt <- data.frame(norm=config$norm_method, norm_quant=config$norm_quantile, 
    impute=config$impute_method, imp_quant=config$impute_quantile, 
    scale=config$impute_scale, test=config$test_method, 
    perm=config$permute_var, nhits=sum(tbl$adj.P.Val < 0.05), ntests=nrow(tbl), 
    time=format(Sys.time(), "%H:%M:%S"), stringsAsFactors=F)
  
  return(rslt)
}

#' Basic tuning loop
#' @description
#'   Run a basic tuning loop according.
#' @details
#'   Run a basic tuning loop: 
#'     \code{f.normalize() -> f.combine_reps() -> f.filter() -> f.impute()} 
#' @param config List with configuration values like those returned by \code{f.new_config()}.
#' @param permute_vars Character vector of variables to permute when generating null distribution. 
#'   Should be in \code{config$test_term}. Do not include \code{""}.
#' @param norm_methods Character vector of methods to try. One or more of:
#'   \code{c("TMM", "TMMwsp", "RLE", "upperquartile", "q50", "q75", "cpm", "vsn", "qquantile", "log2", "none")}.
#' @param impute_methods Character vector of methods to try. One or more of:
#'   \code{c("sample_lod", "unif_sample_lod", "unif_global_lod", "rnorm_feature", "glm_binom", "loess_logit", "glmnet", "rf", "none")}.
#' @param impute_quantiles Numeric vector of quantiles to try for \code{unif_} methods. 
#'   One or more values between \code{0.0} and \code{1.0}.
#' @param impute_scales Numeric vector of scales to try for \code{rnorm_feature}. 
#'   See \code{f.impute_rnorm_feature()} \code{scale}. parameter.
#' @param test_methods Character vector with one or more of: \code{c("trend", "voom")}.
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
#'   \dontrun{
#'     config <- h0testr::f.new_config()
#'     config$frm=~age+sex+age:sex
#'     config$test_term="age:sex"
#'     config$test_method="trend"
#'     tbl <- f.tune(config, permute_vars="age")
#'     head(tbl)
#'   }

f.tune <- function(
    config, 
    permute_vars, 
    norm_methods=c("TMM", "TMMwsp", "RLE", "upperquartile", "q50", "q75", 
      "cpm", "vsn", "qquantile", "log2", "none"),
    impute_methods=c("sample_lod", "unif_sample_lod", "unif_global_lod", 
      "rnorm_feature", "glm_binom", "loess_logit", "glmnet", "rf", "none"),
    impute_quantiles=c(0, 0.01, 0.05, 0.1, 0.25), 
    impute_scales=c(1, 0.5, 0.25, 0.1),
    test_methods=c("voom", "trend")) {

  for(permute_var in c("", permute_vars)) {
    config$permute_var <- permute_var
    
    ## load data:
    f.log_block("loading data", config=config)
    out <- f.load_data(config)            ## overwritten at each iteration
    state1 <- out$state                   ## save for subsequent iterations
    config1 <- out$config                 ## save for subsequent iterations
    rslt <- NULL
    
    for(norm_method in norm_methods) {
      config1$norm_method <- norm_method
      f.log_block("normalize and combine", config=config1)
      out <- f.tune1(state1, config1, norm_method=norm_method)
      state2 <- out$state                 ## save for subsequent iterations
      config2 <- out$config               ## save for subsequent iterations
      
      for(test_method in test_methods) {
        config2$test_method <- test_method
        for(impute_method in impute_methods) {
          config2$impute_method <- impute_method
          if(impute_method %in% c("unif_global_lod", "unif_sample_lod")) {
            for(impute_quantile in impute_quantiles) {
              config2$impute_quantile <- impute_quantile
              f.log_block("impute and test", config=config2)
              rslt_i <- f.tune2(state2, config2)
              rslt <- rbind(rslt, rslt_i)
              f.log_obj(rslt, config=config2)
            }
          } else if(impute_method %in% c("rnorm_feature")) {
            for(impute_scale in impute_scales) {
              config2$impute_scale <- impute_scale
              f.log_block("impute and test", config=config2)
              rslt_i <- f.tune2(state2, config2)
              rslt <- rbind(rslt, rslt_i)
              f.log_obj(rslt, config=config2)
            }
          } else if(impute_method %in% c("sample_lod")) {
            f.log_block("impute and test", config=config2)
            rslt_i <- f.tune2(state2, config2)
            rslt <- rbind(rslt, rslt_i)
            f.log_obj(rslt, config=config2)
          }
        }
      }
    }
  }
  f.log_block("returning result", config=config)
  return(rslt)
}
