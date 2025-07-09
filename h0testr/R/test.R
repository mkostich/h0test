#' Hypothesis testing using \code{limma::voom}
#' @description
#'   Tests for differential expression using the \code{limma::voom()} function.
#' @details
#'   Tests for differential expression using the \code{limma::voom()}. 
#'     Model is fit to \code{config$frm} and an F-test is performed for 
#'     whether the effect of \code{config$test_term} on 
#'     \code{state$expression} is zero. 
#' @param state List with elements formatted like the list returned by \code{f.read_data()}:
#'   \tabular{ll}{
#'     \code{expression} \cr \tab Numeric matrix with non-negative expression values. \cr
#'     \code{features}   \cr \tab A data.frame with feature meta-data for rows of expression. \cr
#'     \code{samples}    \cr \tab A data.frame with observation meta-data for columns of expression. \cr
#'   } 
#' @param config List with configuration values. Uses the following keys:
#'   \tabular{ll}{
#'     \code{feat_col}       \cr \tab Name of column in \code{feature_file_in} that corresponds to rows of \code{data_file_in}. \cr
#'     \code{obs_col}        \cr \tab Name of column in \code{sample_file_in} that corresponds to columns of \code{data_file_in}. \cr
#'     \code{frm}            \cr \tab Formula (formula) to be fit \cr
#'     \code{test_term}      \cr \tab Term (character) to be tested for non-zero coefficient. \cr
#'     \code{sample_factors} \cr \tab List specifying levels of factor variables in \code{config$frm} (see examples). \cr
#'   }
#' @param normalize.method Character in 
#'   \code{c("TMM", "TMMwsp", "RLE", "upperquartile", "none")}.
#' @return A data.frame containing results of test.
#' @examples
#' set.seed(101)
#' ## no missing values: mnar_c0=-Inf, mnar_c1=0, mcar_p=0
#' exprs <- h0testr::f.sim2(n_samps1=6, n_samps2=6, n_genes=25, 
#'   n_genes_signif=5, fold_change=2, mnar_c0=-Inf, mnar_c1=0, mcar_p=0)
#' exprs <- log2(exprs + 1)
#' feats <- data.frame(feature_id=rownames(exprs))
#' samps <- data.frame(observation_id=colnames(exprs), 
#'   condition=c(rep("placebo", 6), rep("drug", 6)))
#' state <- list(expression=exprs, features=feats, samples=samps)
#' 
#' config <- h0testr::f.new_config()
#' config$feat_col <- config$feat_id_col <- config$gene_id_col <- "feature_id"
#' config$obs_col <- config$obs_id_col <- config$sample_id_col <- "observation_id"
#' config$frm <- ~condition
#' config$test_term <- "condition"
#' config$sample_factors <- list(condition=c("placebo", "drug"))
#'
#' ## set up and check configuration, including covariates:
#' out <- h0testr::f.initialize(state, config)
#' 
#' tbl <- h0testr::f.test_voom(out$state, out$config)
#' print(tbl)

f.test_voom <- function(state, config, normalize.method="none") {

  f.check_config(config)
  f.check_state(state, config)
  
  exprs <- state$expression
  i <- apply(exprs, 1, function(v) any(is.na(v)))
  exprs <- exprs[!i, ]
  
  design <- stats::model.matrix(config$frm, data=state$samples)
  obj <- limma::voom(exprs, design, plot=F, normalize.method=normalize.method)
  fit <- limma::lmFit(obj, design)
  fit <- limma::eBayes(fit, trend=F)
  
  frm0 <- stats::as.formula(paste("~0 + ", config$test_term))
  cols <- colnames(stats::model.matrix(frm0, data=state$samples))
  cols <- cols[cols %in% colnames(design)]
  i <- colnames(design) %in% cols
  tbl <- limma::topTable(fit, coef=which(i), number=Inf)
  
  f.msg("tested", nrow(exprs), "genes", config=config)
  f.msg("found", sum(tbl$adj.P.Val < 0.05, na.rm=T), "hits", config=config)
  
  return(tbl)
}

#' Hypothesis testing using \code{limma} trend.
#' @description
#'   Test for differential expression using the \code{limma} 
#'   package trend approach.
#' @details
#'   Tests for differential expression using the \code{limma} pakcage by running 
#'     \code{lmFit()}, then \code{eBayes()}, then \code{topTable()}. Model is fit 
#'     to \code{config$frm} and an F-test is performed for whether the effect of 
#'     \code{config$test_term} on \code{state$expression} is zero. 
#' @param state List with elements formatted like the list returned by \code{f.read_data()}:
#'   \tabular{ll}{
#'     \code{expression} \cr \tab Numeric matrix with non-negative expression values. \cr
#'     \code{features}   \cr \tab A data.frame with feature meta-data for rows of expression. \cr
#'     \code{samples}    \cr \tab A data.frame with observation meta-data for columns of expression. \cr
#'   } 
#' @param config List with configuration values. Uses the following keys:
#'   \tabular{ll}{
#'     \code{feat_col}       \cr \tab Name of column in \code{feature_file_in} corresponding to \code{rownames(state$expression)}. \cr
#'     \code{obs_col}        \cr \tab Name of column in \code{sample_file_in} corresponding to \code{colnames(state$expression)}. \cr
#'     \code{frm}            \cr \tab Formula (formula) to be fit \cr
#'     \code{test_term}      \cr \tab Term (scalar character) to be tested for non-zero coefficient. \cr
#'     \code{sample_factors} \cr \tab List specifying levels of factor variables in \code{config$frm} (see examples). \cr
#'   }
#' @return A \code{data.frame} containing results of test.
#' @examples
#' set.seed(101)
#' ## no missing values: mnar_c0=-Inf, mnar_c1=0, mcar_p=0
#' exprs <- h0testr::f.sim2(n_samps1=6, n_samps2=6, n_genes=25, 
#'   n_genes_signif=5, fold_change=2, mnar_c0=-Inf, mnar_c1=0, mcar_p=0)
#' exprs <- log2(exprs + 1)
#' feats <- data.frame(feature_id=rownames(exprs))
#' samps <- data.frame(observation_id=colnames(exprs), 
#'   condition=c(rep("placebo", 6), rep("drug", 6)))
#' state <- list(expression=exprs, features=feats, samples=samps)
#' 
#' config <- h0testr::f.new_config()
#' config$feat_col <- config$feat_id_col <- config$gene_id_col <- "feature_id"
#' config$obs_col <- config$obs_id_col <- config$sample_id_col <- "observation_id"
#' config$frm <- ~condition
#' config$test_term <- "condition"
#' config$sample_factors <- list(condition=c("placebo", "drug"))
#'
#' ## set up and check covariates and parameters:
#' out <- h0testr::f.initialize(state, config)
#' 
#' tbl <- h0testr::f.test_trend(state, config)
#' print(tbl)

f.test_trend <- function(state, config) {

  f.check_config(config)
  f.check_state(state, config)

  if(!is.matrix(state$expression)) {
    f.err("f.test_trend: !is.matrix(state$expression)", config=config)
  }
  
  design <- stats::model.matrix(config$frm, data=state$samples)
  fit <- limma::lmFit(state$expression, design)
  fit <- limma::eBayes(fit, trend=T)
  
  frm0 <- stats::as.formula(paste("~0 + ", config$test_term))
  cols <- colnames(stats::model.matrix(frm0, data=state$samples))
  cols <- cols[cols %in% colnames(design)]
  i <- colnames(design) %in% cols
  tbl <- limma::topTable(fit, coef=which(i), number=Inf)
  
  f.msg("tested", nrow(state$expression), "genes", config=config)
  f.msg("found", sum(tbl$adj.P.Val < 0.05, na.rm=T), "hits", config=config)
  
  return(tbl)
}

#' Hypothesis testing
#' @description
#'   Test hypotheses using .
#' @details
#'   Tests for differential expression using method specified in config. 
#' @param state List with elements formatted like the list returned by `f.read_data()`:
#'   \tabular{ll}{
#'     \code{expression} \cr \tab Numeric matrix with non-negative expression values. \cr
#'     \code{features}   \cr \tab A data.frame with feature meta-data for rows of expression. \cr
#'     \code{samples}    \cr \tab A data.frame with observation meta-data for columns of expression. \cr
#'   } 
#' @param config List with configuration values. Uses the following keys:
#'   \tabular{ll}{
#'     \code{feat_col}       \cr \tab Name of column in \code{state$fetaures} matching \code{rownames(state$expression)}. \cr
#'     \code{obs_col}        \cr \tab Name of column in \code{state$samples} matching \code{colnames(state$expression)}. \cr
#'     \code{frm}            \cr \tab Formula (formula) to be fit \cr
#'     \code{test_term}      \cr \tab Term (character scalar) to be tested for non-zero coefficient. \cr
#'     \code{sample_factors} \cr \tab List specifying levels of factor variables in \code{config$frm} (see examples). \cr
#'     \code{test_method}    \cr \tab Character scalar in \code{c("trend", "voom")}. \cr
#'   }
#' @return A data.frame containing results of test.
#' @examples
#' set.seed(101)
#' ## no missing values: mnar_c0=-Inf, mnar_c1=0, mcar_p=0
#' exprs <- h0testr::f.sim2(n_samps1=6, n_samps2=6, n_genes=25, 
#'   n_genes_signif=5, fold_change=2, mnar_c0=-Inf, mnar_c1=0, mcar_p=0)
#' exprs <- log2(exprs + 1)
#' feats <- data.frame(feature_id=rownames(exprs))
#' samps <- data.frame(observation_id=colnames(exprs), 
#'   condition=c(rep("placebo", 6), rep("drug", 6)))
#' state <- list(expression=exprs, features=feats, samples=samps)
#' 
#' config <- h0testr::f.new_config()
#' config$feat_id_col <- config$gene_id_col <- config$feat_col <- "feature_id"
#' config$obs_id_col <- config$sample_id_col <- config$obs_col <- "observation_id"
#' config$frm <- ~condition
#' config$test_term <- "condition"
#' config$sample_factors <- list(condition=c("placebo", "drug"))
#' config$test_method <- "trend"
#' config$save_state <- FALSE
#' 
#' ## set up and check covariates and parameters:
#' out <- h0testr::f.initialize(state, config)
#' 
#' out <- h0testr::f.test(state, config)
#' print(out)

f.test <- function(state, config) {

  f.report_config(config)

  if(config$test_method %in% "voom") {
    tbl <- f.test_voom(state, config)
  } else if(config$test_method %in% "trend") {
    tbl <- f.test_trend(state, config)
  } else if(config$test_method %in% "none") {
    f.msg("skipping testing: TEST_METHOD %in% 'none'", config=config)
    return(NULL)
  } else f.err("f.test: unexpected TEST_METHOD:", config$test_method, config=config)
  
  rownames(state$features) <- state$features[[config$feat_col]]
  if(!all(rownames(tbl) %in% rownames(state$features))) {
    f.err("f.test: !all(rownames(tbl) %in% rownames(state$features))", config=config)
  }
  
  tbl0 <- state$features[rownames(tbl), , drop=F]
  tbl <- cbind(tbl0, tbl)
  rownames(tbl) <- NULL
  
  if(config$save_state) {
    file_out <- paste0(config$dir_out, "/", length(config$run_order) + 3, 
      config$result_mid_out, config$suffix_out)
    f.log("writing results to", file_out, config=config)
    f.save_tsv(tbl, file_out)
  }
  
  return(tbl)
}

