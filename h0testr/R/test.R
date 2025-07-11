#' Hypothesis testing using the \code{DEqMS} package
#' @description
#'   Tests for differential expression using the \code{DEqMS::spectraCounteBayes()} function.
#' @details
#'   The \code{DEqMS::spectraCounteBayes()} model is fit to \code{config$frm} 
#'     and a moderated t-test is performed for whether the effect of 
#'     \code{config$test_term} on \code{state$expression} is zero. 
#'   Flow is:
#'     \tabular{l}{
#'       1. for each gene, count number of associated peptides. \cr
#'       2. Fit linear model to \code{config$frm} using \code{limma::lmFit()}. \cr
#'       3. Calculate statistics using \code{limma::eBayes()} on fitted model. \cr
#'       4. Append peptide counts to model returned by \code{limma::eBayes()}. \cr
#'       5. Adjust statistics using \code{DEqMS::spectraCounteBayes()}. \cr
#'       6. Generate hit table with \code{DEqMS::outputResult()}. \cr
#'     }
#' @param state List with elements like those returned by \code{f.read_data()}:
#'   \tabular{ll}{
#'     \code{expression} \cr \tab Numeric matrix with non-negative expression values. \cr
#'     \code{features}   \cr \tab A data.frame with feature meta-data for rows of expression. \cr
#'     \code{samples}    \cr \tab A data.frame with observation meta-data for columns of expression. \cr
#'   } 
#' @param config List with configuration values. Uses the following keys:
#'   \tabular{ll}{
#'     \code{gene_id_col}    \cr \tab Name of column in \code{state$features} with unique gene/protein group ids. \cr
#'     \code{feat_col}       \cr \tab Name of column in \code{state$features} corresponding to \code{rownames(state$expression)}. \cr
#'     \code{obs_col}        \cr \tab Name of column in \code{state$samples} corresponding to \code{colnames(state$expression)}. \cr
#'     \code{frm}            \cr \tab Formula (formula) to be fit \cr
#'     \code{test_term}      \cr \tab Term (character) to be tested for non-zero coefficient. \cr
#'   }
#' @param trend Logical scalar. Whether \code{limma::eBayes()} should use trended dispersion estimate.
#' @return
#'   A data.frame with results of test. Columns include: 
#'     \code{c("logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B", "gene", "count", "sca.t", "sca.P.Value", "sca.adj.pval")}.
#'   Initial statistics from \code{limma}. Columns beginning with \code{sca.} 
#'     added by \code{DEqMS}.
#' @examples
#' ## lengthy setup of expression data:
#' set.seed(101)
#' nsamps <- 6
#' sim <- h0testr::f.sim2(
#'   n_samps1=nsamps, n_samps2=nsamps, n_genes=100, n_genes_signif=20, 
#'   fold_change=1, peps_per_gene=10, reps_per_sample=1, 
#'   p_drop=0.33, mnar_c0=-Inf, mnar_c1=0, mcar_p=0
#' )
#' exprs <- sim$mat
#' gene <- strsplit(rownames(exprs), "_")
#' gene <- sapply(gene, function(v) unlist(v)[1])
#' feats <- data.frame(pep=rownames(exprs), gene=gene)
#' samps <- data.frame(
#'   obs=colnames(exprs), 
#'   grp=c(rep("ctl", nsamps), rep("trt", nsamps)),
#'   sex=rep(c("M", "F"), round(ncol(exprs) / 2))
#' )
#' state <- list(expression=exprs, features=feats, samples=samps)
#' rm(sim, exprs, gene, feats, samps)
#'
#' ## setup config and prep variables of interest for testing:
#' config <- list(
#'   obs_id_col="obs",
#'   sample_id_col="obs",
#'   feat_id_col="pep",
#'   gene_id_col="gene",
#'   frm=~grp+sex+grp:sex, 
#'   test_term="grp",
#'   sample_factors=list( 
#'     grp=c("ctl", "trt"), 
#'     sex=c("F", "M")
#'   )
#' )
#' out <- h0testr::f.initialize(state, config, minimal=TRUE)
#' 
#' ## actual test:
#' tbl <- h0testr::f.test_deqms(out$state, out$config)
#' head(tbl)

f.test_deqms <- function(state, config, trend=FALSE) {
  
  f.check_config(config)
  f.check_state(state, config)
  
  tbl <- table(state$features[, config$gene_id_col], useNA="ifany")
  tbl <- as.data.frame(tbl)
  counts <- tbl$Freq
  names(counts) <- tbl$Var1

  out <- h0testr::f.combine_peps(state, config, method="medianPolish", rescale=TRUE)
  state <- out$state
  config <- out$config
  design <- stats::model.matrix(config$frm, data=state$samples)
  cols_des <- colnames(design)
  
  frm0 <- stats::as.formula(paste("~0 + ", config$test_term))
  cols_pick <- colnames(stats::model.matrix(frm0, data=state$samples))
  cols_pick <- cols_pick[cols_pick %in% cols_des]
  idx <- which(cols_des %in% cols_pick)
  if(length(idx) != 1) f.err("length(idx) != 1", config=config)

  fit <- limma::lmFit(state$expression, design)
  fit <- limma::eBayes(fit, trend=trend)
  fit$count <- counts[rownames(fit$coefficients)]
  fit <- DEqMS::spectraCounteBayes(fit, fit.method="loess")
  tbl <- DEqMS::outputResult(fit, coef_col=idx)

  return(tbl)
}

#' Hypothesis testing using the \code{msqrob2} package
#' @description
#'   Tests for differential expression using the \code{msqrob2::msqrob()} function.
#' @details
#'   The \code{DEqMS::spectraCounteBayes()} model is fit to \code{config$frm} 
#'     and a moderated t-test is performed for whether the effect of 
#'     \code{config$test_term} on \code{state$expression} is zero. 
#'   Flow is:
#'     \tabular{l}{
#'       1. Create a \code{QFeatures} object with \code{QFeatures::readQFeatures()}. 
#'       2. Convert \code{NA}s to zero using \code{QFeatures::zeroIsNA()}. 
#'       3. Make gene/protein-group expression with \code{QFeatures::aggregateFeatures()}.
#'       4. Estimate model parameters using \code{msqrob2::msqrob()}.
#'       5. Calculate test statistics with \code{msqrob2::hypothesisTest()}.
#'     }
#' @param state List with elements like those returned by \code{f.read_data()}:
#'   \tabular{ll}{
#'     \code{expression} \cr \tab Numeric matrix with non-negative expression values. \cr
#'     \code{features}   \cr \tab A data.frame with feature meta-data for rows of expression. \cr
#'     \code{samples}    \cr \tab A data.frame with observation meta-data for columns of expression. \cr
#'   } 
#' @param config List with configuration values. Uses the following keys:
#'   \tabular{ll}{
#'     \code{gene_id_col}    \cr \tab Name of column in \code{state$features} with unique gene/protein group ids. \cr
#'     \code{feat_col}       \cr \tab Name of column in \code{state$features} corresponding to \code{rownames(state$expression)}. \cr
#'     \code{obs_col}        \cr \tab Name of column in \code{state$samples} corresponding to \code{colnames(state$expression)}. \cr
#'     \code{frm}            \cr \tab Formula (formula) to be fit \cr
#'     \code{test_term}      \cr \tab Term (character) to be tested for non-zero coefficient. \cr
#'   }
#' @param maxit Integer scalar >= 1. How many iterations to use for \code{rlm} fitting.
#' @return
#'   A data.frame with test results. Columns include \code{config$gene_id_col} and: 
#'     \code{c("nNonZero .n", "logFC", "se", "df", "t", "pval", "adjPval")}.
#' @examples
#' ## lengthy setup of expression data:
#' set.seed(101)
#' nsamps <- 6
#' sim <- h0testr::f.sim2(
#'   n_samps1=nsamps, n_samps2=nsamps, n_genes=100, n_genes_signif=20, 
#'   fold_change=1, peps_per_gene=10, reps_per_sample=1, 
#'   p_drop=0.33, mnar_c0=-Inf, mnar_c1=0, mcar_p=0
#' )
#' exprs <- sim$mat
#' gene <- strsplit(rownames(exprs), "_")
#' gene <- sapply(gene, function(v) unlist(v)[1])
#' feats <- data.frame(pep=rownames(exprs), gene=gene)
#' samps <- data.frame(
#'   obs=colnames(exprs), 
#'   grp=c(rep("ctl", nsamps), rep("trt", nsamps)),
#'   sex=rep(c("M", "F"), round(ncol(exprs) / 2))
#' )
#' state <- list(expression=exprs, features=feats, samples=samps)
#' rm(sim, exprs, gene, feats, samps)
#'
#' ## setup config and prep variables of interest for testing:
#' config <- list(
#'   obs_id_col="obs",
#'   sample_id_col="obs",
#'   feat_id_col="pep",
#'   gene_id_col="gene",
#'   frm=~grp+sex+grp:sex, 
#'   test_term="grp",
#'   sample_factors=list( 
#'     grp=c("ctl", "trt"), 
#'     sex=c("F", "M")
#'   )
#' )
#' out <- h0testr::f.initialize(state, config, minimal=TRUE)
#' 
#' ## actual test:
#' tbl <- h0testr::f.test_msqrob(out$state, out$config)
#' head(tbl)

f.test_msqrob <- function(state, config, maxit=20) {

  f.check_config(config)
  f.check_state(state, config)

  exprs <- as.data.frame(state$expression)
  n_non0 <- apply(exprs, 1, function(v) sum(v > 0, na.rm=T))
  exprs <- cbind(fnames=rownames(exprs), exprs)
  rownames(exprs) <- NULL
  obj <- QFeatures::readQFeatures(table=exprs, ecol=2:ncol(exprs), 
    fnames="fnames", name="features")  

  ## SummarizedExperiment::rowData(obj[["features"]]) <- S4Vectors::DataFrame(state$features)
  feats <- state$features
  for(nom in names(feats)) {
    if(is.factor(feats[[nom]])) {
      SummarizedExperiment::rowData(obj[["features"]])[[nom]] <- as.character(feats[[nom]])
    } else {
      SummarizedExperiment::rowData(obj[["features"]])[[nom]] <- feats[[nom]]
    }
  }
  SummarizedExperiment::rowData(obj[["features"]])$nNonZero <- n_non0

  samps <- state$samples
  for(nom in names(samps)) {
    if(is.factor(samps[[nom]])) {
      SummarizedExperiment::colData(obj)[[nom]] <- as.character(samps[[nom]])
    } else {
      SummarizedExperiment::colData(obj)[[nom]] <- samps[[nom]]
    }
  }
  
  ## convert NA to 0:
  obj <- QFeatures::zeroIsNA(obj, i="features")
  
  ## aggregate peptides into genes:
  obj <- QFeatures::aggregateFeatures(obj, i="features", 
    fcol=config$gene_id_col, na.rm=T, name="genes", fun=base::colMeans)
  
  obj <- msqrob2::msqrob(object=obj, i="genes", formula=config$frm, maxitRob=maxit)
  
  cols_des <- colnames(stats::model.matrix(config$frm, data=state$samples))
  frm0 <- stats::as.formula(paste("~0 + ", config$test_term))
  cols_pick <- colnames(stats::model.matrix(frm0, data=state$samples))
  cols_pick <- cols_pick[cols_pick %in% cols_des]
  idx <- which(cols_des %in% cols_pick)
  if(length(idx) != 1) f.err("length(idx) != 1", config=config)
  cols_pick <- cols_des[idx]
  
  con <- msqrob2::makeContrast(contrasts=paste0(cols_pick, "=0"), 
    parameterNames=c(cols_pick))
  
  obj <- msqrob2::hypothesisTest(object=obj, i="genes", contrast=con, 
    modelColumn="msqrobModels")
  
  dat <- SummarizedExperiment::rowData(obj[["genes"]])    
  tbl <- dat[[cols_pick]]
  dat <- dat[, c(config$gene_id_col, "nNonZero", ".n")]
  dat <- as.data.frame(dat)
  tbl <- cbind(dat[rownames(tbl), ], tbl)
  tbl <- tbl[order(tbl$adjPval, -abs(tbl$logFC)), ]
  rownames(tbl) <- NULL
  
  return(tbl)
}

#' Hypothesis testing using \code{limma::voom}
#' @description
#'   Tests for differential expression using the \code{limma::voom()} function.
#' @details
#'   The \code{limma::voom()} model is fit to \code{config$frm} and an F-test 
#'     is performed for whether the effect of \code{config$test_term} on 
#'     \code{state$expression} is zero. 
#' @param state List with elements like those returned by \code{f.read_data()}:
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
#'   n_genes_signif=5, fold_change=2, mnar_c0=-Inf, mnar_c1=0, mcar_p=0)$mat
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
#'   package \code{trend} approach.
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
#'   n_genes_signif=5, fold_change=2, mnar_c0=-Inf, mnar_c1=0, mcar_p=0)$mat
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
#'   See invididual \code{f.test_*} methods for more details.
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
#'   n_genes_signif=5, fold_change=2, mnar_c0=-Inf, mnar_c1=0, mcar_p=0)$mat
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

