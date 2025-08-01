#' Hypothesis testing using the \code{DEqMS} package
#' @description
#'   Tests for differential expression using the 
#'     \code{DEqMS::spectraCounteBayes()} function.
#' @details
#'   The \code{DEqMS::spectraCounteBayes()} model is fit to \code{config$frm} 
#'     and a moderated t-test is performed for whether the effect of 
#'     \code{config$test_term} on \code{state$expression} is zero. 
#'   If the number of features per gene/protein-group is the same for all
#'     features, returns same result as \code{h0testr::test_trend()}.
#'   Returns gene-level hypothesis testing results based on 
#'     peptide/precursor-level input.
#'   Flow is:
#'     \tabular{l}{
#'       1. for each gene, count number of associated peptides. \cr
#'       2. Fit linear model to \code{config$frm} using \code{limma::lmFit()}. \cr
#'       3. Calculate statistics using \code{limma::eBayes()} on fitted model. \cr
#'       4. Append peptide counts to model returned by \code{limma::eBayes()}. \cr
#'       5. Adjust statistics using \code{DEqMS::spectraCounteBayes()}. \cr
#'       6. Generate hit table with \code{DEqMS::outputResult()}. \cr
#'     }
#'   See documentation for \code{h0testr::new_config()} 
#'     for more detailed description of configuration parameters. 
#' @param state List with elements like those returned by \code{read_data()}:
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
#' sim <- h0testr::sim2(
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
#' out <- h0testr::initialize(state, config, minimal=TRUE)
#' 
#' ## actual test:
#' tbl <- h0testr::test_deqms(out$state, out$config)
#' head(tbl)

test_deqms <- function(state, config, trend=FALSE) {
  
  check_config(config)
  f.check_state(state, config)

  tbl <- table(state$features[, config$gene_id_col], useNA="ifany")
  tbl <- as.data.frame(tbl)
  counts <- tbl$Freq
  names(counts) <- tbl$Var1
  
  out <- h0testr::combine_features(state, config, method="medianPolish", rescale=TRUE)
  design <- stats::model.matrix(out$config$frm, data=out$state$samples)
  cols_des <- colnames(design)
  
  frm0 <- stats::as.formula(paste("~0 + ", out$config$test_term))
  cols_pick <- colnames(stats::model.matrix(frm0, data=out$state$samples))
  cols_pick <- cols_pick[cols_pick %in% cols_des]
  idx <- which(cols_des %in% cols_pick)
  if(length(idx) != 1) f.err("test_deqms: length(idx) != 1", config=out$config)
  
  fit <- limma::lmFit(out$state$expression, design)
  fit <- limma::eBayes(fit, trend=trend)
  fit$count <- counts[rownames(fit$coefficients)]
  
  if(length(unique(fit$count)) >= 2) {
    fit <- DEqMS::spectraCounteBayes(fit, fit.method="loess")
    tbl <- DEqMS::outputResult(fit, coef_col=idx)
    return(tbl)
  } else {
    f.msg("WARNING: length(unique(fit$count)) < 2;", 
      "falling back to h0testr::test_trend(); unique(fit$count):", 
      unique(fit$count),
      config=config
    )      
    tbl <- test_trend(state, config)
  }
  
  return(tbl)
}

#' Hypothesis testing using the \code{msqrob2} package
#' @description
#'   Tests for differential expression using the \code{msqrob2::msqrob()} function.
#' @details
#'   Returns gene-level hypothesis testing results based on 
#'     peptide/precursor-level input.
#'   Flow is:
#'     \tabular{l}{
#'       1. Create a \code{QFeatures} object with \code{QFeatures::readQFeatures()}. \cr
#'       2. Convert \code{NA}s to zero using \code{QFeatures::zeroIsNA()}. \cr
#'       3. Make gene/protein-group expression with \code{QFeatures::aggregateFeatures()}. \cr
#'       4. Estimate model parameters using \code{msqrob2::msqrob()}. \cr
#'       5. Calculate test statistics with \code{msqrob2::hypothesisTest()}. \cr
#'     }
#'   See documentation for \code{h0testr::new_config()} 
#'     for more detailed description of configuration parameters. 
#' @param state List with elements like those returned by \code{read_data()}:
#'   \tabular{ll}{
#'     \code{expression} \cr \tab Numeric matrix with non-negative expression values. \cr
#'     \code{features}   \cr \tab A data.frame with feature meta-data for rows of expression. \cr
#'     \code{samples}    \cr \tab A data.frame with observation meta-data for columns of expression. \cr
#'   } 
#' @param config List with configuration values. Requires the following keys:
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
#' sim <- h0testr::sim2(
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
#' out <- h0testr::initialize(state, config, minimal=TRUE)
#' 
#' ## actual test:
#' tbl <- h0testr::test_msqrob(out$state, out$config)
#' head(tbl)

test_msqrob <- function(state, config, maxit=100) {

  check_config(config)
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
  if(length(idx) != 1) f.err("test_msqrob: length(idx) != 1", config=config)
  cols_pick <- cols_des[idx]
  
  con <- msqrob2::makeContrast(contrasts=paste0(cols_pick, "=0"), 
    parameterNames=c(cols_pick))
  
  obj <- msqrob2::hypothesisTest(object=obj, i="genes", contrast=con, 
    modelColumn="msqrobModels")
  
  dat <- SummarizedExperiment::rowData(obj[["genes"]])    
  tbl <- dat[[cols_pick]]
  dat <- dat[, c(config$gene_id_col, "nNonZero", ".n"), drop=F]
  dat <- as.data.frame(dat)
  tbl <- cbind(dat[rownames(tbl), , drop=F], tbl)
  tbl <- tbl[order(tbl$adjPval, -abs(tbl$logFC)), ]
  rownames(tbl) <- NULL
  
  return(tbl)
}

#' Hypothesis testing using the \code{proDA} package
#' @description
#'   Tests for differential expression using the \code{proDA::proDA()} function.
#' @details
#'   Uses the \code{proDA::proDA()} function. Returned results sorted by p-value.
#'   Returns peptide/precursor/gene-level hypothesis testing results based on 
#'     peptide/precursor/gene-level input. That is, testing is on features of 
#'     the input, so have to aggregate data to the desired level for 
#'     hypothesis testing first. Main feature of this method is its native 
#'     handling of missing values.
#'   See documentation for \code{h0testr::new_config()} 
#'     for more detailed description of configuration parameters. 
#' @param state List with elements like those returned by \code{read_data()}:
#'   \tabular{ll}{
#'     \code{expression} \cr \tab Numeric matrix with non-negative expression values. \cr
#'     \code{features}   \cr \tab A data.frame with feature meta-data for rows of expression. \cr
#'     \code{samples}    \cr \tab A data.frame with observation meta-data for columns of expression. \cr
#'   } 
#' @param config List with configuration values. Uses the following keys:
#'   \tabular{ll}{
#'     \code{gene_id_col}          \cr \tab Name of column in \code{state$features} with unique gene/protein group ids. \cr
#'     \code{feat_col}             \cr \tab Name of column in \code{state$features} corresponding to \code{rownames(state$expression)}. \cr
#'     \code{obs_col}              \cr \tab Name of column in \code{state$samples} corresponding to \code{colnames(state$expression)}. \cr
#'     \code{frm}                  \cr \tab Formula (formula) to be fit \cr
#'     \code{test_term}            \cr \tab Term (character) to be tested for non-zero coefficient. \cr
#'     \code{normalization_method} \cr \tab If present and \code{is_log_transformed} unset, used to infer it. \cr
#'   }
#' @param is_log_transformed Logical scalar indicating if \code{state$expression} has been log transformed.
#' @param prior_df Strictly positive count (\code{location_prior_df}) indicating number of dfs for prior.
#' @param maxit Strictly positive count indicating maximum number of iterations for \code{proDA::proDA()} algorithm.
#' @return
#'   A data.frame with test results. Columns include \code{config$gene_id_col} and: 
#'     \code{c("nNonZero .n", "logFC", "se", "df", "t", "pval", "adjPval")}.
#' @examples
#' ## lengthy setup of expression data:
#' set.seed(101)
#' nsamps <- 6
#' sim <- h0testr::sim2(
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
#' out <- h0testr::initialize(state, config, minimal=TRUE)
#' 
#' ## actual test:
#' tbl <- h0testr::test_proda(out$state, out$config, is_log_transformed=FALSE)
#' head(tbl)

test_proda <- function(state, config, is_log_transformed=NULL, prior_df=3, maxit=20) {
  
  if(is.null(is_log_transformed) || is_log_transformed %in% "") {
    if(is.null(config$normalization_method) || config$normalization_method %in% "") {
      f.err("test_proda: is_log_transformed and config$normalization_method both unset", 
        config=config)
    }
    if(config$normalization_method %in% c("none")) {
      is_log_transformed <- FALSE 
    } else {
      is_log_transformed <- TRUE
    }
  }
  
  fit <- proDA::proDA(state$expression, design=config$frm, col_data=state$samples, 
    data_is_log_transformed=is_log_transformed, location_prior_df=prior_df, max_iter=maxit)
  
  frm0 <- stats::as.formula(paste("~0 + ", config$test_term))
  cols <- colnames(stats::model.matrix(frm0, data=state$samples))
  i <- grepl(":", cols)
  if(any(i)) cols[i] <- paste0("`", cols[i], "`")
  
  cols2 <- proDA::result_names(fit)
  i <- cols %in% cols2
  if(sum(i) != 1) {
    f.err("test_proda: multi-column match; cols: ", cols, "; cols2: ", 
      cols2, config=config)
  } 
  col_pick <- cols[i]
  
  tbl <- proDA::test_diff(fit, contrast=col_pick, sort_by="pval")
  return(tbl)
}

#' Hypothesis testing using the \code{prolfq} package
#' @description
#'   Tests for differential expression using the \code{prolfq::build_model()} function.
#' @details
#'   Uses the \code{proDA::proDA()} function. Returned results sorted by p-value.
#'     Only works if all terms in \code{config$frm} are \code{factor} or 
#'     \code{character}. Error otherwise.
#'   Flow is:
#'     \tabular{l}{
#'       1. Reshape data into long format with intensities, feature meta, and sample meta. \cr
#'       2. Make and populate \code{prolfqua::AnalysisTableAnnotation} object. \cr
#'       3. Make \code{prolfqua::LFQData} object from data and 
#'            \code{prolfqua::AnalysisTableAnnotation} object. \cr
#'       4. Make \code{prolfqua::strategy_lm} object from \code{config$frm}. \cr
#'       5. Build \code{prolfqua} model from \code{prolfqua::LFQData} and 
#'            \code{prolfqua::strategy_lm} objects. \cr
#'       5. Return sub-table of ANOVA results corresponding to \code{config$test_term}. \cr
#'     }
#'   See documentation for \code{h0testr::new_config()} 
#'     for more detailed description of configuration parameters. 
#' @param state List with elements like those returned by \code{read_data()}:
#'   \tabular{ll}{
#'     \code{expression} \cr \tab Numeric matrix with non-negative expression values. \cr
#'     \code{features}   \cr \tab A data.frame with feature meta-data for rows of expression. \cr
#'     \code{samples}    \cr \tab A data.frame with observation meta-data for columns of expression. \cr
#'   }
#' @param config List with configuration values. Uses the following keys:
#'   \tabular{ll}{
#'     \code{gene_id_col}           \cr \tab Name of column in \code{state$features} with unique gene/protein group ids. \cr
#'     \code{feat_col}              \cr \tab Name of column in \code{state$features} corresponding to \code{rownames(state$expression)}. \cr
#'     \code{obs_col}               \cr \tab Name of column in \code{state$samples} corresponding to \code{colnames(state$expression)}. \cr
#'     \code{frm}                   \cr \tab Formula (formula) to be fit \cr
#'     \code{test_term}             \cr \tab Term (character) to be tested for non-zero coefficient. \cr
#'     \code{sample_factors}        \cr \tab List with one character vector per variable, with factor level ordering. \cr
#'     \code{normalization_method}  \cr \tab If present and \code{is_log_transformed} unset, used to infer it. \cr
#'   }
#' @param is_log_transformed Logical scalar indicating if \code{state$expression} has been log transformed.
#' @return
#'   A data.frame with test results. Columns include \code{config$gene_id_col} and: 
#'     \code{c("nNonZero .n", "logFC", "se", "df", "t", "pval", "adjPval")}.
#' @examples
#' ## lengthy setup of expression data:
#' set.seed(101)
#' nsamps <- 6
#' sim <- h0testr::sim2(
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
#'   sex=rep(c("M", "F"), nsamps)
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
#' out <- h0testr::initialize(state, config, minimal=TRUE)
#' 
#' ## actual test:
#' tbl <- h0testr::test_prolfqua(out$state, out$config, is_log_transformed=FALSE)
#' head(tbl)

test_prolfqua <- function(state, config, is_log_transformed=NULL) {
  
  if(is.null(is_log_transformed) || is_log_transformed %in% "") {
    if(is.null(config$normalization_method) || config$normalization_method %in% "") {
      f.err("test_prolfqua: is_log_transformed and config$normalization_method both unset", 
        config=config)
    }
    if(config$normalization_method %in% c("none")) {
      is_log_transformed <- FALSE 
    } else {
      is_log_transformed <- TRUE
    }
  }
  
  dat <- data.frame(
    state$features[, c(config$gene_id_col, config$feat_id_col)], 
    state$expression
  )
  
  dat <- stats::reshape(dat, direction="long", varying=colnames(state$expression), 
    v.names="intensity", idvar=c(config$gene_id_col, config$feat_id_col), 
    timevar="sample", times=colnames(state$expression))
  rownames(dat) <- NULL
  
  samps <- state$samples
  if(any(duplicated(samps[[config$obs_col]]))) {
    f.err("test_prolfqua: duplicated state$samples[, config$obs_col]; obs_col: ", 
      config$obs_col, config=config)
  }
  rownames(samps) <- samps[[config$obs_col]]
  if(!all(dat$sample %in% rownames(samps))) {
    i <- !(dat$sample %in% rownames(samps))
    f.err("test_prolfqua: observation mismatch: ", 
      paste(dat$sample[i], collapse=", "), config=config)
  }
  
  frm <- paste(as.character(config$frm), collapse="")
  frm <- gsub("[[:space:]]", "", frm)
  frm <- sub("~", "", frm)
  trms <- sort(unique(unlist(strsplit(frm, "[+*:]"))))
  
  meta <- prolfqua::AnalysisTableAnnotation$new()
  meta$workIntensity <- "intensity"
  meta$is_response_transformed <- is_log_transformed
  meta$hierarchy[[config$gene_id_col]] <- config$gene_id_col
  meta$hierarchy[[config$feat_id_col]] <- config$feat_id_col
  
  for(trm in trms) {
    if(!(trm %in% names(config$sample_factors))) {
      f.err("test_prolfqua: !(trm %in% config$sample_factors); trm: ", 
        trm, config=config)
    }
    meta$factors[[trm]] <- trm
    dat[, trm] <- as.character(samps[dat$sample, trm, drop=T])
  }
  obj <- prolfqua::LFQData$new(data=dat, config=meta)
  
  frm <- paste(c("intensity", as.character(config$frm)), collapse=" ")
  strgy <- prolfqua::strategy_lm(frm)
  
  model <- prolfqua::build_model(data=obj$data, model_strategy=strgy, 
    subject_Id=obj$config$hierarchy_keys())
  
  tbl <- as.data.frame(model$get_anova())
  
  return(tbl[tbl$factor %in% config$test_term, , drop=F])
}

#' Hypothesis testing using \code{limma::voom}
#' @description
#'   Tests for differential expression using the \code{limma::voom()} function.
#' @details
#'   The \code{limma::voom()} model is fit to \code{config$frm} and an F-test 
#'     is performed for whether the effect of \code{config$test_term} on 
#'     \code{state$expression} is zero. 
#'   See documentation for \code{h0testr::new_config()} 
#'     for more detailed description of configuration parameters. 
#' @param state List with elements like those returned by \code{read_data()}:
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
#' exprs <- h0testr::sim2(n_samps1=6, n_samps2=6, n_genes=25, 
#'   n_genes_signif=5, fold_change=2, mnar_c0=-Inf, mnar_c1=0, mcar_p=0)$mat
#' exprs <- log2(exprs + 1)
#' feats <- data.frame(feature_id=rownames(exprs))
#' samps <- data.frame(observation_id=colnames(exprs), 
#'   condition=c(rep("placebo", 6), rep("drug", 6)))
#' state <- list(expression=exprs, features=feats, samples=samps)
#' 
#' config <- h0testr::new_config()    ## defaults
#' config$save_state <- FALSE           ## default is TRUE
#' config$feat_col <- config$feat_id_col <- config$gene_id_col <- "feature_id"
#' config$obs_col <- config$obs_id_col <- config$sample_id_col <- "observation_id"
#' config$frm <- ~condition
#' config$test_term <- "condition"
#' config$sample_factors <- list(condition=c("placebo", "drug"))
#'
#' ## set up and check configuration, including covariates:
#' out <- h0testr::initialize(state, config, minimal=TRUE)
#' 
#' tbl <- h0testr::test_voom(out$state, out$config)
#' print(tbl)

test_voom <- function(state, config, normalize.method="none") {

  check_config(config)
  f.check_state(state, config)
  
  exprs <- state$expression
  i <- apply(exprs, 1, function(v) any(is.na(v)))
  exprs <- exprs[!i, , drop=F]
  
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

#' Hypothesis testing using \code{limma} trend
#' @description
#'   Test for differential expression using \code{limma::eBayes(trend=TRUE)}.
#' @details
#'   Wrapper for the \code{limma} pakcage.
#'   Flow is:
#'     \tabular{l}{
#'       1. Fit linear model with \code{lmFit()}. \cr
#'       2. Compute moderated statistics with \code{eBayes(trend=TRUE)}. \cr
#'       3. Generate a \code{data.frame} with results using \code{topTable()}. \cr
#'     }
#'  The \code{lmFit()} model is fit to \code{config$frm} and an F-test is 
#'    performed on each \code{config$feat_col} for whether the effect of 
#'    \code{config$test_term} on \code{state$expression} is zero. 
#'   See documentation for \code{h0testr::new_config()} 
#'     for more detailed description of configuration parameters. 
#' @param state List with elements formatted like the list returned by \code{read_data()}:
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
#' exprs <- h0testr::sim2(n_samps1=6, n_samps2=6, n_genes=25, 
#'   n_genes_signif=5, fold_change=2, mnar_c0=-Inf, mnar_c1=0, mcar_p=0)$mat
#' exprs <- log2(exprs + 1)
#' feats <- data.frame(feature_id=rownames(exprs))
#' samps <- data.frame(observation_id=colnames(exprs), 
#'   condition=c(rep("placebo", 6), rep("drug", 6)))
#' state <- list(expression=exprs, features=feats, samples=samps)
#' 
#' config <- h0testr::new_config()    ## defaults
#' config$save_state <- FALSE           ## default is TRUE
#' config$feat_col <- config$feat_id_col <- config$gene_id_col <- "feature_id"
#' config$obs_col <- config$obs_id_col <- config$sample_id_col <- "observation_id"
#' config$frm <- ~condition
#' config$test_term <- "condition"
#' config$sample_factors <- list(condition=c("placebo", "drug"))
#'
#' ## set up and check covariates and parameters:
#' out <- h0testr::initialize(state, config, minimal=TRUE)
#' 
#' tbl <- h0testr::test_trend(state, config)
#' print(tbl)

test_trend <- function(state, config) {

  check_config(config)
  f.check_state(state, config)

  if(!is.matrix(state$expression)) {
    f.err("test_trend: !is.matrix(state$expression)", config=config)
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

## helper for test():

f.format_msqrob <- function(tbl, config) {
  
  if(!is.data.frame(tbl)) {
    f.err("f.format_msqrob: !is.data.frame(tbl); class(tbl): ", 
      class(tbl), config=config)
  }
  
  if(!(config$gene_id_col %in% names(tbl))) {
    f.err("f.format_msqrob: expected names not %in% names(tbl); names(tbl):", 
      names(tbl), config=config)
  }
  
  nom <- c("logFC", "t", "pval", "adjPval")
  if(!all(nom %in% names(tbl))) {
    f.err("f.format_msqrob: expected names not %in% names(tbl); names(tbl):", 
      names(tbl), config=config)
  }
  
  tbl <- data.frame(feature=tbl[[config$gene_id_col]], expr=as.numeric(NA), 
    logfc=tbl$logFC, stat=tbl$t, lod=as.numeric(NA), 
    pval=tbl$pval, adj_pval=tbl$adjPval)
  
  tbl <- tbl[order(tbl$pval, decreasing=F), , drop=F]
  rownames(tbl) <- NULL
  
  return(tbl)
}

## helper for test():

f.format_proda <- function(tbl, config) {
  
  if(!is.data.frame(tbl)) {
    f.err("f.format_proda: !is.data.frame(tbl); class(tbl): ", 
      class(tbl), config=config)
  }

  nom <- c("name", "avg_abundance", "diff", "t_statistic", "pval", "adj_pval")
  if(!all(nom %in% names(tbl))) {
    f.err("f.format_proda: expected names not %in% names(tbl); names(tbl):", 
      names(tbl), config=config)
  }
  
  tbl <- data.frame(feature=tbl$name, expr=tbl$avg_abundance, 
    logfc=tbl$diff, stat=tbl$t_statistic, lod=as.numeric(NA), 
    pval=tbl$pval, adj_pval=tbl$adj_pval)
  
  tbl <- tbl[order(tbl$pval, decreasing=F), , drop=F]
  rownames(tbl) <- NULL
  
  return(tbl)
}

## helper for test():

f.format_prolfqua <- function(tbl, config) {
  
  if(!is.data.frame(tbl)) {
    f.err("f.format_prolfqua: !is.data.frame(tbl); class(tbl): ", 
      class(tbl), config=config)
  }
  
  if(!(config$feat_col %in% names(tbl))) {
    f.err("f.format_prolfqua: config$feat_col:", config$feat_col, 
      "not %in% names(tbl); names(tbl):", names(tbl), config=config)
  }
  
  nom <- c("F.value", "p.value", "FDR")
  if(!all(nom %in% names(tbl))) {
    f.err("f.format_prolfqua: expected names not %in% names(tbl); names(tbl):", 
      names(tbl), config=config)
  }
  
  tbl <- data.frame(feature=tbl[[config$feat_col]], expr=as.numeric(NA), 
    logfc=as.numeric(NA), stat=tbl$F.value, lod=as.numeric(NA), 
    pval=tbl$p.value, adj_pval=tbl$FDR)
  
  tbl <- tbl[order(tbl$pval, decreasing=F), , drop=F]
  rownames(tbl) <- NULL
  
  return(tbl)
}

## helper for test():

f.format_limma <- function(tbl, config) {
  
  if(!is.data.frame(tbl)) {
    f.err("f.format_limma: !is.data.frame(tbl); class(tbl): ", 
      class(tbl), config=config)
  }
    
  nom <- c("AveExpr", "logFC", "t", "B", "P.Value", "adj.P.Val")
  if(!all(nom %in% names(tbl))) {
    f.err("f.format_limma: expected names not %in% names(tbl); names(tbl):", 
      names(tbl), config=config)
  }
  
  tbl <- data.frame(feature=rownames(tbl), expr=tbl$AveExpr, 
    logfc=tbl$logFC, stat=tbl$t, lod=tbl$B, 
    pval=tbl$P.Value, adj_pval=tbl$adj.P.Val)
  
  tbl <- tbl[order(tbl$pval, decreasing=F), , drop=F]
  rownames(tbl) <- NULL
  
  return(tbl)
}

#' Get vector of \code{test(method=)} method options
#' @description
#'   Get a vector with acceptable values of \code{method} parameter for \code{h0testr::test()}.
#' @return
#'   Character vector with names of acceptable values for \code{h0testr::normalize(method=)}.
#' @examples
#' test_methods <- h0testr::test_methods()
#' cat("Available test methods:\n")
#' for(method in test_methods) {
#'   cat("method:", method, "\n")
#' }

test_methods <- function() {
  return(
    c("trend", "deqms", "msqrob", "proda", "prolfqua", "voom")
  )
}

#' Hypothesis testing
#' @description
#'   Wrapper for various hypothesis testing methods.
#' @details
#'   Tests for differential expression using method specified in config. 
#'   See invididual \code{test_*} methods for more details. 
#'   The \code{method} setting meanings are: 
#'   \tabular{ll}{
#'     \code{trend}  \cr \tab Use \code{limma::eBayes(trend=TRUE)}. \cr
#'     \code{deqms}  \cr \tab Use \code{DEqMS::spectraCounteBayes()}. \cr
#'     \code{msqrob} \cr \tab Use \code{msqrob2::msqrob()}. \cr
#'     \code{proda}  \cr \tab Use \code{proDA::proDA()}. \cr
#'     \code{voom}   \cr \tab Use \code{limma::voom()}. \cr
#'   } 
#'   See documentation for \code{h0testr::new_config()} 
#'     for more detailed description of configuration parameters. 
#' @param state List with elements formatted like the list returned by \code{read_data()}:
#'   \tabular{ll}{
#'     \code{expression} \cr \tab Numeric matrix with non-negative expression values. \cr
#'     \code{features}   \cr \tab A data.frame with feature meta-data for rows of expression. \cr
#'     \code{samples}    \cr \tab A data.frame with observation meta-data for columns of expression. \cr
#'   } 
#' @param config List with configuration values. Uses the following keys:
#'   \tabular{ll}{
#'     \code{feat_col}       \cr \tab Name of column in \code{state$fetaures} matching \code{rownames(state$expression)}. \cr
#'     \code{obs_col}        \cr \tab Name of column in \code{state$samples} matching \code{colnames(state$expression)}. \cr
#'     \code{gene_id_col}    \cr \tab Name of column in \code{state$fetaures} with gene/protein-group ids. \cr
#'     \code{frm}            \cr \tab Formula (formula) to be fit \cr
#'     \code{test_term}      \cr \tab Term (character scalar) to be tested for non-zero coefficient. \cr
#'     \code{sample_factors} \cr \tab List specifying levels of factor variables in \code{config$frm} (see examples). \cr
#'     \code{test_method}    \cr \tab Character scalar in \code{c("voom", "trend", "deqms", "msqrob", "proda"}. \cr
#'   }
#' @param method Name of test method where 
#'   \code{method \%in\% h0testr::test_methods()}.
#' @param is_log_transformed Logical scalar: if \code{state$expression} has 
#'   been log transformed. Required if 
#'   \code{method \%in\% c("proda", "prolfqua")}.
#' @param prior_df Prior degrees of freedom for method \code{proda}; 
#'   where \code{2 <= prior_df <= n_features}.
#' @return A list with the following elements: \cr
#'   \tabular{ll}{
#'     \code{original} \cr \tab A \code{data.frame} with results in native format returned by test. \cr
#'     \code{standard} \cr \tab A \code{data.frame} with results in a standardized format \cr
#'   } \cr
#'   The \code{standard} \code{data.frame} has the following fields: \cr
#'   \tabular{ll}{
#'     \code{feature}   \cr \tab Name of feature tested. \cr
#'     \code{expr}      \cr \tab Average feature expression. \cr
#'     \code{logfc}     \cr \tab Estimated log fold-change between conditions. \cr
#'     \code{stat}      \cr \tab Value of test statistic. \cr
#'     \code{lod}       \cr \tab Log-odds of differential expression. \cr
#'     \code{pval}      \cr \tab Raw p-value resulting from test. \cr
#'     \code{adj_pval}  \cr \tab Adjusted (for multiple testing) p-value. \cr
#'   }
#' @examples
#' set.seed(101)
#' ## no missing values: mnar_c0=-Inf, mnar_c1=0, mcar_p=0
#' exprs <- h0testr::sim2(n_samps1=6, n_samps2=6, n_genes=25, 
#'   n_genes_signif=5, fold_change=2, mnar_c0=-Inf, mnar_c1=0, mcar_p=0)$mat
#' exprs <- log2(exprs + 1)
#' feats <- data.frame(feature_id=rownames(exprs))
#' samps <- data.frame(observation_id=colnames(exprs), 
#'   condition=c(rep("placebo", 6), rep("drug", 6)))
#' state <- list(expression=exprs, features=feats, samples=samps)
#' 
#' config <- list(feat_id_col="feature_id", gene_id_col="feature_id", 
#'   obs_id_col="observation_id", sample_id_col="observation_id", 
#'   frm=~condition, test_term="condition",
#'   sample_factors=list(condition=c("placebo", "drug"))
#' )
#' 
#' ## set up and check covariates and parameters:
#' out <- h0testr::initialize(state, config, minimal=TRUE)
#' state <- out$state
#' config <- out$config
#' 
#' out <- h0testr::test(state, config, method="trend")
#' head(out$original)
#' head(out$standard)

test <- function(state, config, method=NULL, 
    is_log_transformed=NULL, prior_df=NULL) {
  
  if(is.null(method) || method %in% "") method <- config$test_method
  if(is.null(method) || method %in% "") {
    f.err("test: method and config$test_method both unset", config=config)
  }
  if(method %in% c("proda", "prolfqua") && !is.logical(is_log_transformed)) {
    f.err("test: method %in% c('proda', 'prolfqua') && !is.logical(is_log_transformed)", 
      config=config)
  }
  if(is.null(prior_df)) prior_df <- config$test_prior_df
  
  if(method %in% "trend") {
    tbl <- test_trend(state, config)
    tbl2 <- f.format_limma(tbl, config)
  } else if(method %in% "deqms") {
    tbl <- test_deqms(state, config)
    tbl2 <- f.format_limma(tbl, config)
  } else if(method %in% "msqrob") {
    tbl <- test_msqrob(state, config)
    tbl2 <- f.format_msqrob(tbl, config)
  } else if(method %in% "proda") {
    tbl <- test_proda(state, config, 
      is_log_transformed=is_log_transformed, prior_df=prior_df)
    tbl2 <- f.format_proda(tbl, config)
  } else if(method %in% "prolfqua") {
    tbl <- test_prolfqua(state, config, 
      is_log_transformed=is_log_transformed)
    tbl2 <- f.format_prolfqua(tbl, config)
  } else if(method %in% "voom") {
    tbl <- test_voom(state, config)
    tbl2 <- f.format_limma(tbl, config)
  } else if(method %in% "none") {
    f.msg("skipping testing: method %in% 'none'", config=config)
    return(NULL)
  } else f.err("test: unexpected method:", config$test_method, config=config)
  
  feats <- state$features
  if(method %in% c("deqms", "msqrob")) {
    test_col <- config$gene_id_col
    feats <- feats[!duplicated(feats[[test_col]]), , drop=F]
    if(!(config$feat_id_col %in% config$gene_id_col)) {
      feats[[config$feat_id_col]] <- NULL
    }
  } else {
    test_col <- config$feat_col
  }
  rownames(feats) <- feats[[test_col]]
  if(!all(tbl2$feature %in% rownames(feats))) {
    f.err("test: !all(tbl2$feature %in% rownames(feats))", config=config)
  }
  
  tbl0 <- feats[rownames(tbl), , drop=F]
  tbl <- cbind(tbl0, tbl)
  rownames(tbl) <- NULL
  
  if((!is.null(config$save_state)) && config$save_state) {
    file_out <- paste0(config$dir_out, "/", length(config$run_order) + 3, 
      config$result_mid_out, config$suffix_out)
    f.log("writing results to", file_out, config=config)
    f.save_tsv(tbl, file_out)
  }
  
  return(list(original=tbl, standard=tbl2))
}
