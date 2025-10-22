#' Combine replicate observations
#' @description
#'   Combine expression signals from replicate observations, such as technical replicates.
#' @details Combines signals for each feature across technical replicates by taking median.
#'   Deletes \code{state$samples[, config$obs_id_col} if not same as 
#'     \code{config$sample_id_col}.
#'   Sets \code{config$obs_col} and \code{config$obs_id_col} to 
#'     \code{config$sample_id_col}.
#'   See documentation for \code{h0testr::new_config()} 
#'     for more detailed description of configuration parameters. 
#' @param state List with elements formatted like the list returned by 
#'   \code{read_data()}:
#'   \tabular{ll}{
#'     \code{expression} \cr \tab Numeric matrix with non-negative expression values. \cr
#'     \code{features}   \cr \tab A data.frame with feature meta-data for rows of expression. \cr
#'     \code{samples}    \cr \tab A data.frame with observation meta-data for columns of expression. \cr
#'   } 
#' @param config List with configuration values. Uses the following keys:
#'   \tabular{ll}{
#'     \code{feat_col}        \cr \tab Name of column in \code{state$features} matching \code{rownames(state$expression)}. \cr
#'     \code{obs_id_col}      \cr \tab Name of column in \code{state$samples} matching \code{colnames(state$expression)}. \cr
#'     \code{sample_id_col}   \cr \tab Name of column in \code{state$samples} with unique sample (case) ids. \cr
#'     \code{save_state}      \cr \tab Logical scalar indicating whether to save new state to disk. \cr
#'     \code{dir_out}         \cr \tab Output directory path (scalar character); only needed if \code{save_state == TRUE}. \cr
#'     \code{data_mid_out}    \cr \tab Midfix of expression matrix filename; only needed if \code{save_state == TRUE}.\cr
#'     \code{feature_mid_out} \cr \tab Midfix of feature metadata filename; only needed if \code{save_state == TRUE}.\cr
#'     \code{sample_mid_out}  \cr \tab Midfix of sample metadata filename; only needed if \code{save_state == TRUE}.\cr
#'     \code{suffix_out}      \cr \tab Suffix of output files; only needed if \code{save_state == TRUE}.\cr
#'   }
#' @param fn Function for combining replicated measurements.
#' @return A list (the processed state) with the following elements:
#'   \tabular{ll}{
#'     \code{expression} \cr \tab Numeric matrix with non-negative expression values. \cr
#'     \code{features}   \cr \tab A data.frame with feature meta-data for rows of expression. \cr
#'     \code{samples}    \cr \tab A data.frame with observation meta-data for columns of expression. \cr
#'   } 
#' @examples
#' set.seed(101)
#' exprs <- h0testr::sim1(n_obs=6, n_feats=8)$mat
#' feats <- data.frame(feature_id=rownames(exprs))
#' samps <- data.frame(observation_id=colnames(exprs))
#' samps$sample_id=c("samp1", "samp1", "samp2", "samp2", "samp3", "samp3")
#' state <- list(expression=exprs, features=feats, samples=samps)
#' config <- list(feat_col="feature_id", obs_id_col="observation_id", sample_id_col="sample_id", save_state=FALSE)
#' out <- h0testr::combine_replicates(state, config)
#' print(state)
#' print(out$state)
#' str(config)
#' str(out$config)

combine_replicates <- function(state, config, fn=stats::median) {
  
  check_config(config)
  
  if(is.null(config$sample_id_col) || !(config$sample_id_col %in% names(state$samples))) {
    f.err("combine_replicates: !(config$sample_id_col %in% names(state$samples));",
      "config$sample_id_col:", config$sample_id_col, 
      "; names(state$samples):", names(state$samples), config=config)
  }
  
  if(config$sample_id_col %in% config$obs_id_col) {
    f.msg("combine_replicates: config$sample_id_col %in% config$obs_id_col; ", 
      "returning unchanged state and updated config.", config=config)
    config$obs_col <- config$obs_id_col <- config$sample_id_col
  } else {
  
    f <- function(v, s) tapply(v, s, fn, na.rm=T)
    sample_ids <- state$samples[[config$sample_id_col]]
    state$expression <- t(apply(state$expression, 1, f, sample_ids))
    state$samples <- state$samples[!duplicated(sample_ids), , drop=F]
    
    if(!is.null(config$obs_id_col)) {
      if(config$obs_id_col != config$sample_id_col) {
        state$samples[[config$obs_id_col]] <- NULL
      }
    }
    config$obs_col <- config$obs_id_col <- config$sample_id_col
    
    sample_ids <- as.character(state$samples[[config$sample_id_col]])
    if(!all(sample_ids %in% colnames(state$expression))) {
      f.err("combine_replicates: !all(samples[[config$sample_id_col]] %in% colnames(expression))", 
        config=config)
    }  
    state$expression <- state$expression[, as.character(sample_ids), drop=F]
  }
  
  f.check_state(state, config)
  f.report_state(state, config)
  
  if(!is.null(config$run_order)) {
    i <- config$run_order %in% "combine_replicates"
    if(any(i)) {
      prfx <- paste0(which(i)[1] + 2, ".combined_replicates")
    } else {
      prfx <- "combined_replicates"
    }
  }
  f.save_state(state, config, prefix=prfx)
  
  return(list(state=state, config=config))
}

## helper for combine_features(); uses config$gene_id_col; 
##   sets config$feat_col and config$feat_id_col to config$gene_id_col:

f.combine_features_median_polish <- function(state, config, maxit=30) {
  
  feats <- state$features
  genes <- feats[[config$gene_id_col]]
  i <- is.na(genes) | genes %in% ""
  feats[[config$gene_id_col]][i] <- paste0("unknown_", feats[[config$feat_id_col]][i])
  genes <- feats[[config$gene_id_col]]
  feats <- feats[!duplicated(genes), , drop=F]
  
  if(!is.null(config$feat_id_col)) {
    if(config$feat_id_col != config$gene_id_col) {
      feats[[config$feat_id_col]] <- NULL
    }
  }
  config$feat_col <- config$feat_id_col <- config$gene_id_col
  
  f <- function(idxs) {
    v <- MsCoreUtils::medianPolish(state$expression[idxs, , drop=F], 
      na.rm=T, maxiter=maxit)
    if(any(v <= 0, na.rm=T)) v <- v + 2 * abs(min(v, na.rm=T)) + 1
    return(v)
  }
  exprs <- tapply(1:nrow(state$expression), genes, f)
  exprs <- do.call(rbind, exprs)
  
  genes <- as.character(feats[[config$gene_id_col]])
  exprs <- exprs[as.character(genes), , drop=F]
  state <- list(expression=exprs, features=feats, samples=state$samples)
  
  return(list(state=state, config=config))
}

## helper for combine_features(); uses config$gene_id_col;
##   sets config$feat_col and config$feat_id_col to config$gene_id_col:

f.combine_features_robust_summary <- function(state, config) {
  
  feats <- state$features
  genes <- feats[[config$gene_id_col]]
  i <- is.na(genes) | genes %in% ""
  feats[[config$gene_id_col]][i] <- paste0("unknown_", feats[[config$feat_id_col]][i])
  genes <- feats[[config$gene_id_col]]
  feats <- feats[!duplicated(genes), , drop=F]
  
  if(!is.null(config$feat_id_col)) {
    if(config$feat_id_col != config$gene_id_col) {
      feats[[config$feat_id_col]] <- NULL
    }
  }
  config$feat_col <- config$feat_id_col <- config$gene_id_col
  
  f <- function(idxs) {
    x <- state$expression[idxs, , drop=F]
    i <- apply(x, 1, function(v) !all(is.na(v) | v %in% 0))
    v <- MsCoreUtils::robustSummary(x[i, , drop=F], na.rm=T)
    if(any(v <= 0, na.rm=T)) v <- v + 2 * abs(min(v, na.rm=T)) + 1
    return(v)
  }
  exprs <- tapply(1:nrow(state$expression), genes, f)
  nom <- names(exprs)
  exprs <- do.call(rbind, exprs)
  rownames(exprs) <- nom
  
  genes <- as.character(feats[[config$gene_id_col]])
  exprs <- exprs[as.character(genes), , drop=F]
  
  state <- list(expression=exprs, features=feats, samples=state$samples)
  
  return(list(state=state, config=config))
}

#' Combine low-level features into high-level features
#' @description
#'   Combine e.g. peptide signals into gene group signals
#' @details
#'   Combines signals from lower-level features, like peptides, 
#'     into higher level features, like gene groups or protein groups.
#'   If \code{config$gene_id_col == config$feat_col}, then no changes
#'     made to \code{state} or \code{config}.
#'   Sets \code{config$feat_col} and \code{config$feat_id_col} to 
#'     \code{config$gene_id_col}.
#'   If \code{rescale \%in\% TRUE}, and \code{!(method \%in\% "none")}, 
#'     peptides rescaled by dividing by peptide mean (after removing 
#'     \code{NA}s).
#'   For \code{method="medianPolish"}, for each unique gene in 
#'     \code{state$features[, config$gene_id_col]}, the submatrix of 
#'       corresponding peptide signals across all samples is decomposed into: 
#'       \code{pep_exprs == median_column_effect + median_row_effect + overall_median},
#'     Then \code{median_column_effect} is returned, after shifted to 
#'       ensure all values strictly positive: \code{all(expression > 0, na.rm=TRUE)}.
#'   Uses \code{MsCoreUtils::medianPolish} and \code{MsCoreUtils::robustSummary}.
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
#'     \code{obs_col}         \cr \tab Name of column in \code{state$samples} matching \code{colnames(state$expression)}. \cr
#'     \code{feat_col}        \cr \tab Name of column in \code{state$features} matching \code{rownames(state$expression)}. \cr
#'     \code{gene_id_col}     \cr \tab Name of column in \code{state$features} specifying high-level feature for aggregation. \cr
#'     \code{save_state}      \cr \tab Logical scalar indicating whether to save new state to disk. \cr
#'     \code{dir_out}         \cr \tab Output directory path (scalar character); only needed if \code{save_state == TRUE}. \cr
#'     \code{data_mid_out}    \cr \tab Midfix of expression matrix filename; only needed if \code{save_state == TRUE}.\cr
#'     \code{feature_mid_out} \cr \tab Midfix of feature metadata filename; only needed if \code{save_state == TRUE}.\cr
#'     \code{sample_mid_out}  \cr \tab Midfix of sample metadata filename; only needed if \code{save_state == TRUE}.\cr
#'     \code{suffix_out}      \cr \tab Suffix of output files; only needed if \code{save_state == TRUE}.\cr
#'   }
#' @param method Name (character scalar) of method to use for combining, 
#'   where \code{method \%in\% c("medianPolish", "robustSummary", "none")}. Default: \code{"medianPolish"}.
#' @param rescale Logical scalar indicating whether to rescale peptides prior to aggregation. Default: \code{FALSE}.
#' @return A list (the processed state) with the following elements:
#'   \tabular{ll}{
#'     \code{expression} \cr \tab Numeric matrix with non-negative expression values. \cr
#'     \code{features}   \cr \tab A data.frame with feature meta-data for rows of expression. \cr
#'     \code{samples}    \cr \tab A data.frame with observation meta-data for columns of expression. \cr
#'   } 
#' @examples
#' ## set up data for examples:
#' set.seed(101)
#' exprs <- h0testr::sim2(n_samps1=3, n_samps2=3, n_genes=5, n_genes_signif=1, peps_per_gene=3, reps_per_sample=1)$mat
#' tmp <- strsplit(rownames(exprs), "_")
#' feats <- data.frame(pep=rownames(exprs), gene=sapply(tmp, "[", 1))
#' tmp <- strsplit(colnames(exprs), "_")
#' samps <- data.frame(obs=colnames(exprs), grp=sapply(tmp, "[", 1))
#' state <- list(expression=exprs, features=feats, samples=samps)
#' config <- list(feat_id_col="pep", gene_id_col="gene", obs_col="obs", feat_col="pep",
#'   feature_aggregation="medianPolish", save_state=FALSE)
#' rm(tmp, exprs, feats, samps)
#' print(state)
#' str(config)
#'
#' ## combine peps using config$feature_aggregation ("medianPolish"):
#' out <- h0testr::combine_features(state, config)
#' print(out$state)
#' str(out$config)
#' 
#' ## combine peps, overriding config$feature_aggregation:
#' out <- h0testr::combine_features(state, config, method="robustSummary")
#' print(out$state)
#' str(out$config)
#'
#' ## combine peps with rescaling:
#' out <- h0testr::combine_features(state, config, rescale=TRUE)
#' print(out$state)
#' str(out$config)
#'
#' ## with method="none":
#' out <- h0testr::combine_features(state, config, method="none")
#' print(out$state)
#' str(out$config)

combine_features <- function(state, config, method=NULL, rescale=FALSE) {
  
  check_config(config)
  
  if(is.null(config$gene_id_col) || !(config$gene_id_col %in% names(state$features))) {
    f.err("combine_features: config$gene_id_col %in% names(state$features);",
      "config$gene_id_col:", config$gene_id_col, 
      "; names(state$features):", names(state$features), config=config)
  }
  
  if(is.null(method) || method %in% "") method <- config$feature_aggregation
  if(is.null(method) || method %in% "") method <- "medianPolish"
  if(is.null(rescale)) rescale <- config$feature_aggregation_scaled
  if(is.null(rescale)) rescale <- FALSE
  
  if(config$gene_id_col %in% config$feat_col) {
    f.msg("combine_features: config$gene_id_col %in% config$feat_col; ", 
      "returning unchanged state and updated config.", config=config)
    config$feat_col <- config$gene_id_col <- config$feat_id_col
    method <- "none"
  }     
  
  if(rescale && !(method %in% "none")) {
    f.msg("combine_features: rescaling state$expression", config=config)
    f <- function(v) v / mean(v, na.rm=T)
    state$expression <- t(apply(state$expression, 1, f))
  }
  
  f.msg("combine_features: method:", method, "; rescale:", rescale, config=config)
  
  if(method %in% "medianPolish") {
    out <- f.combine_features_median_polish(state, config)
  } else if(method %in% "robustSummary") {
    out <- f.combine_features_robust_summary(state, config)
  } else if(method %in% "none") {
    config$feat_col <- config$gene_id_col <- config$feat_id_col
    out <- list(state=state, config=config)
  } else {
    f.err("combine_features: unexpected method: ", method, 
      '; should be one of: c("medianPolish", "robustSummary", "none")', 
      config=config
    )
  }
  state <- out$state
  config <- out$config
  
  f.check_state(state, config)
  f.report_state(state, config)
  
  if(!is.null(config$run_order)) {
    i <- config$run_order %in% "combine_features"
    if(any(i)) {
      prfx <- paste0(which(i)[1] + 2, ".combined_features")
    } else {
      prfx <- "combined_features"
    }
  }
  f.save_state(state, config, prefix=prfx)
  
  return(list(state=state, config=config))
}
