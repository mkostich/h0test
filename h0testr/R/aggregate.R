#' Combine replicate observations
#' @description
#'   Combine expression signals from replicate observations, such as technical replicates.
#' @details Combines signals for each feature across technical replicates by taking median.
#'   \code{fn} is applied with \code{na.rm=TRUE}, so a replicate group in which a
#'     feature was never measured is combined from no values at all, and what
#'     \code{fn} returns for that is up to \code{fn}: \code{stats::median} gives
#'     \code{NA}, but \code{sum} gives \code{0} and \code{mean} gives \code{NaN}.
#'     Such groups are therefore identified before \code{fn} runs and set to
#'     \code{NA} afterwards, whatever \code{fn} returned, so that a missing value
#'     is never recorded as a measurement of zero. \code{NA} is the only indicator
#'     of a missing value; see \code{h0testr::initialize()}.
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
  
    ## fn is applied with na.rm=T, so what it returns for a replicate group whose
    ##   every value is missing is up to fn, and only stats::median gets it right:
    ##   sum() returns 0 and mean() returns NaN. A 0 there is a missing value
    ##   written as a real one, and the most extreme value in the matrix on a log
    ##   scale, which every downstream step then reads as a measurement. NA is the
    ##   only indicator of a missing value; see f.zeros_to_na(). So the all-missing
    ##   groups are found before fn runs and restored after, which does not depend
    ##   on fn returning anything in particular:

    f <- function(v, s) {
      out <- tapply(v, s, fn, na.rm=T)
      n_obs <- tapply(v, s, function(x) sum(!is.na(x)))
      out[n_obs %in% 0] <- NA
      return(out)
    }
    sample_ids <- state$samples[[config$sample_id_col]]

    ## metadata are collapsed by keeping the first observation of each sample
    ##   below, so a covariate that varies between the replicates of one sample
    ##   would be silently taken from that first observation, quietly changing
    ##   the design; that is an error instead:

    if(!is.null(config$frm)) {
      vars <- sort(unique(f.parse_frm(config$frm, config)$vars))
      vars <- vars[vars %in% names(state$samples)]
      for(nom in vars) {
        v <- as.character(state$samples[[nom]])
        n_vals <- tapply(v, sample_ids, function(x) length(unique(x)))
        i <- !is.na(n_vals) & n_vals > 1
        if(any(i)) {
          sid <- names(n_vals)[i][1]
          j <- sample_ids %in% sid
          f.err("combine_replicates: covariate", nom, "varies between the",
            "replicates of a sample, so cannot be combined;", "\n",
            "config$sample_id_col:", config$sample_id_col, "; n samples",
            "affected:", sum(i), "\n", "first affected sample:", sid,
            "; its", nom, "values:", v[j], config=config)
        }
      }
    }

    state$expression <- t(apply(state$expression, 1, f, sample_ids))
    state$samples <- state$samples[!duplicated(sample_ids), , drop=F]
    
    if(!is.null(config$obs_id_col)) {
      if(config$obs_id_col != config$sample_id_col) {
        state$samples[[config$obs_id_col]] <- NULL
      }
    }
    config$obs_col <- config$obs_id_col <- config$sample_id_col
    
    sample_ids <- state$samples[[config$sample_id_col]]
    if(!all(sample_ids %in% colnames(state$expression))) {
      f.err("combine_replicates: !all(samples[[config$sample_id_col]] %in% colnames(expression))", 
        config=config)
    }  
    state$expression <- state$expression[, sample_ids, drop=F]
  }
  
  f.check_state(state, config)
  f.report_state(state, config)
  
  prfx <- "combined_replicates"
  if(!is.null(config$run_order)) {
    i <- config$run_order %in% "combine_replicates"
    if(any(i)) {
      prfx <- paste0(which(i)[1] + 2, ".combined_replicates")
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
  
  ## medianPolish() decomposes into overall, feature and sample effects, so a
  ##   sample sitting far enough below the overall level comes back at or below
  ##   zero. On the log scale combine_features() requires, that is an ordinary
  ##   small value and nothing needs doing about it. Earlier code added a
  ##   per-group constant to lift the whole vector positive, which is a
  ##   per-protein rescaling: harmless to a per-feature coefficient, since it is
  ##   absorbed by the intercept, but it moves the feature's average expression,
  ##   which limma's trended prior and the imputers' p(missing|intensity) both
  ##   read. On rdtc_seer2 precursors it displaced 23% of protein groups:

  f <- function(idxs) {
    MsCoreUtils::medianPolish(state$expression[idxs, , drop=F],
      na.rm=T, maxiter=maxit)
  }
  exprs <- tapply(1:nrow(state$expression), genes, f)
  exprs <- do.call(rbind, exprs)
  
  genes <- feats[[config$gene_id_col]]
  exprs <- exprs[genes, , drop=F]
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
    ## NA is the only indicator of a missing value; see f.zeros_to_na(). Treating
    ##   0 as missing here would drop a feature whose values are all legitimately
    ##   0, which on a log scale is an ordinary measurement:

    i <- apply(x, 1, function(v) !all(is.na(v)))
    x <- x[i, , drop=F]

    ## MsCoreUtils::robustSummary() sheds rank deficient columns of its design by
    ##   dropping those whose coefficient came back exactly 0:
    ##
    ##     repeat { fit <- .lm.fit(X, expression); id <- fit$coefficients != 0
    ##              X <- X[, id, drop=FALSE]; if(all(id)) break }
    ##
    ##   A sample whose values in this gene group are every one of them exactly 0
    ##   has a coefficient of exactly 0 for the honest reason, so its column is
    ##   discarded as though it were rank deficient and the lookup of that sample
    ##   at the end of robustSummary() finds nothing: the summary comes back NA
    ##   for it, silently, and indistinguishably from a sample where the gene was
    ##   never measured. When every column qualifies the design empties out
    ##   altogether and MASS::rlm() stops with a dimnames error instead. Checked
    ##   for here because neither outcome is recoverable downstream, and because
    ##   no upstream filter covers it: the features involved need not be constant,
    ##   so prefilter()'s distinct value rule and filter_features(remove_constant)
    ##   both pass them. Only groups of two or more features are affected, since
    ##   robustSummary() returns a single feature as it stands without fitting.
    ##   Tested against MsCoreUtils 1.12.0:

    if(nrow(x) >= 2) {
      j <- apply(x, 2, function(v) any(!is.na(v)) && all(v[!is.na(v)] %in% 0))
      if(any(j)) {
        f.err("f.combine_features_robust_summary: gene", genes[idxs][1], "has",
          sum(j), "sample(s) whose every measured value in the group is exactly",
          "0, which MsCoreUtils::robustSummary() silently summarizes as NA;",
          "\n", "affected samples:", colnames(x)[j], "\n",
          "features in the group:", rownames(x), "\n",
          "use config$feature_aggregation 'medianPolish', which is unaffected,",
          "or set those values to NA if they are meant to be missing",
          config=config)
      }
    }

    ## a sample estimated below the overall level of the group comes back at or
    ##   below zero, which on the log scale combine_features() requires is an
    ##   ordinary small value; see f.combine_features_median_polish() for why the
    ##   per-group constant that used to lift these is not harmless:

    return(MsCoreUtils::robustSummary(x, na.rm=T))
  }
  exprs <- tapply(1:nrow(state$expression), genes, f)
  nom <- names(exprs)
  exprs <- do.call(rbind, exprs)
  rownames(exprs) <- nom
  
  genes <- feats[[config$gene_id_col]]
  exprs <- exprs[genes, , drop=F]
  
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
#'   \code{rescale=TRUE} (or \code{config$feature_aggregation_scaled=TRUE}) is
#'     refused unless \code{method} is \code{"none"}, where nothing is
#'     aggregated and the setting is reported as ignored. It divided each feature
#'     by its own mean, which is a raw scale operation, and aggregation requires
#'     log scale data (see below), where a per-feature mean near zero explodes
#'     the feature and a negative one flips the sign of its contrasts. The log
#'     scale form of the same idea, subtracting the per-feature mean, is absorbed
#'     exactly by \code{medianPolish()}'s own per-feature effect and so changes
#'     no per-sample effect, which is why the option is removed rather than
#'     corrected.
#'   \code{method="robustSummary"} stops if a gene group of two or more features
#'     contains a sample whose every measured value in that group is exactly
#'     \code{0}. \code{MsCoreUtils::robustSummary()} drops such a sample from its
#'     design and reports it as \code{NA}, which cannot afterwards be told apart
#'     from a sample where the gene was never measured, so it is refused rather
#'     than summarized. This cannot arise from raw input, where
#'     \code{h0testr::initialize()} converts zeros to \code{NA}; it is reachable
#'     only for input that was already transformed and holds genuine zeros.
#'     \code{method="medianPolish"} is unaffected.
#'   For \code{method="medianPolish"}, for each unique gene in 
#'     \code{state$features[, config$gene_id_col]}, the submatrix of 
#'       corresponding peptide signals across all samples is decomposed into: 
#'       \code{pep_exprs == median_column_effect + median_row_effect + overall_median},
#'     Then \code{median_column_effect} is returned as it stands. Values at or
#'       below zero are ordinary on the log scale this function requires, and are
#'       returned unaltered; earlier versions added a per-group constant to lift
#'       each group's vector strictly positive, which left a per-feature
#'       coefficient untouched but moved the feature's average expression, and so
#'       the trended prior of \code{limma::eBayes(trend=TRUE)} and the
#'       \code{p(missing|intensity)} models of several imputers.
#'   Both aggregation methods fit an additive model of an overall level plus
#'     per-feature and per-sample effects, which holds for mass spectrometry
#'     signal only after a log transform, since sample loading and precursor
#'     response are multiplicative. Stops with an error unless
#'     \code{config$is_log_transformed} is \code{TRUE}; run
#'     \code{h0testr::normalize()} first (the default \code{config$run_order}
#'     does), declare an already-transformed input with
#'     \code{config$is_log_transformed <- TRUE}, or use \code{method="none"},
#'     which fits nothing and is exempt.
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
#'     \code{is_log_transformed} \cr \tab Logical scalar; must be \code{TRUE} unless \code{method} is \code{"none"}, since both aggregators fit an additive model. \cr
#'     \code{save_state}      \cr \tab Logical scalar indicating whether to save new state to disk. \cr
#'     \code{dir_out}         \cr \tab Output directory path (scalar character); only needed if \code{save_state == TRUE}. \cr
#'     \code{data_mid_out}    \cr \tab Midfix of expression matrix filename; only needed if \code{save_state == TRUE}.\cr
#'     \code{feature_mid_out} \cr \tab Midfix of feature metadata filename; only needed if \code{save_state == TRUE}.\cr
#'     \code{sample_mid_out}  \cr \tab Midfix of sample metadata filename; only needed if \code{save_state == TRUE}.\cr
#'     \code{suffix_out}      \cr \tab Suffix of output files; only needed if \code{save_state == TRUE}.\cr
#'   }
#' @param method Name (character scalar) of method to use for combining, 
#'   where \code{method \%in\% c("medianPolish", "robustSummary", "none")}. Default: \code{"medianPolish"}.
#' @param rescale Logical scalar indicating whether to rescale peptides prior to aggregation. Must be
#'   \code{FALSE} unless \code{method} is \code{"none"}; see Details. Default: \code{FALSE}.
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
#'
#' ## both aggregators fit an additive model, so they need log scale data;
#' ##   normalize() would do this and set the flag in a full workflow:
#' state <- list(expression=log2(exprs + 1), features=feats, samples=samps)
#' config <- list(feat_id_col="pep", gene_id_col="gene", obs_col="obs", feat_col="pep",
#'   feature_aggregation="medianPolish", is_log_transformed=TRUE, save_state=FALSE)
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
#' ## rescaling before aggregation is refused, since it divides by a per-feature
#' ##   mean, which has no meaning on the log scale aggregation requires:
#' out <- try(h0testr::combine_features(state, config, rescale=TRUE), silent=TRUE)
#' cat(class(out), "\n")
#'
#' ## with method="none":
#' out <- h0testr::combine_features(state, config, method="none")
#' print(out$state)
#' str(out$config)

combine_features <- function(state, config, method=NULL, rescale=FALSE) {
  
  check_config(config)
  
  if(is.null(config$gene_id_col) || !(config$gene_id_col %in% names(state$features))) {
    f.err("combine_features: !(config$gene_id_col %in% names(state$features));",
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

  ## both aggregators fit an additive model: an overall level plus a per-feature
  ##   and a per-sample effect. Mass spectrometry effects are multiplicative
  ##   instead, since loading a sample hot scales every feature in it and a
  ##   precursor's response factor is a gain rather than an offset, so the
  ##   additive form only holds after a log transform. Fit to raw intensities the
  ##   model is misspecified: the per-sample effect is not a shared constant, the
  ##   residual scale is set by the most abundant features in the group so the
  ##   rest are down-weighted for being small rather than for disagreeing, and
  ##   the summary that comes out is an arithmetic contrast where the quantity
  ##   wanted is a ratio. Refused rather than warned about because the resulting
  ##   bias grows with the dynamic range within each group, so it is largest
  ##   exactly where aggregation matters most. method "none" does no fitting and
  ##   is exempt, including when it was selected above:

  if(!(method %in% "none") && !isTRUE(config$is_log_transformed)) {
    f.err("combine_features: aggregation method", method, "fits an additive",
      "model, which requires log scale data, but config$is_log_transformed is",
      "not TRUE;", "\n",
      "  to fix, either (1) run normalize() before combine_features(), which",
      "log transforms and sets the flag -- the default config$run_order",
      "already does this, so check yours if you set it;", "\n",
      "  or (2) if state$expression is already on a log-like scale, say so with",
      "config$is_log_transformed <- TRUE before this step;", "\n",
      "  or (3) set config$feature_aggregation (or the method argument) to",
      "'none' to skip aggregation altogether;", "\n",
      "  config$run_order:", config$run_order,
      "; config$normalization_method:", config$normalization_method,
      config=config)
  }

  ## rescaling divided each feature by its own mean, which is a multiplicative
  ##   operation and so belongs to the raw scale. The check above means it can
  ##   only ever have run on log data, where it has no coherent reading: a
  ##   feature whose mean lands near zero, as happens routinely once RLE has
  ##   centered the data, is divided by a near-zero divisor and explodes, and a
  ##   feature with a negative mean has its sign flipped, reversing the direction
  ##   of every contrast it contributes to. On simulated data with unbalanced
  ##   missingness it moved medianPolish()'s per-sample effects by up to 0.93
  ##   log2 units. Refused rather than corrected because the log scale form of
  ##   the same idea, subtracting the per-feature mean, is absorbed exactly by
  ##   medianPolish()'s own per-feature effect: it leaves the per-sample effects
  ##   identical to 1e-15 and only shifts each group's overall level, which is
  ##   the average expression limma's trended prior reads. So there is nothing
  ##   the option can do that is both meaningful and useful:

  if(isTRUE(rescale)) {
    if(method %in% "none") {
      f.msg("combine_features: rescale is TRUE but method is 'none', so nothing",
        "is aggregated and no rescaling is done", config=config)
      rescale <- FALSE
    } else {
      f.err("combine_features: rescaling before aggregation is no longer",
        "supported; method:", method, "\n",
        "  it divided each feature by its own mean, which is a raw scale",
        "operation, but aggregation requires log scale data, where dividing by a",
        "mean near zero explodes the feature and dividing by a negative mean",
        "flips the sign of its contrasts;", "\n",
        "  done correctly on the log scale it would subtract the per-feature",
        "mean, which medianPolish() already absorbs into its own per-feature",
        "effect, so removing it changes no per-sample effect;", "\n",
        "  to fix, set rescale=FALSE and config$feature_aggregation_scaled to",
        "FALSE;", "\n",
        "  config$feature_aggregation_scaled:",
        config$feature_aggregation_scaled, config=config)
    }
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
  
  prfx <- "combined_features"
  if(!is.null(config$run_order)) {
    i <- config$run_order %in% "combine_features"
    if(any(i)) {
      prfx <- paste0(which(i)[1] + 2, ".combined_features")
    }
  }
  f.save_state(state, config, prefix=prfx)
  
  return(list(state=state, config=config))
}
