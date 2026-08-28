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
#'     of a missing value; see \code{h0testr::init_state()}.
#'   Deletes \code{state$samples[, config$obs_id_col} if not same as
#'     \code{config$sample_id_col}.
#'   Sets \code{config$obs_col} and \code{config$obs_id_col} to 
#'     \code{config$sample_id_col}.
#'   See documentation for \code{h0testr::new_config()} 
#'     for more detailed description of configuration parameters. 
#' @param state List with elements formatted like the list returned by 
#'   \code{read_data()}:
#'   \tabular{ll}{
#'     \code{expression} \tab Numeric matrix with non-negative expression values. \cr
#'     \code{features}   \tab A data.frame with feature meta-data for rows of expression. \cr
#'     \code{samples}    \tab A data.frame with observation meta-data for columns of expression. \cr
#'   } 
#' @param config List with configuration values. Uses the following keys:
#'   \tabular{ll}{
#'     \code{feat_col}        \tab Name of column in \code{state$features} matching \code{rownames(state$expression)}. \cr
#'     \code{obs_id_col}      \tab Name of column in \code{state$samples} matching \code{colnames(state$expression)}. \cr
#'     \code{sample_id_col}   \tab Name of column in \code{state$samples} with unique sample (case) ids. \cr
#'     \code{save_state}      \tab Logical scalar indicating whether to save new state to disk. \cr
#'     \code{dir_out}         \tab Output directory path (scalar character); only needed if \code{save_state == TRUE}. \cr
#'     \code{data_mid_out}    \tab Midfix of expression matrix filename; only needed if \code{save_state == TRUE}.\cr
#'     \code{feature_mid_out} \tab Midfix of feature metadata filename; only needed if \code{save_state == TRUE}.\cr
#'     \code{sample_mid_out}  \tab Midfix of sample metadata filename; only needed if \code{save_state == TRUE}.\cr
#'     \code{suffix_out}      \tab Suffix of output files; only needed if \code{save_state == TRUE}.\cr
#'   }
#' @param fn Function for combining replicated measurements.
#' @return A list (the processed state) with the following elements:
#'   \tabular{ll}{
#'     \code{expression} \tab Numeric matrix with non-negative expression values. \cr
#'     \code{features}   \tab A data.frame with feature meta-data for rows of expression. \cr
#'     \code{samples}    \tab A data.frame with observation meta-data for columns of expression. \cr
#'   } 
#' @examples
#' set.seed(101)
#' exprs <- h0testr::sim1(n_obs=6, n_feats=8)$mat
#' feats <- data.frame(feature_id=rownames(exprs))
#' samps <- data.frame(observation_id=colnames(exprs))
#' samps$sample_id=c("samp1", "samp1", "samp2", "samp2", "samp3", "samp3")
#' state <- list(expression=exprs, features=feats, samples=samps)
#' config <- list(feat_col="feature_id", obs_id_col="observation_id",
#'   sample_id_col="sample_id", save_state=FALSE)
#' out <- h0testr::combine_replicates(state, config)
#' print(state)
#' print(out$state)
#' str(config)
#' str(out$config)
#' @export

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

    f <- function(v, s) {
      out <- tapply(v, s, fn, na.rm=T)
      n_obs <- tapply(v, s, function(x) sum(!is.na(x)))
      out[n_obs %in% 0] <- NA
      return(out)
    }
    sample_ids <- state$samples[[config$sample_id_col]]

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

f.n_feats_col <- function(config) {
  nom <- config$n_feats_col
  if(length(nom) != 1 || is.na(nom) || !nzchar(nom)) return("n_feats")
  return(nom)
}

f.combine_method_col <- function(config) {
  nom <- config$combine_method_col
  if(length(nom) != 1 || is.na(nom) || !nzchar(nom)) return("combine_method")
  return(nom)
}

## helper for the combine_features() aggregators; records which aggregator
##   summarized each gene in the gene level feature table:

f.set_combine_method <- function(feats, route, config, caller) {
  nom <- f.combine_method_col(config)
  if(nom %in% names(feats)) {
    f.msg(caller, ": state$features already has a column named", nom,
      "; leaving it as it is rather than recording this call in it", config=config)
    return(feats)
  }
  if(length(route) %in% 1 && is.null(names(route))) {   ## one route for every gene
    feats[[nom]] <- rep(as.character(route), nrow(feats))
  } else {
    feats[[nom]] <- as.character(route[as.character(feats[[config$gene_id_col]])])
  }
  return(feats)
}

f.gene_ids <- function(feats, config, caller="f.gene_ids") {

  if(is.null(config$gene_id_col) || !(config$gene_id_col %in% names(feats))) {
    f.err(caller, ": !(config$gene_id_col %in% names(feats));",
      "config$gene_id_col:", config$gene_id_col,
      "; names(feats):", names(feats), config=config)
  }

  genes <- as.character(feats[[config$gene_id_col]])
  i <- is.na(genes) | genes %in% ""
  if(any(i)) genes[i] <- paste0("unknown_", feats[[config$feat_id_col]][i])

  return(genes)
}

f.gene_features <- function(feats, config, caller="f.gene_features") {

  genes <- f.gene_ids(feats, config, caller)
  feats[[config$gene_id_col]] <- genes

  keep <- rep(TRUE, ncol(feats))
  for(idx in seq_len(ncol(feats))) {
    if(names(feats)[idx] %in% config$gene_id_col) next
    v <- as.character(feats[[idx]])
    v[is.na(v)] <- "\001na\001"          ## so that NA and NA count as agreeing
    keep[idx] <- all(tapply(v, genes, function(w) length(unique(w)) %in% 1))
  }

  if(any(!keep)) {
    f.msg(caller, ": dropping", sum(!keep), "of", ncol(feats), "feature metadata",
      "columns whose values vary among the features of a gene, and so have no gene",
      "level value:", names(feats)[!keep], config=config)
  }

  out <- feats[!duplicated(genes), keep, drop=F]

  ## the number of features aggregated into each gene:

  nom <- f.n_feats_col(config)
  if(!(nom %in% names(out))) {
    n <- table(genes)
    out[[nom]] <- as.integer(n[out[[config$gene_id_col]]])
  }

  return(list(features=out, genes=genes))
}

## helper for the combine_features() aggregators; both summarize some or all gene
##   groups with MsCoreUtils::medianPolish(), which calls stats::medpolish(), which
##   warns once per group when it runs out of sweeps. 

f.median_polish <- function(x, maxit) {
  converged <- TRUE
  v <- withCallingHandlers(
    MsCoreUtils::medianPolish(x, na.rm=T, maxiter=maxit),
    warning=function(w) {
      if(grepl("did not converge", conditionMessage(w), fixed=T)) {
        converged <<- FALSE
        invokeRestart("muffleWarning")
      }
    })
  return(list(x=v, converged=converged))
}

## helper for the combine_features() aggregators; reports the gene groups whose
##   median polish ran out of sweeps, in the terms established above:

f.msg_median_polish <- function(noms, n_polished, caller, maxit, config) {
  if(length(noms) < 1) return(invisible(NULL))
  f.msg(caller, ":", length(noms), "of", n_polished, "gene group(s) summarized by",
    "MsCoreUtils::medianPolish() did not meet stats::medpolish()'s convergence",
    "criterion within", maxit, "sweeps. That criterion compares the sum of absolute",
    "residuals between successive sweeps, so it describes the sweeps rather than the",
    "fit, and a residual sum that decays by a constant fraction or cycles never",
    "meets it however many sweeps are allowed, even where the summary itself has",
    "stopped moving; raising the sweep limit is usually not the remedy the wording",
    "of medpolish()'s own warning implies. The summaries are returned as they stand,",
    "which for most such groups is the fully iterated value to within 1e-8, but a",
    "minority are genuinely cut short: check these genes if a summary looks wrong.",
    "First affected gene(s):", utils::head(noms, 5), config=config)
  return(invisible(NULL))
}

## helper for combine_features(); uses config$gene_id_col;
##   sets config$feat_col and config$feat_id_col to config$gene_id_col:

f.combine_features_median_polish <- function(state, config, maxit=30) {

  out <- f.gene_features(state$features, config, "f.combine_features_median_polish")
  feats <- out$features
  genes <- out$genes
  config$feat_col <- config$feat_id_col <- config$gene_id_col

  ## medianPolish() decomposes into overall, feature and sample effects; a
  ##   sample far below the overall level can come back <= 0:

  f <- function(idxs) {
    f.median_polish(state$expression[idxs, , drop=F], maxit)
  }
  out <- tapply(1:nrow(state$expression), genes, f)
  nom <- names(out)
  exprs <- do.call(rbind, lapply(out, function(v) v$x))
  rownames(exprs) <- nom

  unconv <- !vapply(out, function(v) v$converged, logical(1))
  f.msg_median_polish(nom[unconv], length(out),
    "f.combine_features_median_polish", maxit, config)

  genes <- feats[[config$gene_id_col]]
  exprs <- exprs[genes, , drop=F]
  feats <- f.set_combine_method(feats, "medianPolish", config,
    "f.combine_features_median_polish")
  state <- list(expression=exprs, features=feats, samples=state$samples)

  return(list(state=state, config=config))
}

## helper for combine_features(); uses config$gene_id_col;
##   sets config$feat_col and config$feat_id_col to config$gene_id_col:

f.combine_features_robust_summary <- function(state, config, maxit=30) {

  out <- f.gene_features(state$features, config, "f.combine_features_robust_summary")
  feats <- out$features
  genes <- out$genes
  config$feat_col <- config$feat_id_col <- config$gene_id_col

  f <- function(idxs) {
    x <- state$expression[idxs, , drop=F]
    ## NA is the only indicator of a missing value:

    i <- apply(x, 1, function(v) !all(is.na(v)))
    x <- x[i, , drop=F]

    ## MsCoreUtils::robustSummary() sheds rank deficient columns of its design by
    ##   dropping those whose coefficient came back exactly 0.

    if(nrow(x) >= 2) {
      j <- apply(x, 2, function(v) any(!is.na(v)) && all(v[!is.na(v)] %in% 0))
      if(any(j)) {
        p <- f.median_polish(x, maxit)
        return(list(x=p$x, route="medianPolish", n_samps=sum(j),
          converged=p$converged))
      }
    }

    ## a sample estimated below the overall level of the group comes back at or
    ##   below zero:

    return(list(x=MsCoreUtils::robustSummary(x),
      route="robustSummary", n_samps=0, converged=TRUE))
  }
  out <- tapply(1:nrow(state$expression), genes, f)
  nom <- names(out)
  route <- stats::setNames(vapply(out, function(v) v$route, character(1)), nom)
  n_samps <- vapply(out, function(v) v$n_samps, numeric(1))
  exprs <- do.call(rbind, lapply(out, function(v) v$x))
  rownames(exprs) <- nom

  i <- route %in% "medianPolish"
  if(any(i)) {
    f.msg("f.combine_features_robust_summary:", sum(i), "of", length(route),
      "gene group(s) contain a sample whose every measured value in the group is",
      "exactly 0,", sum(n_samps), "such (gene, sample) cell(s) in all;",
      "MsCoreUtils::robustSummary() summarizes those samples as NA, silently, so",
      "these groups were summarized with MsCoreUtils::medianPolish() instead;",
      "recorded per gene in", paste0("state$features$", f.combine_method_col(config)),
      "; first affected gene(s):", utils::head(nom[i], 5), config=config)
  }

  unconv <- !vapply(out, function(v) v$converged, logical(1))
  f.msg_median_polish(nom[unconv], sum(i),
    "f.combine_features_robust_summary", maxit, config)

  genes <- feats[[config$gene_id_col]]
  exprs <- exprs[genes, , drop=F]
  feats <- f.set_combine_method(feats, route, config,
    "f.combine_features_robust_summary")

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
#'   A feature whose gene id is missing or blank becomes a gene of its own, named
#'     \code{unknown_} followed by its feature id, rather than all such features
#'     being pooled into one group.
#'   The returned \code{state$features} has one row per gene and keeps only those
#'     columns whose value is the same for every feature of a gene. 
#'   Adds \code{config$n_feats_col} (default \code{"n_feats"}), the number of
#'     features aggregated into each gene. This is the covariate
#'     \code{h0testr::test_deqms()} moderates against. 
#'   \code{rescale=TRUE} is refused unless \code{method} is \code{"none"}, where
#'     nothing is aggregated and the setting is reported as ignored. 
#'   Adds \code{config$combine_method_col} (default \code{"combine_method"}),
#'     naming the aggregator that summarized each gene. 
#'   \code{method="robustSummary"} summarizes a gene group with
#'     \code{MsCoreUtils::medianPolish()} instead if the group holds two or more
#'     features and contains a sample whose every measured value in the group is
#'     exactly \code{0}. 
#'   For \code{method="medianPolish"}, for each unique gene in 
#'     \code{state$features[, config$gene_id_col]}, the submatrix of 
#'       corresponding peptide signals across all samples is decomposed into: 
#'       \code{pep_exprs == median_column_effect + median_row_effect + overall_median},
#'     Then \code{median_column_effect} is returned as it stands. Values at or
#'       below zero are ordinary on the log scale this function requires, and are
#'       returned unaltered.
#'   \code{stats::medpolish()}, which \code{MsCoreUtils::medianPolish()} calls,
#'     warns when it runs out of sweeps. That warning is intercepted and reported
#'     through the log instead, once per call, naming the affected genes.
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
#'     \code{expression} \tab Numeric matrix with non-negative expression values. \cr
#'     \code{features}   \tab A data.frame with feature meta-data for rows of expression. \cr
#'     \code{samples}    \tab A data.frame with observation meta-data for columns of expression. \cr
#'   } 
#' @param config List with configuration values. Uses the following keys:
#'   \tabular{ll}{
#'     \code{obs_col}         \tab Name of column in \code{state$samples} matching \code{colnames(state$expression)}. \cr
#'     \code{feat_col}        \tab Name of column in \code{state$features} matching \code{rownames(state$expression)}. \cr
#'     \code{gene_id_col}     \tab Name of column in \code{state$features} specifying high-level feature for aggregation. \cr
#'     \code{n_feats_col}     \tab Name of new column in \code{state$features} for the number of features aggregated into each gene; preserved if already present. \cr
#'     \code{combine_method_col} \tab Name of new column in \code{state$features} for the aggregator that summarized each gene; preserved if already present. \cr
#'     \code{is_log_transformed} \tab Logical scalar; must be \code{TRUE} unless \code{method} is \code{"none"}, since both aggregators fit an additive model. \cr
#'     \code{save_state}      \tab Logical scalar indicating whether to save new state to disk. \cr
#'     \code{dir_out}         \tab Output directory path (scalar character); only needed if \code{save_state == TRUE}. \cr
#'     \code{data_mid_out}    \tab Midfix of expression matrix filename; only needed if \code{save_state == TRUE}.\cr
#'     \code{feature_mid_out} \tab Midfix of feature metadata filename; only needed if \code{save_state == TRUE}.\cr
#'     \code{sample_mid_out}  \tab Midfix of sample metadata filename; only needed if \code{save_state == TRUE}.\cr
#'     \code{suffix_out}      \tab Suffix of output files; only needed if \code{save_state == TRUE}.\cr
#'   }
#' @param method Name (character scalar) of method to use for combining, 
#'   where \code{method \%in\% c("medianPolish", "robustSummary", "none")}. Default: \code{"medianPolish"}.
#' @param rescale Logical scalar indicating whether to rescale peptides prior to aggregation. Must be
#'   \code{FALSE} unless \code{method} is \code{"none"}; see Details. Default: \code{FALSE}.
#' @return A list (the processed state) with the following elements:
#'   \tabular{ll}{
#'     \code{expression} \tab Numeric matrix with non-negative expression values. \cr
#'     \code{features}   \tab A data.frame with feature meta-data for rows of expression; one row per gene unless \code{method} is \code{"none"}, and carrying \code{config$n_feats_col} and \code{config$combine_method_col}. \cr
#'     \code{samples}    \tab A data.frame with observation meta-data for columns of expression. \cr
#'   }
#' @examples
#' ## set up data for examples; p_drop=0 so that every gene keeps all three of its
#' ##   peptides and so has something to aggregate:
#' set.seed(101)
#' samps <- h0testr::sim_samples(factors=list(grp=c("ctl", "trt")), n_per_cell=3)
#' sim <- h0testr::sim_design(samps, frm=~grp, test_term="grp", n_genes=5, n_genes_signif=1,
#'   peps_per_gene=3, p_drop=0)
#'
#' ## both aggregators fit an additive model, so they need log scale data;
#' ##   normalize() would do this and set the flag in a full workflow:
#' state <- sim$state
#' state$expression <- log2(state$expression + 1)
#' config <- sim$config
#' config$is_log_transformed <- TRUE
#' config$feature_aggregation <- "medianPolish"
#' config$feat_col <- config$feat_id_col   ## init_state() sets these two
#' config$obs_col <- config$obs_id_col
#' rm(samps, sim)
#' print(state)
#' str(config)
#'
#' ## combine peps using config$feature_aggregation ("medianPolish"):
#' out <- h0testr::combine_features(state, config)
#' print(out$state)
#' str(out$config)
#' 
#' ## combine peps, overriding config$feature_aggregation; state$features records
#' ##   which aggregator summarized each gene:
#' out <- h0testr::combine_features(state, config, method="robustSummary")
#' print(out$state)
#' str(out$config)
#'
#' ## a gene group holding a sample whose every measured value is exactly 0 is
#' ##   summarized with medianPolish() instead, robustSummary() reporting such a
#' ##   sample as NA; the route is recorded per gene:
#' state2 <- state
#' gene1 <- state2$features$gene_id[1]
#' state2$expression[state2$features$gene_id %in% gene1, 1] <- 0
#' out <- h0testr::combine_features(state2, config, method="robustSummary")
#' print(out$state$features)
#' print(out$state$expression[, 1])
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
#' @export

combine_features <- function(state, config, method=NULL, rescale=FALSE) {
  
  check_config(config)
  
  if(is.null(config$gene_id_col) || !(config$gene_id_col %in% names(state$features))) {
    f.err("combine_features: !(config$gene_id_col %in% names(state$features));",
      "config$gene_id_col:", config$gene_id_col, 
      "; names(state$features):", names(state$features), config=config)
  }
  
  if(is.null(method) || method %in% "") method <- config$feature_aggregation
  if(is.null(method) || method %in% "") method <- "medianPolish"
  if(is.null(rescale)) rescale <- FALSE
  
  if(config$gene_id_col %in% config$feat_col) {
    f.msg("combine_features: config$gene_id_col %in% config$feat_col; ", 
      "returning unchanged state and updated config.", config=config)
    config$feat_col <- config$gene_id_col <- config$feat_id_col
    method <- "none"
  }

  ## both aggregators fit an additive model: an overall level plus a per-feature
  ##   and a per-sample effect. Mass spectrometry effects are multiplicative
  ##   instead, so the additive form only holds after a log transform:

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
  ##   operation and so belongs to the raw scale:

  if(isTRUE(rescale)) {
    if(method %in% "none") {
      f.msg("combine_features: rescale is TRUE but method is 'none', so nothing",
        "is aggregated and no rescaling is done", config=config)
      rescale <- FALSE
    } else {
      f.err("combine_features: rescaling before aggregation is not",
        "supported; method:", method, "\n",
        " divided each feature by its mean, which is a raw scale",
        "operation, but aggregation requires log scale data\n",
        "  to fix, drop the argument or set rescale=FALSE; the config key that used",
        "to set it, config$feature_aggregation_scaled, has been removed",
        config=config)
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
