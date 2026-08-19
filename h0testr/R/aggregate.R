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

## The gene level form of a feature metadata table: one row per value of
##   config$gene_id_col, carrying only what is actually a property of the gene. Used by
##   both aggregators below and by test() for the two methods that aggregate
##   internally, so that a gene level result table reports the same metadata whichever
##   route produced it.
##   A column whose value differs among the features of a gene describes the feature
##   and not the gene, so it has no reading on a gene row; carried forward it silently
##   reports whichever feature happened to sort first, which is what this used to do.
##   Dropped rather than collapsed into a list or a joined string, since every
##   downstream step expects one value per gene. Nothing is lost that matters:
##   config$feat_id_col always varies within a gene and so goes automatically, and the
##   per-feature statistics add_filter_stats() writes are recomputed by filter().
##   Returns the table and the per-input-row gene ids, since the caller needs the
##   latter to group the expression matrix and the two must agree about which features
##   fell into the unknown_* genes below:

## The gene id of each row of a feature metadata table, with a feature that has no
##   gene assignment made its own gene rather than all such features being pooled into
##   one meaningless group. Separate from f.gene_features() because f.gene_counts()
##   needs these ids and nothing else, and the two must agree about them:

## The name of the column holding the number of features behind each gene. Defaulted
##   rather than required, since a hand-built config carrying only the keys a function
##   documents is supported usage throughout the package, and the aggregators would
##   otherwise fail on one that predates this column:

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
##   summarized each gene in the gene level feature table. An existing column of
##   that name is preserved rather than overwritten, for the same reason
##   f.gene_features() preserves the feature counts: a table that has been
##   aggregated once already carries the record of how its values were combined
##   from features, and a second pass has one feature per gene and so nothing to
##   say about it:

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

  ## the number of features aggregated into each gene, which is a property of the gene
  ##   and is what DEqMS moderates on; see f.gene_counts(). Preserved rather than
  ##   recomputed when already present, since a table that has been here before is
  ##   already one row per gene and would recount every gene as 1, discarding the
  ##   counts of the features that actually went into it:

  nom <- f.n_feats_col(config)
  if(!(nom %in% names(out))) {
    n <- table(genes)
    out[[nom]] <- as.integer(n[out[[config$gene_id_col]]])
  }

  return(list(features=out, genes=genes))
}

## helper for the combine_features() aggregators; both summarize some or all gene
##   groups with MsCoreUtils::medianPolish(), which calls stats::medpolish(), which
##   warns once per group when it runs out of sweeps. Called here so that warning
##   can be reported by the caller, once per call and in terms of what it means,
##   rather than escaping one gene at a time.
##
##   medpolish() stops when the sum of absolute residuals changes between
##   successive sweeps by less than eps=0.01 of itself, which is a statement about
##   the sweeps and not about the fit, and it is routinely unreachable: a residual
##   sum that decays by a constant fraction each sweep changes by that same
##   fraction forever, and one that cycles with period 2 never settles either.
##   Neither is pathological. Measured over the rdtc_seer2 precursors, 3913 protein
##   groups of two or more: 84 groups do not converge in 30 sweeps, and raising
##   maxiter to 1000 leaves 82 of those 84 still unconverged while tripling the
##   time spent aggregating (39 s to 127 s) and changing a summary by more than
##   0.01 for a single group. Over 4512 random groups the split was exact: where
##   the residual sum cycled, the returned summary was a fixed point, identical at
##   30 and at 1000 sweeps for all 87 such groups; where it was still creeping, 60
##   of 60 moved, one by 8.7 on the log2 scale. So more sweeps is neither the
##   remedy the warning implies nor free, and silence would hide the creeping
##   minority that really are cut short. Hence a message rather than either.
##
##   Only the non-convergence warning is intercepted; any other warning from the
##   fit is left to propagate:

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
##   median polish ran out of sweeps, in the terms established above. Both
##   aggregators say the same thing, so both call this:

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
    ##   altogether and MASS::rlm() stops with a dimnames error instead. Detected
    ##   here because neither outcome is recoverable downstream, and because no
    ##   upstream filter covers it: the features involved need not be constant,
    ##   so prefilter()'s distinct value rule and filter_features(remove_constant)
    ##   both pass them. Only groups of two or more features are affected, since
    ##   robustSummary() returns a single feature as it stands without fitting.
    ##
    ##   Such a group is summarized by MsCoreUtils::medianPolish() instead, which
    ##   fits the same additive model by medians and so has no design to shed and
    ##   no rank to lose. The alternatives were refusing the run, which was the
    ##   previous behaviour and stops a whole analysis over a handful of genes,
    ##   and writing NA for the affected samples, which cannot be right: the
    ##   failing unit is a (gene, sample) cell, so removing features can never
    ##   resolve it, and a value of 0 on the log scale is an ordinary measurement
    ##   here (see f.zeros_to_na()), not a missing one. The route taken is
    ##   recorded per gene rather than left to the log, since the two aggregators
    ##   do not return the same numbers.
    ##
    ##   The rule is deliberately conservative: with an odd number of features the
    ##   offending coefficient sometimes comes back at 1e-16 rather than exactly
    ##   0, survives the != 0 test and is summarized correctly. Those genes are
    ##   routed to medianPolish() as well, which costs a small difference in the
    ##   summary rather than a wrong value. Tested against MsCoreUtils 1.12.0:

    if(nrow(x) >= 2) {
      j <- apply(x, 2, function(v) any(!is.na(v)) && all(v[!is.na(v)] %in% 0))
      if(any(j)) {
        p <- f.median_polish(x, maxit)
        return(list(x=p$x, route="medianPolish", n_samps=sum(j),
          converged=p$converged))
      }
    }

    ## a sample estimated below the overall level of the group comes back at or
    ##   below zero, which on the log scale combine_features() requires is an
    ##   ordinary small value; see f.combine_features_median_polish() for why the
    ##   per-group constant that used to lift these is not harmless.
    ##
    ##   No na.rm here: robustSummary()'s formals are (x, ...), so an na.rm passed
    ##   in falls through to MASS::rlm(), which has no such argument either and
    ##   warns "some of ... do not match" on every group summarized. It drops NA
    ##   itself, by the !is.na(x) mask it builds before fitting, so the values are
    ##   identical either way; only the warning goes. Checked against MsCoreUtils
    ##   1.12.0 and MASS 7.3-60:

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
#'     columns whose value is the same for every feature of a gene. A column that
#'     varies within a gene describes the feature rather than the gene, so it has
#'     no gene level value; carried forward it would report whichever feature
#'     happened to come first. The dropped columns are named in the log. This
#'     always removes \code{config$feat_id_col}, and removes the per-feature
#'     statistics \code{h0testr::add_filter_stats()} writes, which
#'     \code{h0testr::filter()} recomputes for the aggregated features.
#'   Adds \code{config$n_feats_col} (default \code{"n_feats"}), the number of
#'     features aggregated into each gene. This is the covariate
#'     \code{h0testr::test_deqms()} moderates against, so recording it here lets
#'     that method run on an already aggregated state. An existing column of that
#'     name is preserved rather than recomputed, since a table that has been
#'     aggregated once has one feature per gene and recounting would report
#'     \code{1} for every gene.
#'   \code{rescale=TRUE} is refused unless \code{method} is \code{"none"}, where
#'     nothing is aggregated and the setting is reported as ignored. The config
#'     key that used to set it, \code{config$feature_aggregation_scaled}, has been
#'     removed, having only ever been refused. It divided each feature
#'     by its own mean, which is a raw scale operation, and aggregation requires
#'     log scale data (see below), where a per-feature mean near zero explodes
#'     the feature and a negative one flips the sign of its contrasts. The log
#'     scale form of the same idea, subtracting the per-feature mean, is absorbed
#'     exactly by \code{medianPolish()}'s own per-feature effect and so changes
#'     no per-sample effect, which is why the option is removed rather than
#'     corrected.
#'   Adds \code{config$combine_method_col} (default \code{"combine_method"}),
#'     naming the aggregator that summarized each gene. For
#'     \code{method="medianPolish"} that is \code{"medianPolish"} throughout; for
#'     \code{method="robustSummary"} it is \code{"robustSummary"} except for the
#'     gene groups described next. Nothing is added by \code{method="none"},
#'     which aggregates nothing. An existing column of that name is preserved
#'     rather than overwritten, for the same reason the feature counts are.
#'   \code{method="robustSummary"} summarizes a gene group with
#'     \code{MsCoreUtils::medianPolish()} instead if the group holds two or more
#'     features and contains a sample whose every measured value in the group is
#'     exactly \code{0}. \code{MsCoreUtils::robustSummary()} sheds such a sample
#'     from its design, since its coefficient comes back exactly \code{0} as
#'     though the design were rank deficient, and reports the sample as
#'     \code{NA}, which cannot afterwards be told apart from a sample where the
#'     gene was never measured; when every sample qualifies, the design empties
#'     and \code{MASS::rlm()} stops instead. \code{medianPolish()} fits the same
#'     additive model by medians, so it has no design to shed, and the affected
#'     genes are summarized rather than lost. Earlier versions refused the whole
#'     run here. Writing \code{NA} for those samples was the other option and is
#'     wrong: the failing unit is a (gene, sample) cell, so no choice of features
#'     to keep can resolve it, and \code{0} on the log scale this function
#'     requires is an ordinary measurement, not a missing one. The count of
#'     affected genes and cells is logged and the route recorded per gene, since
#'     the two aggregators do not return the same numbers. The rule is
#'     deliberately conservative: the offending coefficient sometimes comes back
#'     at \code{1e-16} rather than exactly \code{0} and is then summarized
#'     correctly, and such a group is routed to \code{medianPolish()} anyway.
#'     None of this can arise from raw input, where
#'     \code{h0testr::initialize()} converts zeros to \code{NA}; it is reachable
#'     only for input that was already transformed and holds genuine zeros.
#'     \code{method="medianPolish"} is unaffected throughout.
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
#'   \code{stats::medpolish()}, which \code{MsCoreUtils::medianPolish()} calls,
#'     warns when it runs out of sweeps. That warning is intercepted and reported
#'     through the log instead, once per call, naming the affected genes. Its own
#'     wording implies that more sweeps would help, and usually they would not:
#'     the convergence criterion compares the sum of absolute residuals between
#'     successive sweeps, so a residual sum that decays by a constant fraction, or
#'     cycles with period 2, never satisfies it however long it runs, even once the
#'     summary itself has stopped moving. On 3913 rdtc_seer2 protein groups, 84 did
#'     not converge in the 30 sweeps used here, and allowing 1000 left 82 of them
#'     unconverged while tripling the time spent aggregating and changing one
#'     group's summary by more than \code{0.01}. A minority are genuinely cut
#'     short, which is why this is reported rather than suppressed. Any other
#'     warning from the fit propagates as usual.
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
#'     \code{n_feats_col}     \cr \tab Name of new column in \code{state$features} for the number of features aggregated into each gene; preserved if already present. \cr
#'     \code{combine_method_col} \cr \tab Name of new column in \code{state$features} for the aggregator that summarized each gene; preserved if already present. \cr
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
#'     \code{features}   \cr \tab A data.frame with feature meta-data for rows of expression; one row per gene unless \code{method} is \code{"none"}, and carrying \code{config$n_feats_col} and \code{config$combine_method_col}. \cr
#'     \code{samples}    \cr \tab A data.frame with observation meta-data for columns of expression. \cr
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
#' config$feat_col <- config$feat_id_col   ## initialize() sets these two
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

combine_features <- function(state, config, method=NULL, rescale=FALSE) {
  
  check_config(config)
  
  if(is.null(config$gene_id_col) || !(config$gene_id_col %in% names(state$features))) {
    f.err("combine_features: !(config$gene_id_col %in% names(state$features));",
      "config$gene_id_col:", config$gene_id_col, 
      "; names(state$features):", names(state$features), config=config)
  }
  
  if(is.null(method) || method %in% "") method <- config$feature_aggregation
  if(is.null(method) || method %in% "") method <- "medianPolish"
  ## no config key sets this any more: config$feature_aggregation_scaled was removed,
  ##   having only ever been refused, so rescaling is reachable only by asking for it
  ##   here, and only to be told why it is gone:

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
