## Helper for test_lm(). Numeric vector y with one feature's expression across
##   all observations, design matrix X for config$frm over those same
##   observations, integer vector cols_test indexing the columns of X that carry
##   the test of config$test_term, integer vector cols_report indexing the
##   columns whose coefficients are to be returned; returns a named numeric
##   vector with the $pval and the $stat of an F test of the two models, followed
##   by the coefficient estimates of the full model for cols_report.
##   Both models are fitted from columns of one design matrix rather than from
##   formula text. Re-deriving the reduced model as a formula lets
##   stats::model.matrix() re-code the remaining factors to full rank, restoring
##   the span that was meant to be removed and leaving nothing to test: dropping
##   the intercept from ~grp gives ~0 + grp, whose three indicator columns span
##   the same space as the intercept plus two contrasts. Subsetting columns
##   instead keeps the reduced model a strict sub-model of the full one, and
##   makes this test agree with filter_features_by_estimability(), which screens
##   features by the rank of these same two matrices.
##   X_red is the reduced design f.design_test_cols() derived, taken as a matrix
##   rather than as columns to drop, so that a config$contrast run needs no separate
##   code here: the model constrained so that the contrast is zero is nested in the
##   full model in exactly the same way, and the same F test of the two is the test
##   of the contrast:

f.test_lm_feat <- function(y, X, X_red, cols_report) {

  i <- !is.na(y)
  yy <- y[i]
  xf <- X[i, , drop=F]
  xr <- X_red[i, , drop=F]

  noms <- colnames(X)[cols_report]
  out <- c(pval=as.numeric(NA), stat=as.numeric(NA),
    stats::setNames(rep(as.numeric(NA), length(noms)), noms))

  ## the two ranks filter_features_by_estimability() screens on, from the same helper
  ##   and restricted the same way to the observations where this feature was measured:
  ##   its df_test is rank_full - rank_red and its df_resid is n - rank_full, and both
  ##   have to be at least 1 for there to be a test. Neither is guaranteed here, since
  ##   that filter is optional and this function is reached without it. Checked ahead of
  ##   the fits rather than after them because the failures are not all recoverable: for
  ##   a feature measured in no observation at all stats::lm() does not return a fit, it
  ##   errors with "0 (non-NA) cases", and that error propagates out of the apply() in
  ##   test_lm() and takes every other feature's result with it. A saturated full model
  ##   does return, but with no residual variance to test against, so stats::anova()
  ##   reports the F and its p-value as NaN. Returning NA for such a feature says what
  ##   happened instead; see test_lm() for what becomes of it:

  rank_full <- f.design_rank(xf)

  if(rank_full - f.design_rank(xr) < 1 || length(yy) - rank_full < 1) return(out)

  fit_full <- stats::lm(yy ~ xf + 0)

  ## a reduced model with no columns left is the null model with no parameters,
  ##   which lm() will not express as a matrix with zero columns:

  if(ncol(xr) %in% 0) {
    fit_reduced <- stats::lm(yy ~ 0)
  } else {
    fit_reduced <- stats::lm(yy ~ xr + 0)
  }

  ## the exact F test of the two nested least squares fits, rather than the likelihood
  ##   ratio chi-square lmtest::lrtest() reported. The chi-square is the large-sample
  ##   limit of this test, and it reaches that limit by treating the residual variance
  ##   as known rather than estimated; with the handful of observations per condition
  ##   that a proteomics design usually has, that makes it markedly anti-conservative.
  ##   On a 12 observation two-group fixture the two disagreed by an order of magnitude
  ##   on the same data, p = 7.8e-04 by F against p = 6.9e-05 by chi-square, and in
  ##   that direction for every feature. This F is also the zero prior degrees of freedom
  ##   limit of the moderated F the limma family here reports, so it makes "lm" the
  ##   unmoderated member of the same family rather than a different test.
  ##   stats::anova() takes the models smaller first, so the comparison is in row 2:

  tbl <- stats::anova(fit_reduced, fit_full)

  out[["pval"]] <- tbl[["Pr(>F)"]][2]
  out[["stat"]] <- tbl[["F"]][2]

  ## lm() prefixes coefficient names with the name of the matrix it was given:

  coefs <- stats::coef(fit_full)
  names(coefs) <- sub("^xf", "", names(coefs))

  out[noms] <- coefs[noms]          ## missing (aliased) coefficients give NA

  return(out)
}

#' Hypothesis testing using the \code{stats::lm()} function
#' @description
#'   Tests for differential expression by fitting full and reduced linear
#'     models, then comparing them with an F test.
#' @details
#'   This test is suitable for use with missing values without imputation. The
#'     test does not use moderated standard error estimates, so should only be 
#'     relied upon when there is plenty of replication (at least 5 observations 
#'     per condition, preferably more).
#'   The design matrix for \code{config$frm} is built once, with
#'     \code{stats::model.matrix()}. For each feature, \code{stats::lm()} fits a
#'     full model to that matrix and a reduced model to the same matrix with the
#'     columns of \code{config$test_term} removed, both restricted to the
#'     observations where the feature was measured. Removing columns rather than
#'     rebuilding a reduced formula matters: re-deriving the reduced model from
#'     formula text lets \code{stats::model.matrix()} re-code the remaining
#'     factors to full rank and restore the span that was to be removed, leaving
#'     nothing to test. The columns removed are those of \code{config$test_term}
#'     together with those of every higher-order term containing it, so testing
#'     a variable involved in an interaction tests the interaction too.
#'   The two models are compared with the exact F test of \code{stats::anova()},
#'     rather than with the likelihood ratio chi-square of
#'     \code{lmtest::lrtest()} that earlier versions used. The chi-square is the
#'     large-sample limit of the same test, reached by treating the residual
#'     variance as known rather than estimated, and is anti-conservative at the
#'     replication these designs usually have. Raw p-values are adjusted for
#'     multiple testing using \code{stats::p.adjust()}.
#'   A feature measured in too few observations to support the test yields
#'     \code{NA} for its p-value, its adjusted p-value and its statistic, which
#'     earlier versions reported as a p-value of 1. Too few means either that the
#'     term under test has no estimable degrees of freedom left once the
#'     unmeasured observations are dropped, or that the full model is left with
#'     no residual degrees of freedom; these are the \code{df_test} and
#'     \code{df_resid} of \code{h0testr::filter_features_by_estimability()},
#'     computed here from the same two matrices by the same helper. Use that
#'     function beforehand to drop such features instead. Since
#'     \code{stats::p.adjust()} takes its \code{n} from the p-values that are not
#'     \code{NA}, an untested feature does not count against the features that
#'     were tested.
#'   Any covariate type that \code{stats::model.matrix()} accepts is supported,
#'     categorical or continuous, in any combination, as are interactions and
#'     formulas without an intercept.
#' @param state List with elements like those returned by \code{read_data()}:
#'   \tabular{ll}{
#'     \code{expression} \cr \tab Numeric matrix with non-negative expression values. \cr
#'     \code{features}   \cr \tab A data.frame with feature meta-data for rows of expression. \cr
#'     \code{samples}    \cr \tab A data.frame with observation meta-data for columns of expression. \cr
#'   } 
#' @param config List with configuration values. Uses the following keys:
#'   \tabular{ll}{
#'     \code{feat_col}       \cr \tab Name of column in \code{state$features} corresponding to \code{rownames(state$expression)}. \cr
#'     \code{frm}            \cr \tab Formula (formula object) to be fit \cr
#'     \code{test_term}      \cr \tab Term (character scalar) to be tested for non-zero coefficient. \cr
#'     \code{contrast}      \cr \tab Weighted sum (character scalar) of coefficients of \code{config$frm} to test instead of \code{config$test_term}; "" for none. \cr
#'   }
#' @param fdr.method Character scalar specifying method to use for multiple
#'   testing adjustment of p-values. See \code{stats::p.adjust.methods()} for
#'   latest list of valid choices. Currently one of:
#'   \code{c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY")}.
#'   Defaults to \code{"BH"}, which is what every other test method here uses;
#'   earlier versions defaulted this one method to the considerably more
#'   conservative \code{"BY"}, so the same data gave \code{"lm"} a larger adjusted
#'   p-value than the other methods for reasons that had nothing to do with the
#'   test. Pass \code{"BY"} to control the false discovery rate under arbitrary
#'   dependence between features. Reachable only by calling this function directly:
#'   \code{h0testr::test()} passes no value for it and no \code{config} key carries one,
#'   so a workflow run gets \code{"BH"}, which is the point of that being the default.
#' @return
#'   A list with a \code{$hits} data.frame of results and a \code{$fit} that is
#'     always \code{NULL}, since a separate model is fitted to every feature.
#'     Columns of \code{$hits} are \code{c("feature", "p.adj", "pval", "stat")},
#'     where \code{stat} is the F statistic of the model comparison, then the
#'     full-model coefficient estimates for the columns under test, preceded by
#'     \code{Intercept} when \code{config$frm} has one, then the feature metadata
#'     columns from \code{state$features}. Coefficient columns that are
#'     \code{NA} for every feature are dropped, and \code{:} in an interaction
#'     coefficient name becomes \code{.}.
#'   The coefficient columns are effect sizes on the scale of \code{state$expression}: for a
#'     two level factor, the difference between its levels, so a log fold change when the input
#'     is log transformed; for a \strong{continuous} covariate, the change \strong{per unit} of
#'     it, whose size depends on the units the covariate is recorded in. This is what
#'     \code{h0testr::test()} reports as \code{logfc} when one design matrix column carries the
#'     test; see \code{h0testr::test()} for what it reports for a joint test or a contrast.
#' @examples
#' set.seed(101)
#' samps <- h0testr::sim_samples(factors=list(condition=c("placebo", "drug")),
#'   n_per_cell=6)
#' sim <- h0testr::sim_design(samps, frm=~condition, test_term="condition", n_genes=25,
#'   n_genes_signif=5, effects=2)
#' state <- sim$state
#' state$expression <- log2(state$expression + 1)
#' 
#' config <- sim$config   ## frm, test_term, the id columns and reference_levels, set
#'                        ##   to match what was simulated; save_state is FALSE
#' config$is_log_transformed <- TRUE   ## normalize() would; logged just above
#' rm(samps, sim)
#'
#' ## set up and check covariates and parameters:
#' out <- h0testr::initialize(state, config, minimal=TRUE)
#' out$state <- h0testr::filter_features_by_formula(out$state, out$config)
#' 
#' tbl <- h0testr::test_lm(out$state, out$config)
#' print(tbl)

test_lm <- function(state, config, fdr.method="BH") {

  check_config(config)
  f.check_state(state, config)

  ## the levels of each factor covariate ordered with the declared reference level
  ##   first, before any design is built: f.design_X() hands state$samples to
  ##   stats::model.matrix() as they are, and a character column there is levelled by
  ##   sorting, so without this a direct caller gets coefficients named for a
  ##   reference level it did not ask for, and the sign of the effect flipped, with
  ##   nothing said about it. A no-op for a run that came through test(), which has
  ##   already done exactly this; kept because this function is exported and a direct
  ##   caller of it reaches no other place that does. Same call, and for the same
  ##   reason, as in test_proda() and test_prolfqua():

  state <- f.relevel_state_covariates(state, config, caller="test_lm")

  ## the full design and the columns of it that carry the test; shared with
  ##   filter_features_by_estimability(), so the models compared here are the
  ##   ones that were screened for estimability. Errors if config$test_term does
  ##   not fit config$frm, or if the test it names is vacuous:

  design <- f.design_test_cols(state, config)
  X <- design$X
  cols_test <- design$cols_test

  ## coefficients reported: those under test, plus the intercept when the design
  ##   has one, which is what earlier versions of this function reported:

  asgn <- attr(X, "assign")
  cols_report <- sort(unique(c(which(asgn %in% 0), cols_test)))

  f.msg("test_lm:", f.test_label(design, config), "; fdr.method:", fdr.method,
    "; design columns:", ncol(X), "; test columns:", length(cols_test),
    "; df:", design$df_intend, config=config)

  hits <- t(apply(state$expression, 1, f.test_lm_feat, X, design$X_red,
    cols_report))

  ## f.test_lm_feat() returns the p-value and the statistic first, in that order, and
  ##   the coefficients after them; taken by position rather than by name so that a
  ##   covariate whose design column happens to be called pval or stat cannot be
  ##   mistaken for one of them:

  pvals <- hits[, 1, drop=T]
  stats_f <- hits[, 2, drop=T]
  coefs <- hits[, -(1:2), drop=F]
  n <- apply(coefs, 2, function(v) sum(!is.na(v)))
  coefs <- coefs[, n > 0, drop=F]

  ## a feature the design could not support keeps its NA p-value, where it used to be
  ##   reported as 1.0. Of all the values that could stand in for a test that was never
  ##   run, 1.0 is the one that is certainly wrong: it is a claim about the data, and it
  ##   is indistinguishable in the result table from a feature that really was tested
  ##   and really showed nothing. It also inflated the count that stats::p.adjust()
  ##   divides by, which takes n from the p-values that are not NA, so every other
  ##   feature's adjusted p-value was penalized for tests that never happened. The NA
  ##   stays in place instead, and carries through p.adjust() to p.adj and into the
  ##   standardized table, where it reads as "not tested". Use
  ##   h0testr::filter_features_by_estimability() beforehand to drop such features:

  fdrs <- stats::p.adjust(pvals, method=fdr.method)

  hits <- data.frame(feature=rownames(hits), p.adj=fdrs, pval=pvals, stat=stats_f,
    coefs)
  rownames(hits) <- NULL
  i <- names(hits) %in% "X.Intercept."
  if(any(i)) names(hits)[i] <- "Intercept"
  
  tmp <- state$features
  rownames(tmp) <- tmp[[config$feat_col]]
  hits <- cbind(hits, tmp[hits$feature, , drop=F])
  rownames(hits) <- NULL
  
  return(list(hits=hits, fit=NULL))
}

