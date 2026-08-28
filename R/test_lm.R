## Helper for test_lm(). Numeric vector y with one feature's expression across
##   all observations, design matrix X for config$frm over those same
##   observations, integer vector cols_test indexing columns of X that support
##   test of config$test_term, integer vector cols_report indexing 
##   columns whose coefficients are to be returned; returns a named numeric
##   vector with $pval and $stat of an F test of the two models, followed
##   by the coefficient estimates of the full model for cols_report:

f.test_lm_feat <- function(y, X, X_red, cols_report) {

  i <- !is.na(y)
  yy <- y[i]
  xf <- X[i, , drop=F]
  xr <- X_red[i, , drop=F]

  noms <- colnames(X)[cols_report]
  out <- c(pval=as.numeric(NA), stat=as.numeric(NA),
    stats::setNames(rep(as.numeric(NA), length(noms)), noms))

  rank_full <- f.design_rank(xf)

  if(rank_full - f.design_rank(xr) < 1 || length(yy) - rank_full < 1) return(out)

  fit_full <- stats::lm(yy ~ xf + 0)

  if(ncol(xr) %in% 0) {
    fit_reduced <- stats::lm(yy ~ 0)
  } else {
    fit_reduced <- stats::lm(yy ~ xr + 0)
  }

  ## exact F test of two nested least squares fits, rather than a likelihood
  ##   ratio chi-square lmtest::lrtest() reported. The chi-square is the large-sample
  ##   limit of this test:

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
#'   Suitable for use with missing values without imputation. Does not use moderated 
#'     standard error estimates, so only rely upon when there is plenty of replication 
#'     (at least 5 observations per condition, preferably more).
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
#'   Defaults to \code{"BH"}.
#' @return
#'   List with a \code{$hits} data.frame of results and a \code{$fit} that is
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
#'     \code{h0testr::test_h0()} reports as \code{logfc} when one design matrix column carries the
#'     test; see \code{h0testr::test_h0()} for what it reports for a joint test or a contrast.
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
#' out <- h0testr::init_state(state, config, minimal=TRUE)
#' out$state <- h0testr::filter_features_by_formula(out$state, out$config)
#' 
#' tbl <- h0testr::test_lm(out$state, out$config)
#' print(tbl)
#' @export

test_lm <- function(state, config, fdr.method="BH") {

  check_config(config)
  f.check_state(state, config)

  ## levels of each factor covariate ordered with declared reference level first, 
  ##   before any design is built:

  state <- f.relevel_state_covariates(state, config, caller="test_lm")

  ## full design and columns of it that carry test; shared with
  ##   filter_features_by_estimability():

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
  ##   the coefficients after them; taken by position rather than by name:

  pvals <- hits[, 1, drop=T]
  stats_f <- hits[, 2, drop=T]
  coefs <- hits[, -(1:2), drop=F]
  n <- apply(coefs, 2, function(v) sum(!is.na(v)))
  coefs <- coefs[, n > 0, drop=F]

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

