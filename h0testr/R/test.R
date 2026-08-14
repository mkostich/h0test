## Helper for f.normalize_terms(), which is a helper for f.design_test_cols(),
##   which serves test() and filter_features_by_estimability(). Converts formula
##   config$frm to character, makes intercept explicit (either '0' or '1'), sorts
##   variables in interaction terms (so e.g. 'sex:age' becomes 'age:sex'), then
##   returns formula representation as tokenized character vector. So would take
##   formula e.g. ~age + strain + strain:age, and return:
##   c("1", "age", "strain", "age:strain").

f.formula2terms <- function(config) {

  frm <- config$frm

  if(is.null(frm) || all(as.character(frm) %in% "")) {
    f.err("f.formula2terms: config$frm empty or undefined.", config=config)
  }

  parsed <- f.parse_frm(frm, config)

  ## make intercept explicit, so that it can be dropped or kept downstream;
  ##   '*', '^', and '/' expansion, as well as intercept removal, are handled
  ##   by stats::terms() within f.parse_frm(), so e.g. ~x1*x2 and
  ##   ~x1+x2+x1:x2 both give c("1", "x1", "x2", "x1:x2"):

  if(parsed$intercept %in% 1) {
    terms <- c("1", parsed$labels)
  } else {
    terms <- c("0", parsed$labels)
  }

  return(terms)
}

## Helper for f.normalize_terms(). Given a parsed formula (see f.parse_frm())
##   and a canonical config$test_term, returns the character vector of term
##   labels to be dropped from the full model to form the reduced model.
##   A test_term naming a variable drops that variable's term along with every
##   higher-order term containing it (so testing 'x1' in ~x1*x2 is a 2 df test
##   of 'x1' and 'x1:x2'); this keeps the reduced model hierarchical, so the
##   test does not depend on the contrast coding or on which level of x2 is the
##   reference. A test_term naming an interaction term drops just that term,
##   and is an error if some higher-order term in the formula contains it:

f.test_term_drops <- function(parsed, test_term, config) {

  if(test_term %in% "1") {
    if(parsed$intercept %in% 0) {
      f.err("f.test_term_drops: test_term is intercept ('1'), but frm has no",
        "intercept; frm:", parsed$frm, config=config)
    }
    return("1")
  }

  vars <- rownames(parsed$factors)

  ## test_term names a variable: drop every term it takes part in:

  if(test_term %in% vars) {
    drops <- parsed$labels[parsed$factors[test_term, ] > 0]
    return(drops)
  }

  ## test_term names a term of the model: drop it alone, but only if no
  ##   higher-order term contains it:

  if(test_term %in% parsed$labels) {
    vars_test <- vars[parsed$factors[, test_term] > 0]
    for(lbl in setdiff(parsed$labels, test_term)) {
      vars_lbl <- vars[parsed$factors[, lbl] > 0]
      if(all(vars_test %in% vars_lbl)) {
        f.err("f.test_term_drops: test_term '", test_term,
          "' is contained in higher-order term '", lbl,
          "' of frm, so cannot be tested on its own; test the variables of",
          "'", test_term, "' individually, or test '", lbl, "' instead",
          config=config)
      }
    }
    return(test_term)
  }

  f.err("f.test_term_drops: test_term '", test_term,
    "' is neither a variable nor a term of frm; variables:",
    paste(vars, collapse=" "), "; terms:",
    paste(parsed$labels, collapse=" "), config=config)
}

## Helper for f.design_test_cols(). Character scalar config$test_term, formula
##   config$frm; returns list with canonicalized character scalar $test_term,
##   tokenized character vector $frm_terms carrying either '1' for intercept or
##   '0' for no intercept, and character vector $drop_terms with the terms to be
##   dropped from the full model to form the reduced model. Interaction terms in
##   $frm_terms and $test_term are sorted alphabetically (so 'sex:age' becomes
##   'age:sex'), to facilitate formula/term comparison.

f.normalize_terms <- function(config) {

  test_term <- config$test_term
  frm <- config$frm
  
  if(is.null(test_term) || test_term %in% "") {
    f.err("f.normalize_terms: config$test_term empty or undefined.", config=config)
  }
  
  if(is.null(frm) || all(as.character(frm) %in% "")) {
    f.err("f.normalize_terms: config$frm empty or undefined.", config=config)
  }
  
  if(length(test_term) != 1) {
    f.err("f.normalize_terms: length(test_term) != 1; test_term: '", 
      paste(test_term), "'", config=config)
  }
  
  if(grepl("[\\*\\-\\|\\(\\)\\^/]", test_term)) {
    f.err("f.normalize_terms: cannot handle '*', '-', '|', '^', '/', '(', or ')'",
      " in test_term: '", test_term, "'; test_term must name a single variable",
      " or a single interaction term, e.g. 'age' or 'age:sex'", config=config)
  }

  test_term <- gsub("[[:space:]]", "", test_term)

  if(test_term %in% "0") {
    f.err("f.normalize_terms: invalid test term: '", test_term, "'",
      config=config)
  }

  test_term <- f.canon_label(test_term)

  parsed <- f.parse_frm(frm, config)
  frm_terms <- f.formula2terms(config)    ## returns character vector

  ## terms of the full model to be dropped to form the reduced model; throws
  ##   an error if test_term is not compatible with frm:

  drop_terms <- f.test_term_drops(parsed, test_term, config)

  return(list(test_term=test_term, frm_terms=frm_terms, drop_terms=drop_terms))
}

## Apply config$contrast to a limma fit, for the three engines built on one:
##   limma::contrasts.fit() re-expresses the fit in terms of the contrasts it is
##   given, leaving one coefficient per contrast, so with the single contrast
##   config$contrast names the fit that comes back has one column and the moderated
##   test on it is the test of the contrast. The statistic and the logFC are the
##   engine's own throughout, unlike the joint msqrob case; see f.msqrob_wald().
##   Must run between limma::lmFit() and limma::eBayes(), since the moderation is of
##   the variance of the fit as re-expressed. Returns the fit and the coefficient(s)
##   to report: that one column for a contrast, and the columns carrying the test of
##   config$test_term otherwise, when the fit is returned untouched:

f.limma_contrast_fit <- function(fit, design, config) {

  if(is.null(design$contrast)) return(list(fit=fit, coef=design$cols_test))

  con <- matrix(design$contrast, ncol=1,
    dimnames=list(colnames(design$X), trimws(config$contrast)))

  return(list(fit=limma::contrasts.fit(fit, contrasts=con), coef=1L))
}

## Helper for test_lm(). Numeric vector y with one feature's expression across
##   all observations, design matrix X for config$frm over those same
##   observations, integer vector cols_test indexing the columns of X that carry
##   the test of config$test_term, integer vector cols_report indexing the
##   columns whose coefficients are to be returned; returns a named numeric
##   vector with the $pval of a likelihood ratio test followed by the
##   coefficient estimates of the full model for cols_report.
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
##   full model in exactly the same way, and the same likelihood ratio test of the
##   two is the test of the contrast:

f.test_lm_feat <- function(y, X, X_red, cols_report) {

  i <- !is.na(y)
  yy <- y[i]
  xf <- X[i, , drop=F]
  xr <- X_red[i, , drop=F]

  fit_full <- stats::lm(yy ~ xf + 0)

  ## a reduced model with no columns left is the null model with no parameters,
  ##   which lm() will not express as a matrix with zero columns:

  if(ncol(xr) %in% 0) {
    fit_reduced <- stats::lm(yy ~ 0)
  } else {
    fit_reduced <- stats::lm(yy ~ xr + 0)
  }

  tbl <- lmtest::lrtest(fit_full, fit_reduced)
  pval <- tbl[["Pr(>Chisq)"]][2]

  ## lm() prefixes coefficient names with the name of the matrix it was given:

  coefs <- stats::coef(fit_full)
  names(coefs) <- sub("^xf", "", names(coefs))

  noms <- colnames(X)[cols_report]
  out <- coefs[noms]                ## missing (aliased) coefficients give NA
  names(out) <- noms

  return(c(pval=pval, out))
}

#' Hypothesis testing using the \code{stats::lm()} function
#' @description
#'   Tests for differential expression by fitting full and reduced linear 
#'     models, then statistically compare them using a likelihood ratio test.
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
#'   The two models are compared using a likelihood ratio test, implemented with
#'     the \code{lmtest::lrtest()} function. Raw p-values are adjusted for
#'     multiple testing using \code{stats::p.adjust()}. A feature measured in too
#'     few observations to support the full model yields an \code{NA} p-value,
#'     which is reported as 1; use
#'     \code{h0testr::filter_features_by_estimability()} beforehand to drop such
#'     features instead, which screens the ranks of these same two matrices.
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
#' @return
#'   A list with a \code{$hits} data.frame of results and a \code{$fit} that is
#'     always \code{NULL}, since a separate model is fitted to every feature.
#'     Columns of \code{$hits} are \code{c("feature", "p.adj", "pval")}, then the
#'     full-model coefficient estimates for the columns under test, preceded by
#'     \code{Intercept} when \code{config$frm} has one, then the feature metadata
#'     columns from \code{state$features}. Coefficient columns that are
#'     \code{NA} for every feature are dropped, and \code{:} in an interaction
#'     coefficient name becomes \code{.}.
#' @examples
#' set.seed(101)
#' exprs <- h0testr::sim2(n_samps1=6, n_samps2=6, n_genes=25, 
#'   n_genes_signif=5, fold_change=2)$mat
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
#' config$reference_levels <- c(condition="placebo")
#'
#' ## set up and check covariates and parameters:
#' out <- h0testr::initialize(state, config, minimal=TRUE)
#' out$state <- h0testr::filter_features_by_formula(out$state, out$config)
#' 
#' tbl <- h0testr::test_lm(out$state, out$config)
#' print(tbl)

test_lm <- function(state, config, fdr.method="BY") {

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

  pvals <- hits[, 1, drop=T]
  coefs <- hits[, -1, drop=F]
  n <- apply(coefs, 2, function(v) sum(!is.na(v)))
  coefs <- coefs[, n > 0, drop=F]
  
  i <- is.na(pvals)
  if(any(i)) pvals[i] <- 1.0    ## or maybe runif(sum(i)) or 0.5?
  
  fdrs <- stats::p.adjust(pvals, method=fdr.method)
  
  hits <- data.frame(feature=rownames(hits), p.adj=fdrs, pval=pvals, coefs)
  rownames(hits) <- NULL
  i <- names(hits) %in% "X.Intercept."
  if(any(i)) names(hits)[i] <- "Intercept"
  
  tmp <- state$features
  rownames(tmp) <- tmp[[config$feat_col]]
  hits <- cbind(hits, tmp[hits$feature, , drop=F])
  rownames(hits) <- NULL
  
  return(list(hits=hits, fit=NULL))
}

## The number of features behind each gene, which is the covariate DEqMS moderates the
##   per-gene variance against. Taken from config$n_feats_col when the table carries it,
##   so that a state that has already been through combine_features() still reports the
##   counts of the features that went into it rather than 1 per gene; counted from the
##   gene ids otherwise. Keyed by the ids f.gene_ids() assigns, not by the raw column,
##   so that the counts of the unknown_* genes are found rather than coming back NA:

f.gene_counts <- function(state, config, caller="f.gene_counts") {

  feats <- state$features
  nom <- f.n_feats_col(config)

  if(nom %in% names(feats)) {
    out <- as.integer(feats[[nom]])
    names(out) <- f.gene_ids(feats, config, caller)
    return(out[!duplicated(names(out))])
  }

  n <- table(f.gene_ids(feats, config, caller))
  out <- as.integer(n)
  names(out) <- names(n)

  return(out)
}

#' Hypothesis testing using the \code{DEqMS} package
#' @description
#'   Tests for differential expression using the 
#'     \code{DEqMS::spectraCounteBayes()} function.
#' @details
#'   The \code{DEqMS::spectraCounteBayes()} model is fit to \code{config$frm}
#'     and a moderated t-test is performed for whether the effect of
#'     \code{config$test_term} on \code{state$expression} is zero.
#'   \code{DEqMS} moderates the t-statistic of a single coefficient and has no
#'     F-analogue, so \code{config$test_term} must resolve to a single design matrix
#'     column; it is an
#'     error if it does not. That rules out testing a factor with more than two
#'     levels, and also testing a variable that appears in an interaction, since by
#'     marginality the test then covers every term containing the variable: with
#'     \code{config$frm = ~sex * batch} and \code{config$test_term = "sex"}, the test
#'     is a joint test of \code{sexM} and \code{sexM:batchb2}, which this engine
#'     cannot perform. Use \code{h0testr::test_lm()},
#'     \code{h0testr::test_trend()} or \code{h0testr::test_voom()} for those tests.
#'   \code{config$contrast} is the other way through that restriction, and the only
#'     one that keeps this engine: a contrast is one degree of freedom however many
#'     coefficients it weights, and \code{limma::contrasts.fit()} leaves a fit with
#'     the single coefficient \code{DEqMS} moderates. So a comparison between two
#'     particular levels of a factor with more than two, or between two levels of a
#'     factor that also appears in an interaction, runs here, while the
#'     corresponding \code{config$test_term} does not. What is given up is the joint
#'     hypothesis: a contrast compares the levels it names, holding the other
#'     variables of any higher-order term at their reference level, which
#'     \code{h0testr::new_config()} describes and which is warned about when it
#'     applies.
#'   Stops if every gene/protein-group has the same number of features. The whole
#'     of what \code{DEqMS} adds to \code{limma} is a variance prior fitted against
#'     that count, and \code{DEqMS::spectraCounteBayes()} has no spread to fit its
#'     loess against when the count does not vary. Earlier versions fell through to
#'     \code{h0testr::test_trend()} on the un-aggregated state, which returned
#'     feature level rows under the \code{"deqms"} label together with a \code{fit}
#'     that was not a \code{DEqMS} fit. Use \code{config$test_method="trend"}, which
#'     is the same \code{limma} fit without the count based prior, or supply feature
#'     level data. \code{h0testr::tune()} checks the same condition before testing
#'     and records the combination as untested rather than stopping the sweep.
#'   The counts come from \code{config$n_feats_col} when \code{state$features}
#'     carries it, which \code{h0testr::combine_features()} writes, and are counted
#'     from \code{config$gene_id_col} otherwise. So this method runs on an already
#'     aggregated state as well as on feature level input.
#'   Aggregates peptides internally with \code{h0testr::combine_features()},
#'     which fits an additive model and so requires
#'     \code{config$is_log_transformed} to be \code{TRUE}; \code{DEqMS} is built
#'     on \code{limma} and wants the same scale. Running
#'     \code{h0testr::normalize()} first satisfies both.
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
#'     \code{frm}            \cr \tab Formula (formula) to be fit. \cr
#'     \code{test_term}      \cr \tab Term (character) to be tested for non-zero coefficient. \cr
#'     \code{contrast}      \cr \tab Weighted sum (character scalar) of coefficients of \code{config$frm} to test instead of \code{config$test_term}; "" for none. \cr
#'   }
#' @param trend Logical scalar. Whether \code{limma::eBayes()} should use trended dispersion estimate.
#' @return
#'   A list with components:
#'   \tabular{ll}{
#'     \code{hits}  \cr \tab \code{data.frame} of results; columns: 
#'       \code{c("logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B", "gene", "count", "sca.t", "sca.P.Value", "sca.adj.pval")}.
#'       Initial statistics from \code{limma}. Columns beginning with \code{sca.} added by \code{DEqMS}. \cr
#'     \code{fit}   \cr \tab Model returned by \code{DEqMS::spectraCounteBayes}. \cr
#'   } 
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
#'   ## no grp:sex term here: by marginality, testing "grp" in ~grp+sex+grp:sex is a
#'   ##   joint test of grptrt and grptrt:sexM, and this method tests one
#'   ##   coefficient at a time. Test "grp:sex" to test the interaction itself:
#'   frm=~grp+sex,
#'   test_term="grp",
#'   reference_levels=c(grp="ctl", sex="F")
#' )
#' out <- h0testr::initialize(state, config, minimal=TRUE)
#'
#' ## test_deqms() aggregates peptides internally with combine_features(), which
#' ##   fits an additive model, and DEqMS is built on limma; both want log scale.
#' ##   normalize() does this and sets the flag in a full workflow. Applied after
#' ##   initialize(), so that the raw zeros became NA first:
#' out$state$expression <- log2(out$state$expression + 1)
#' out$config$is_log_transformed <- TRUE
#'
#' ## actual test:
#' result <- h0testr::test_deqms(out$state, out$config)
#' head(result$hits)

test_deqms <- function(state, config, trend=FALSE) {

  check_config(config)
  f.check_state(state, config)

  save_state <- config$save_state
  config$save_state <- FALSE
  ## unqualified, like every other internal call in the package, so that test_deqms()
  ##   also works when the sources are loaded without installing. rescale was TRUE
  ##   here, a raw scale division by the per-feature mean applied to the log scale
  ##   data combine_features() requires; medianPolish() absorbs the log scale form of
  ##   that same centering into its own per-feature effect, so dropping it leaves the
  ##   per-sample effects unchanged while keeping each gene's overall level, which is
  ##   the average expression DEqMS's underlying limma trend reads:
  out <- combine_features(state, config, method="medianPolish", rescale=FALSE)
  config$save_state <- save_state
  
  ## the design and the single column carrying the test. DEqMS::spectraCounteBayes()
  ##   moderates the t-statistic of one coefficient and the package offers no
  ##   F-analogue, so only a test that resolves to one coefficient can be run;
  ##   f.design_test_cols_max() errors otherwise. Selecting by coefficient name
  ##   instead used to hide the shortfall whenever the name match happened to yield
  ##   exactly one column, which is the usual case for a two-level factor or a
  ##   numeric covariate inside an interaction: testing 'sex' in ~sex*batch matched
  ##   sexM alone and quietly dropped sexM:batchb2 from the test:

  design <- f.design_test_cols_max(out$state, out$config, "test_deqms", max_cols=1)

  ## DEqMS's whole contribution is a variance prior fitted against the number of
  ##   features behind each gene, so it has nothing to offer when that number is the
  ##   same for every gene: DEqMS::spectraCounteBayes() has no spread to fit the loess
  ##   against. This used to fall through to test_trend() on the un-aggregated state,
  ##   which returned feature level rows under the deqms label together with a fit that
  ##   was not a DEqMS fit; refused here instead, and f.tune2() checks the same
  ##   condition up front so that a sweep records the combination as untested and
  ##   carries on. After the cap above, so that a hypothesis this engine cannot express
  ##   is still reported as such rather than as a shortage of counts:

  counts <- f.gene_counts(state, config, "test_deqms")

  if(length(unique(counts)) < 2) {
    f.err("test_deqms: every gene has the same number of features, so there is no",
      "spread for DEqMS::spectraCounteBayes() to fit its variance prior against;",
      "\n", "  features per gene:", unique(counts),
      "; genes:", length(counts), "\n",
      "  to fix, either use config$test_method 'trend', which is the same limma fit",
      "without the count based prior, or supply feature level data with varying",
      "numbers of features per gene;", "\n",
      "  config$gene_id_col:", config$gene_id_col,
      "; config$feat_col:", config$feat_col, config=config)
  }

  ## config$contrast is the route to a test DEqMS otherwise cannot run: whatever it
  ##   weights, limma::contrasts.fit() leaves one coefficient to moderate, which is
  ##   all the cap above requires. So a contrast within a multi-level factor, or one
  ##   between two levels of a factor that also appears in an interaction, reaches
  ##   this engine while the corresponding config$test_term does not:

  fit <- limma::lmFit(out$state$expression, design$X)
  lc <- f.limma_contrast_fit(fit, design, out$config)
  idx <- lc$coef

  fit <- limma::eBayes(lc$fit, trend=trend)
  fit$count <- counts[rownames(fit$coefficients)]

  ## the aggregated matrix has one row per gene and f.gene_counts() one count per gene,
  ##   so a miss here means the two disagree about the gene ids rather than that a count
  ##   is genuinely unknown, and DEqMS would fit its prior against NA:

  if(any(is.na(fit$count))) {
    f.err("test_deqms: no feature count for", sum(is.na(fit$count)), "of",
      length(fit$count), "genes of the aggregated expression matrix;", "\n",
      "  first few:", utils::head(rownames(fit$coefficients)[is.na(fit$count)], 5),
      config=config)
  }

  fit  <- DEqMS::spectraCounteBayes(fit, fit.method="loess")
  hits <- DEqMS::outputResult(fit, coef_col=idx)

  return(list(hits=hits, fit=fit))
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
#'       5. Calculate test statistics with \code{msqrob2::hypothesisTest()}, when
#'            \code{config$test_term} resolves to a single design matrix column, and
#'            otherwise with a joint Wald test over the columns carrying the test. \cr
#'     }
#'   \code{msqrob2::hypothesisTest()} returns one table per contrast rather than a
#'     joint test over several, so it can answer only a \code{config$test_term} that
#'     resolves to a single design matrix column. That would rule out testing a factor
#'     with more than two levels, and testing a variable that appears in an
#'     interaction, since by marginality the test then covers every term containing the
#'     variable. The limit is in that function rather than in the fit: each fitted
#'     \code{msqrob2} \code{StatModel} carries the coefficients, the unscaled
#'     covariance of the design, and the moderated variance and posterior degrees of
#'     freedom, so for a test spanning several columns the joint statistic is computed
#'     here from those models: \code{W}, the quadratic form of the coefficients under
#'     test in the inverse of their covariance, that covariance being the unscaled
#'     covariance of those columns times the posterior variance. It is reported as
#'     \code{F = W / length(cols)} on \code{length(cols)} and \code{dfPosterior}
#'     degrees of freedom.
#'   \strong{That joint statistic is computed by \code{h0testr}, not returned by
#'     \code{msqrob2}}, which is the one place in this package where the reported
#'     statistic is not the engine's own. The single column case is untouched and is
#'     still exactly what \code{msqrob2::hypothesisTest()} returns; with one column the
#'     statistic above reduces to the square of the moderated t that function reports,
#'     which is why it is referred to an F rather than to a chi-square.
#'   A gene whose model could not be fit, or whose covariance is singular over the
#'     columns under test, is reported as \code{NA} rather than dropped.
#'   The reported contrast is the tested coefficient under the factor level ordering
#'     \code{config$reference_levels} declares and \code{h0testr::initialize()}
#'     resolves, so for a two level factor \code{logFC} is the non-reference level
#'     minus the reference level. \code{msqrob2::msqrob()} builds its own design, so
#'     that ordering is carried through to it explicitly; a coefficient the fit does
#'     not have is an error here, since \code{msqrob2::hypothesisTest()} answers one
#'     with a table of \code{NA}s rather than a complaint.
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
#'     \code{gene_id_col}      \cr \tab Name of column in \code{state$features} with unique gene/protein group ids. \cr
#'     \code{feat_col}         \cr \tab Name of column in \code{state$features} corresponding to \code{rownames(state$expression)}. \cr
#'     \code{obs_col}          \cr \tab Name of column in \code{state$samples} corresponding to \code{colnames(state$expression)}. \cr
#'     \code{frm}              \cr \tab Formula (formula) to be fit. \cr
#'     \code{test_term}        \cr \tab Term (character) to be tested for non-zero coefficient. \cr
#'     \code{contrast}      \cr \tab Weighted sum (character scalar) of coefficients of \code{config$frm} to test instead of \code{config$test_term}; "" for none. \cr
#'     \code{reference_levels} \cr \tab Named character vector with the reference level of each factor variable in \code{config$frm}. \cr
#'     \code{covariate_types}  \cr \tab Optional; classification of variables in \code{config$frm}, as set by \code{initialize()}. \cr
#'     \code{factor_levels}    \cr \tab Optional; resolved levels of each factor variable, as set by \code{initialize()}. \cr
#'   }
#' @param maxit Integer scalar >= 1. How many iterations to use for \code{rlm} fitting.
#' @return
#'   A list with components:
#'   \tabular{ll}{
#'     \code{hits}  \cr \tab \code{data.frame} of results; columns \code{config$gene_id_col},
#'       \code{c("nNonZero", ".n")}, and then either
#'       \code{c("logFC", "se", "df", "t", "pval", "adjPval")} when one design matrix
#'       column carries the test, or
#'       \code{c("f_statistic", "df1", "df2", "pval", "adjPval")} when several do. \cr
#'     \code{fit}   \cr \tab Model returned by \code{msqrob2::hypothesisTest()}, or by
#'       \code{msqrob2::msqrob()} for a joint test, which does not go through
#'       \code{msqrob2::hypothesisTest()}. \cr
#'   }
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
#'   ## no grp:sex term here, so testing "grp" is the single coefficient grptrt, which
#'   ##   msqrob2::hypothesisTest() reports as a moderated t. Adding grp:sex would make
#'   ##   it a joint test of grptrt and grptrt:sexM, reported as an F computed here:
#'   frm=~grp+sex,
#'   test_term="grp",
#'   reference_levels=c(grp="ctl", sex="F")
#' )
#' out <- h0testr::initialize(state, config, minimal=TRUE)
#'
#' ## actual test:
#' result <- h0testr::test_msqrob(out$state, out$config)
#' head(result$hits)

test_msqrob <- function(state, config, maxit=100) {

  check_config(config)
  f.check_state(state, config)

  exprs <- as.data.frame(state$expression)
  ## msqrob2 reads this as the number of observations a feature was measured in,
  ##   so it counts non-missing values, and NA is the only indicator of a missing
  ##   value; see f.zeros_to_na(). Counting v > 0 instead undercounts on a log
  ##   scale, where a value at or below zero is an ordinary measurement, and
  ##   nNonZero feeds msqrob2's weighting, so the undercount moves real results:

  n_non0 <- apply(exprs, 1, function(v) sum(!is.na(v)))
  exprs <- cbind(fnames=rownames(exprs), exprs)
  rownames(exprs) <- NULL
  obj <- QFeatures::readQFeatures(table=exprs, ecol=2:ncol(exprs), 
    fnames="fnames", name="features")  

  ## SummarizedExperiment::rowData(obj[["features"]]) <- S4Vectors::DataFrame(state$features)
  ## f.gene_ids(), so that QFeatures::aggregateFeatures() below groups by the same gene
  ##   ids combine_features() would have assigned: a feature with a blank or missing
  ##   gene id becomes its own gene rather than all of them pooling into one group, and
  ##   the ids that come back match the gene level metadata table test() builds:

  feats <- state$features
  feats[[config$gene_id_col]] <- f.gene_ids(feats, config, "test_msqrob")

  for(nom in names(feats)) {
    if(is.factor(feats[[nom]])) {
      SummarizedExperiment::rowData(obj[["features"]])[[nom]] <- as.character(feats[[nom]])
    } else {
      SummarizedExperiment::rowData(obj[["features"]])[[nom]] <- feats[[nom]]
    }
  }
  SummarizedExperiment::rowData(obj[["features"]])$nNonZero <- n_non0

  ## msqrob2::msqrob() builds its own design, as model.matrix(config$frm,
  ##   colData(obj)), so the factor covariates of config$frm have to reach colData()
  ##   with the level ordering initialize() resolved: a character column there is
  ##   re-leveled by sorting, and the coefficient this function then asks for by name
  ##   no longer exists, which msqrob2::hypothesisTest() reports as a table of NAs
  ##   rather than as an error. Covariates outside config$frm play no part in the fit,
  ##   so they are passed through as before:

  samps <- state$samples
  types <- f.covariate_types(state, config)
  fvars <- names(types)[types %in% "factor"]

  for(nom in names(samps)) {
    if(nom %in% fvars) {
      SummarizedExperiment::colData(obj)[[nom]] <- f.relevel_covariate(samps[[nom]],
        nom, config, "test_msqrob")
    } else if(is.factor(samps[[nom]])) {
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
  
  ## f.parse_frm()$frm, not config$frm: f.design_test_cols_max() below selects the
  ##   tested column from a design built on the parsed formula, whose interaction
  ##   labels have their variables sorted, and stats::model.matrix() names an
  ##   interaction column in the order the term is written, so ~sex*batch yields
  ##   sexM:batchb2 here and batchb2:sexM there. Same fit either way, but the
  ##   contrast below is matched to the fit by name:

  parsed <- f.parse_frm(config$frm, config)

  obj <- msqrob2::msqrob(object=obj, i="genes", formula=parsed$frm, maxitRob=maxit)

  ## the design and the columns carrying the test. msqrob2::hypothesisTest() returns
  ##   one table per contrast rather than a joint test over several, so it can answer
  ##   only a test that resolves to one coefficient; a test spanning more than one goes
  ##   to f.msqrob_wald() below instead, which computes the joint statistic from the
  ##   same fitted models. See test_deqms() for what selecting by coefficient name used
  ##   to hide:

  design <- f.design_test_cols(state, config)
  cols_pick <- colnames(design$X)[design$cols_test]

  ## msqrob2::hypothesisTest() matches the contrast to the fit by parameter name, and
  ##   answers a name the fit does not have with a column of NAs instead of an error,
  ##   which no caller can distinguish from a test that ran and found nothing. So
  ##   check that the design msqrob2::msqrob() built from colData() above, which is
  ##   the model.matrix() call msqrob2::msqrobLm() makes, does carry the column
  ##   selected here:

  parms <- colnames(stats::model.matrix(parsed$frm,
    data=as.data.frame(SummarizedExperiment::colData(obj))))

  if(!all(cols_pick %in% parms)) {
    f.err("test_msqrob: design matrix column", cols_pick, "carrying the test of",
      f.test_label(design, config), "is not among the parameters",
      "msqrob2::msqrob() fit;", "\n", "parameters fit:",
      paste(parms, collapse=", "), config=config)
  }

  dat <- SummarizedExperiment::rowData(obj[["genes"]])
  dat <- as.data.frame(dat[, c(config$gene_id_col, "nNonZero", ".n"), drop=F])

  ## config$contrast is a weighted sum of coefficients, which is one contrast however
  ##   many it weights, so msqrob2::hypothesisTest() answers it directly and the
  ##   statistic reported is the package's own moderated t rather than the joint F
  ##   below. The weight matrix is the one f.contrast_vector() built, handed over as a
  ##   matrix rather than as text for msqrob2::makeContrast() to re-parse, so that
  ##   there is no second reading of config$contrast to disagree with the first:

  if(!is.null(design$contrast)) {

    nom <- make.names(trimws(config$contrast))

    con <- matrix(design$contrast, ncol=1,
      dimnames=list(colnames(design$X), nom))

    obj <- msqrob2::hypothesisTest(object=obj, i="genes", contrast=con,
      modelColumn="msqrobModels")

    tbl <- SummarizedExperiment::rowData(obj[["genes"]])[[nom]]
    tbl <- cbind(dat[rownames(tbl), , drop=F], tbl)
    tbl <- tbl[order(tbl$adjPval, -abs(tbl$logFC)), ]
    rownames(tbl) <- NULL

    return(list(hits=tbl, fit=obj))
  }

  ## one column: the package's own moderated t, unchanged. More than one: the joint
  ##   Wald test over the same fitted models, whose one column case reproduces the
  ##   branch above exactly, which is asserted in the tests rather than assumed:

  if(length(cols_pick) %in% 1) {

    con <- msqrob2::makeContrast(contrasts=paste0(cols_pick, "=0"),
      parameterNames=c(cols_pick))

    obj <- msqrob2::hypothesisTest(object=obj, i="genes", contrast=con,
      modelColumn="msqrobModels")

    tbl <- SummarizedExperiment::rowData(obj[["genes"]])[[cols_pick]]
    tbl <- cbind(dat[rownames(tbl), , drop=F], tbl)
    tbl <- tbl[order(tbl$adjPval, -abs(tbl$logFC)), ]

  } else {

    tbl <- f.msqrob_wald(obj, cols_pick, config)
    tbl <- cbind(dat[rownames(tbl), , drop=F], tbl)
    tbl <- tbl[order(tbl$adjPval, -tbl$f_statistic), ]
  }

  rownames(tbl) <- NULL

  return(list(hits=tbl, fit=obj))
}

## helper for test_msqrob(): the joint test over several coefficients that
##   msqrob2::hypothesisTest() cannot express. That limit is in the API and not in the
##   fit: hypothesisTest() loops over the columns of a contrast and returns one table
##   per column, but each fitted msqrob2 StatModel carries the coefficients, the
##   unscaled covariance of the design, and the moderated variance and posterior
##   degrees of freedom that hypothesisTest() itself uses, so the joint statistic
##     W = b' (sigma_posterior^2 * V_unscaled[cols, cols])^-1 b
##   follows directly from what the package fit. Reported as F = W / length(cols) on
##   length(cols) and dfPosterior degrees of freedom, not as a chi-square on
##   length(cols): the variance is estimated rather than known, and with one column W
##   is the square of the moderated t msqrob2 reports, whose reference distribution is
##   F(1, dfPosterior) exactly. A chi-square would therefore disagree with msqrob2's
##   own answer for the case msqrob2 handles, by more as dfPosterior falls and as more
##   columns are tested, and always in the anticonservative direction.
##   Unlike every other statistic this package reports, this one is computed here
##   rather than by the engine; test_msqrob() says so in its documentation.
##   Returns one row per gene in rowData() order, with the gene ids as rownames. A
##   model that could not be fit, or whose covariance is singular over the tested
##   columns, gives a row of NAs, which is what hypothesisTest() reports for a fit that
##   failed:

f.msqrob_wald <- function(obj, cols, config) {

  dat <- SummarizedExperiment::rowData(obj[["genes"]])

  if(!("msqrobModels" %in% names(dat))) {
    f.err("f.msqrob_wald: msqrob2::msqrob() left no 'msqrobModels' column of",
      "fitted models; columns present:", paste(names(dat), collapse=", "),
      config=config)
  }

  models <- dat$msqrobModels

  out <- data.frame(f_statistic=rep(as.numeric(NA), length(models)),
    df1=as.numeric(NA), df2=as.numeric(NA), pval=as.numeric(NA))
  rownames(out) <- rownames(dat)

  for(idx in seq_along(models)) {

    beta <- try(msqrob2::getCoef(models[[idx]]), silent=T)
    if(inherits(beta, "try-error") || is.null(names(beta))) next
    if(!all(cols %in% names(beta))) next

    b <- beta[cols]
    if(any(is.na(b))) next

    vcv <- try(msqrob2::getVcovUnscaled(models[[idx]]), silent=T)
    sig <- try(msqrob2::getSigmaPosterior(models[[idx]]), silent=T)
    dfp <- try(msqrob2::getDfPosterior(models[[idx]]), silent=T)

    if(inherits(vcv, "try-error") || inherits(sig, "try-error") ||
      inherits(dfp, "try-error")) next

    if(!is.matrix(vcv) || !all(cols %in% rownames(vcv))) next
    if(length(sig) != 1 || is.na(sig) || !is.finite(sig) || sig <= 0) next
    if(length(dfp) != 1 || is.na(dfp) || dfp <= 0) next    ## dfp may be Inf

    ## solve() rather than a pseudo-inverse: a singular covariance over the tested
    ##   columns means the data cannot separate them for this gene, and reporting NA
    ##   for it matches how the engines treat a coefficient they could not estimate:

    vcv <- vcv[cols, cols, drop=F] * sig^2
    stat <- try(drop(t(b) %*% solve(vcv) %*% b), silent=T)

    if(inherits(stat, "try-error") || !is.finite(stat) || stat < 0) next

    out$f_statistic[idx] <- stat / length(cols)
    out$df1[idx] <- length(cols)
    out$df2[idx] <- dfp
    out$pval[idx] <- stats::pf(stat / length(cols), df1=length(cols), df2=dfp,
      lower.tail=F)
  }

  out$adjPval <- stats::p.adjust(out$pval, method="BH")

  return(out)
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
#'   When \code{config$test_term} resolves to a single design matrix column, that
#'     column is tested by \code{proDA::test_diff()} as a single contrast, and
#'     \code{logFC} is the corresponding coefficient: for a two level factor, the non
#'     reference level minus the reference level declared in
#'     \code{config$reference_levels}. When it resolves to several columns, which by
#'     the marginality rule it does for a factor with more than two levels and for a
#'     variable appearing in an interaction, all of them are tested jointly, as a
#'     likelihood ratio test of the full model against the model with those columns
#'     dropped. That test is reported as an F statistic with \code{logFC} \code{NA},
#'     there being no single difference to report, as for the F tests of
#'     \code{h0testr::test_trend()} and \code{h0testr::test_voom()}.
#'   Testing the intercept (\code{config$test_term} \code{"1"}) is supported:
#'     \code{proDA::proDA()} renames that column \code{Intercept}, and the contrast is
#'     looked up under that name.
#'   The one formula this method cannot take that the others can is a
#'     \code{config$test_term} that leaves no parameters in the reduced model, which
#'     happens when \code{config$frm} suppresses the intercept and
#'     \code{config$test_term} names every remaining term (\code{~0+grp} testing
#'     \code{"grp"}). \code{proDA} cannot fit a model with no parameters, so that is
#'     refused with an informative error; use \code{h0testr::test_lm()},
#'     \code{h0testr::test_trend()}, \code{h0testr::test_voom()} or
#'     \code{h0testr::test_prolfqua()} for it, or keep the intercept in
#'     \code{config$frm}, which is the usual comparison among levels in any case.
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
#'     \code{frm}                  \cr \tab Formula (formula) to be fit. \cr
#'     \code{test_term}            \cr \tab Term (character) to be tested for non-zero coefficient. \cr
#'     \code{contrast}      \cr \tab Weighted sum (character scalar) of coefficients of \code{config$frm} to test instead of \code{config$test_term}; "" for none. \cr
#'     \code{reference_levels}     \cr \tab Reference level of each factor variable in \code{frm}; sets the direction of \code{diff}. \cr
#'     \code{normalization_method} \cr \tab If present and \code{is_log_transformed} unset, used to infer it. \cr
#'   }
#' @param is_log_transformed Logical scalar: whether \code{state$expression} has
#'   been log transformed. Defaults to \code{config$is_log_transformed}, which
#'   \code{h0testr::initialize()} and \code{h0testr::normalize()} maintain;
#'   passing both is an error unless they agree.
#' @param prior_df Strictly positive count (\code{location_prior_df}) indicating number of dfs for prior.
#' @param maxit Strictly positive count indicating maximum number of iterations for \code{proDA::proDA()} algorithm.
#' @return
#'   A list with components:
#'   \tabular{ll}{
#'     \code{hits}  \cr \tab \code{data.frame} of results from \code{proDA::test_diff()}, sorted by
#'       p-value; columns \code{c("name", "pval", "adj_pval", "diff", "t_statistic", "se", "df",
#'       "avg_abundance", "n_approx", "n_obs")} for the single contrast test, and
#'       \code{c("name", "pval", "adj_pval", "f_statistic", "df1", "df2", "avg_abundance",
#'       "n_approx", "n_obs")} for the likelihood ratio test. \cr
#'     \code{fit}   \cr \tab Model returned by \code{proDA::proDA()}. \cr
#'   }
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
#'   ## no grp:sex term here, so testing "grp" is the single coefficient grptrt, and a
#'   ##   fold change is reported. Adding one would, by marginality, make it a joint
#'   ##   test of grptrt and grptrt:sexM, run as a likelihood ratio test reporting an F
#'   ##   statistic and no fold change:
#'   frm=~grp+sex,
#'   test_term="grp",
#'   reference_levels=c(grp="ctl", sex="F")
#' )
#' out <- h0testr::initialize(state, config, minimal=TRUE)
#'
#' ## actual test:
#' result <- h0testr::test_proda(out$state, out$config, is_log_transformed=FALSE)
#' head(result$hits)

test_proda <- function(state, config, is_log_transformed=NULL, prior_df=3, maxit=20) {
  
  is_log_transformed <- f.is_log_transformed(is_log_transformed, config,
    "test_proda")
  
  ## the design columns carrying the test, derived before the fit so that a
  ##   config$test_term that does not fit config$frm is an error before the
  ##   expensive part rather than after it. f.parse_frm()$frm is handed to
  ##   proDA::proDA() below, not config$frm: proDA::proDA() builds its own design,
  ##   as model.matrix(design, col_data), and the reduced model given to
  ##   proDA::test_diff() further down is this design$X with columns removed by
  ##   position, so the two designs have to be the same matrix. They are not for a
  ##   two sided config$frm, whose response model.matrix() would have to find in
  ##   col_data, nor for a formula writing an interaction ahead of its variables,
  ##   which changes the order stats::model.matrix() names the interaction column in:

  design <- f.design_test_cols(state, config)
  cols_pick <- colnames(design$X)[design$cols_test]

  fit <- proDA::proDA(state$expression, design=design$parsed$frm, col_data=state$samples,
    data_is_log_transformed=is_log_transformed, location_prior_df=prior_df, max_iter=maxit)

  ## proDA::proDA() renames '(Intercept)' to 'Intercept'; past that its design
  ##   should be column for column the one above, being the same model.matrix()
  ##   call on the same data. Checked rather than assumed, since a divergence would
  ##   have the reduced model below drop columns by position from one design and be
  ##   compared against the other, which is a different hypothesis than the one
  ##   config$test_term names and than filter_features_by_estimability() screened
  ##   features against. The Wald branch would catch it looking its coefficient up
  ##   by name; the likelihood ratio branch has no name to look up:

  cols_want <- colnames(design$X)
  cols_want[cols_want %in% "(Intercept)"] <- "Intercept"
  cols_have <- colnames(proDA::design(fit))

  if(!identical(cols_want, cols_have)) {
    f.err("test_proda: the design proDA::proDA() built does not match the design",
      "h0testr derived from config$frm;", "\n", "proDA::design(fit):", cols_have,
      "\n", "h0testr:", cols_want, config=config)
  }

  ## config$contrast goes to the likelihood ratio branch below whatever it weights.
  ##   proDA::test_diff() does take a contrast expression, so a one degree of freedom
  ##   Wald test is available, but the constrained design f.design_contrast() has
  ##   already built reaches the same hypothesis through the interface used here for
  ##   every joint test, with no second parsing of the contrast by proDA to disagree
  ##   with f.contrast_vector() about what was weighted:

  if(is.null(design$contrast) && length(design$cols_test) %in% 1) {

    ## one design column carries the test, so the Wald test on that coefficient is
    ##   the test config$test_term names. See test_deqms() for what selecting by
    ##   coefficient name used to hide. proDA::result_names() quotes coefficient
    ##   names containing ':', and knows the intercept under the name proDA::proDA()
    ##   renamed it to, which is the column config$test_term '1' selects:

    col_pick <- cols_pick
    if(col_pick %in% "(Intercept)") col_pick <- "Intercept"
    if(grepl(":", col_pick)) col_pick <- paste0("`", col_pick, "`")

    cols2 <- proDA::result_names(fit)
    if(!(col_pick %in% cols2)) {
      f.err("test_proda: coefficient", col_pick, "carrying the test of",
        "config$test_term", config$test_term, "is not among",
        "proDA::result_names(fit):", cols2, config=config)
    }

    tbl <- proDA::test_diff(fit, contrast=col_pick, sort_by="pval")

  } else {

    ## several design columns carry the test, which no single contrast can express,
    ##   so compare the full model against the reduced model f.design_test_cols()
    ##   formed by dropping them: a likelihood ratio test over all of them at once,
    ##   which is the test the marginality rule makes config$test_term mean.
    ##   Handed over as a matrix rather than as formula text, so that the models
    ##   compared are the two f.design_test_cols() already derived, with no second
    ##   derivation of factor level ordering or of interaction column naming to
    ##   disagree with the first. Reported as an F statistic with no fold change,
    ##   like the other multiple coefficient tests, there being no single
    ##   difference to report:

    x_red <- design$X_red

    ## no columns left in the reduced model, which f.design_test_cols() has already
    ##   warned about: config$frm suppresses the intercept and config$test_term names
    ##   every remaining term, so the test is of whether the means are all zero rather
    ##   than of whether they differ. The OLS and limma based methods fit a model with
    ##   no parameters and run that test; proDA::proDA() cannot, and
    ##   proDA::test_diff() reports the empty reduced model from deep inside itself as
    ##   "'d' must be a nonempty numeric vector", which names neither the reduced model
    ##   nor config$frm:

    if(ncol(x_red) %in% 0) {
      f.err("test_proda: dropping config$test_term '", config$test_term,
        "' leaves a reduced model with no parameters, which proDA cannot fit;", "\n",
        "config$frm:", deparse(design$parsed$frm), "; columns carrying the test:",
        paste0(paste(cols_pick, collapse=", "), ";"), "\n",
        "keep the intercept in config$frm for the usual comparison among levels, or",
        "use test_method 'lm', 'trend', 'voom' or 'prolfqua' to test against zero",
        config=config)
    }

    ## proDA::test_diff() requires a full rank reduced model, and reports a rank
    ##   deficient one as colinear covariates, which says nothing about which part
    ##   of config$frm is responsible:

    rank_red <- f.design_rank(x_red)

    if(rank_red < ncol(x_red)) {
      f.err("test_proda: the reduced model for", f.test_label(design, config),
        "is rank deficient, so the likelihood ratio test",
        "against it is not defined;", "\n", "config$frm:",
        deparse(design$parsed$frm), "; reduced model columns:", colnames(x_red),
        "; rank:", rank_red, config=config)
    }

    f.msg("test_proda:", f.test_label(design, config), "carries",
      length(design$cols_test), "design columns (",
      paste(cols_pick, collapse=", "), "), so testing by likelihood ratio against",
      "the reduced model rather than by a single proDA contrast", config=config)

    tbl <- proDA::test_diff(fit, reduced_model=x_red, sort_by="pval")
  }

  return(list(hits=tbl, fit=fit))
}

## The F statistic and p-value for a nested pair of linear models, given the
##   residual sums of squares and residual degrees of freedom of each. Exact for a
##   linear model, and equal to the Wald F-test of the same hypothesis. The error
##   variance and its degrees of freedom are arguments rather than derived from the
##   full fit inside, so that the same code gives the ordinary F when passed the
##   feature's own residual variance and the moderated F when passed the variance and
##   degrees of freedom f.moderate_var() shrinks across features:

f.nested_f <- function(rss_red, rss_full, df_red, df_full, s2_err, df_err) {

  df_num <- df_red - df_full
  ss <- rss_red - rss_full
  ms <- ss / df_num
  fval <- ms / s2_err

  return(data.frame(Df=df_num, Sum.Sq=ss, Mean.Sq=ms, F.value=fval,
    p.value=stats::pf(fval, df_num, df_err, lower.tail=F)))
}

## The error variance and its degrees of freedom, shrunk across features toward a
##   common prior. prolfqua exposes this as prolfqua::squeezeVarRob(), which with
##   robust=FALSE returns exactly what limma::squeezeVar() returns, so moderating
##   here borrows the estimator test_trend(), test_lm() and test_voom() already use
##   rather than introducing a second one. It is also the piece of prolfqua that
##   nothing else in this package can supply: prolfqua::ContrastsModerated() applies
##   it, but only one contrast at a time.
##   Moderation does not depend on how many degrees of freedom are being tested. A
##   moderated F is the ordinary F with the per-feature residual variance replaced by
##   the posterior variance and the denominator degrees of freedom raised by those of
##   the prior, which at one numerator degree of freedom is the square of the
##   moderated t. So a caller needs no branch on df, and the reported statistic comes
##   from one variance estimator whatever config$test_term turns out to carry.
##   Returns its arguments unchanged when there is nothing to borrow: given a single
##   usable feature, prolfqua returns a prior df equal to that feature's own residual
##   df, which would double the denominator degrees of freedom on the strength of one
##   feature's variance.
##   A covariate makes the prior a function of that covariate instead of one number
##   shared by every feature, fitted as a natural spline of up to four degrees of
##   freedom on the log residual variances. Passed mean feature intensity, that is
##   the mean-variance trend of limma::eBayes(trend=TRUE), which is what
##   test_trend() reports:

f.moderate_var <- function(s2, df_resid, config, covariate=NULL,
  who="f.moderate_var") {

  ok <- !is.na(s2) & !is.na(df_resid) & df_resid > 0

  if(sum(ok) < 2) {
    f.msg("WARNING:", who, ": not moderating the error variance:", sum(ok),
      "of", length(s2), "features have a usable residual variance, so there is",
      "nothing to borrow across features", config=config)
    return(list(s2=s2, df=df_resid, df_prior=0, var_prior=as.numeric(NA),
      trend=FALSE))
  }

  ## a spline needs more features than it has degrees of freedom and a covariate
  ##   without spread gives it nothing to bend to, and prolfqua answers either with
  ##   an all-NA prior rather than by falling back, which would lose every p-value.
  ##   Checked here so the run continues on the flat prior and says why. Features
  ##   with a single residual degree of freedom are left out of the count because
  ##   prolfqua leaves them out of the fit:

  if(!is.null(covariate)) {
    n_uniq <- length(unique(covariate))
    spline_df <- min(4, n_uniq)
    n_fit <- sum(ok & df_resid > 1)

    if(length(covariate) != length(s2) || any(!is.finite(covariate)) ||
        n_uniq < 2 || n_fit <= spline_df) {
      f.msg("WARNING:", who, ": fitting the prior variance against a covariate was",
        "requested but the covariate does not support it, so the prior is flat;", "\n",
        "features:", length(s2), "; usable for the fit:", n_fit, "; covariate values:",
        length(covariate), "; distinct:", n_uniq, "; all finite:",
        all(is.finite(covariate)), "; spline degrees of freedom needed:", spline_df,
        config=config)
      covariate <- NULL
    }
  }

  sv <- prolfqua::squeezeVarRob(s2, df=df_resid, covariate=covariate, robust=FALSE)

  if(any(is.infinite(sv$df.prior))) {
    f.msg("WARNING:", who, ": the prior degrees of freedom of the moderated error",
      "variance are infinite, so every feature is assigned its prior error variance",
      "(", f.prior_txt(sv$var.prior), ") and the denominator degrees of freedom are",
      "unbounded;", "\n", "that happens when the per-feature residual variances are",
      "nearly identical, which is expected of simulated data and worth a look",
      "otherwise", config=config)
  }

  return(list(s2=sv$var.post, df=df_resid + sv$df.prior, df_prior=sv$df.prior,
    var_prior=sv$var.prior, trend=!is.null(covariate)))
}

## the prior variance as one line of a log: a number when it is flat, a range when
##   it is a trend, since then there is one value per feature:

f.prior_txt <- function(x) {

  if(length(unique(x)) > 1) {
    return(paste("range", signif(min(x), 4), "to", signif(max(x), 4)))
  }

  return(as.character(signif(x[1], 4)))
}

## Per-feature F-tests comparing prolfqua models of a full and a reduced design,
##   with the error variance moderated across features unless config$test_moderate is
##   FALSE, against a flat prior unless config$test_trend is TRUE and a covariate to
##   fit the prior against is supplied, keyed by config$feat_id_col so that it can be
##   subset to the features that survive the drop below.
##   Features whose data do not support the requested test are dropped explicitly
##   and reported: a rank shortfall means the reduced fit lost fewer columns than
##   the test needs, which happens when missingness leaves a covariate aliased for
##   that feature. Worth being loud about, because both silent alternatives are
##   worse: stats::anova() omits the term for such a feature so it used to vanish
##   from the results unannounced, and the comparison below would otherwise report
##   it with a zero numerator degrees of freedom and an NaN F:

f.prolfqua_nested_f <- function(fit_full, fit_red, design, config, covariate=NULL) {

  idvars <- unique(c(config$gene_id_col, config$feat_id_col))

  a <- as.data.frame(fit_full$modelDF)
  a <- a[, c(idvars, "isSingular", "nrcoef", "df.residual", "sigma"), drop=F]
  b <- as.data.frame(fit_red$modelDF)
  b <- b[, c(idvars, "df.residual", "sigma"), drop=F]
  names(b)[names(b) %in% "df.residual"] <- "df_red"
  names(b)[names(b) %in% "sigma"] <- "sigma_red"

  tbl <- merge(a, b, by=idvars, all=F, sort=F)
  n_in <- nrow(tbl)

  ## sigma is NaN with no residual degrees of freedom left, so that is a shortfall
  ##   too, and an NA in either fit means prolfqua could not fit that feature:

  df_num <- tbl$df_red - tbl$df.residual
  ok <- !is.na(tbl$sigma) & !is.na(tbl$sigma_red) & !is.na(df_num) &
    tbl$df.residual > 0 & df_num %in% design$df_intend

  if(any(!ok)) {
    f.msg("WARNING: f.prolfqua_nested_f: dropping", sum(!ok), "of", n_in,
      "features whose data do not support the test of",
      if(is.null(design$contrast)) {
        paste0("config$test_term '", config$test_term, "';")
      } else {
        paste0("config$contrast '", trimws(config$contrast), "';")
      }, "\n", "the reduced model is not",
      design$df_intend, "degree(s) of freedom below the full model for these",
      "features, so the hypothesis is not estimable for them;", "\n",
      "first few:", paste(utils::head(tbl[[config$feat_id_col]][!ok], 5),
      collapse=", "), config=config)
  }

  if(!any(ok)) {
    f.err("f.prolfqua_nested_f: no feature supports the test of",
      "config$test_term '", config$test_term, "';", "features considered:", n_in,
      "; design columns carrying the test:", length(design$cols_test),
      "; degrees of freedom intended:", design$df_intend, config=config)
  }

  tbl <- tbl[ok, , drop=F]

  ## prolfqua reports stats::sigma(), so the residual sum of squares of each fit is
  ##   sigma^2 times its residual degrees of freedom. The numerator of the F is a
  ##   property of the two designs and is the same either way; only the denominator
  ##   changes under moderation, so both statistics come from one pair of fits:

  rss_red <- tbl$sigma_red^2 * tbl$df_red
  rss_full <- tbl$sigma^2 * tbl$df.residual

  ord <- f.nested_f(rss_red=rss_red, rss_full=rss_full, df_red=tbl$df_red,
    df_full=tbl$df.residual, s2_err=tbl$sigma^2, df_err=tbl$df.residual)

  ## the moderated statistic is the one reported, because borrowing variance across
  ##   features is what prolfqua contributes over a plain per-feature fit and is what
  ##   the limma-family methods here already report. The unmoderated statistic is
  ##   carried alongside it so the two can be compared without refitting, and
  ##   config$test_moderate=FALSE reports the unmoderated one instead:

  moderate <- is.null(config$test_moderate) || isTRUE(config$test_moderate)
  trend <- isTRUE(config$test_trend)

  if(moderate) {

    ## subset after the drop above, so that features the test cannot be run for do
    ##   not contribute to the prior either. A name that is not in the covariate
    ##   gives NA, which f.moderate_var() reports and falls back from:

    cov <- NULL

    if(trend) {
      if(is.null(covariate)) {
        f.err("f.prolfqua_nested_f: config$test_trend is TRUE but no covariate was",
          "supplied to fit the prior variance against", config=config)
      }
      cov <- unname(covariate[as.character(tbl[[config$feat_id_col]])])
    }

    mod <- f.moderate_var(tbl$sigma^2, tbl$df.residual, config, covariate=cov,
      who="f.prolfqua_nested_f")
    stat <- f.nested_f(rss_red=rss_red, rss_full=rss_full, df_red=tbl$df_red,
      df_full=tbl$df.residual, s2_err=mod$s2, df_err=mod$df)
    f.msg("f.prolfqua_nested_f: moderated the error variance across", nrow(tbl),
      "features against a", if(mod$trend) "trended" else "flat", "prior; prior df:",
      signif(mod$df_prior, 4), "; prior variance:", f.prior_txt(mod$var_prior),
      config=config)

  } else {
    mod <- list(s2=tbl$sigma^2, df=tbl$df.residual, df_prior=as.numeric(NA),
      var_prior=as.numeric(NA), trend=FALSE)
    stat <- ord
    f.msg("f.prolfqua_nested_f: config$test_moderate is FALSE, so reporting the",
      "unmoderated F-test", if(trend) paste("and ignoring config$test_trend, which",
      "sets the prior of a moderation that is not being done") else "", config=config)
  }

  ## the term label stays config$test_term, rather than naming the comparison, so
  ##   that result tables and anything keyed on it read the same as before:

  out <- data.frame(tbl[, idvars, drop=F], isSingular=tbl$isSingular,
    nrcoef=tbl$nrcoef, factor=config$test_term, stat, moderated=moderate,
    trend=mod$trend, s2.denom=mod$s2, df.denom=mod$df, df.prior=mod$df_prior,
    s2.prior=rep_len(mod$var_prior, nrow(tbl)),
    F.value.unmod=ord$F.value, p.value.unmod=ord$p.value,
    FDR=stats::p.adjust(stat$p.value, method="BH"), stringsAsFactors=F)

  out <- out[order(out$p.value, decreasing=F), , drop=F]
  rownames(out) <- NULL

  return(out)
}

#' Hypothesis testing using the \code{prolfq} package
#' @description
#'   Tests for differential expression using the \code{prolfq::build_model()} function.
#' @details
#'   Uses the \code{prolfqua::build_model()} function. Returned results sorted by
#'     p-value.
#'   The test is a comparison of nested models: \code{prolfqua::build_model()} fits
#'     the full design and the design with the columns carrying
#'     \code{config$test_term} removed, and the reported F-test is of whether those
#'     columns are all zero. This is exact for a linear model, is the same hypothesis
#'     \code{test_proda()} tests by likelihood ratio, and equals the Wald F-test of
#'     the same contrast. Any \code{config$test_term} is therefore testable, including
#'     one that the marginality rule spreads over several terms: testing \code{"grp"}
#'     in \code{~grp*sex} covers \code{grp} and \code{grp:sex} together, as it does
#'     under \code{test_method} \code{"lm"}, \code{"trend"} and \code{"voom"}.
#'   The design is built by \code{stats::model.matrix()} and its columns are handed to
#'     the fit as numeric columns, so the fit contains no factors: categorical and
#'     continuous covariates go through the same code, and a covariate of either kind
#'     may appear anywhere in \code{config$frm}, including in an interaction. The
#'     consequence for the returned \code{fit} is that it has no \code{xlevels} and
#'     its \code{model} holds indicator columns rather than the original covariates.
#'     Nothing is lost, since \code{colnames(design$X)} are exactly the coefficient
#'     names of the fit: the level ordering set by \code{initialize()} is visible
#'     there as which level of each factor is absent, that being the reference level
#'     declared in \code{config$reference_levels}. The F-test itself does not depend
#'     on the choice of reference level.
#'   The error variance of each feature is shrunk toward a common prior across
#'     features before the F is formed, unless \code{config$test_moderate} is
#'     \code{FALSE}. The shrinkage is \code{prolfqua::squeezeVarRob(robust=FALSE)},
#'     which returns exactly what \code{limma::squeezeVar()} returns, so the reported
#'     statistic is the moderated F of \code{limma}: the ordinary F with the
#'     denominator replaced by the posterior variance and the denominator degrees of
#'     freedom raised by those of the prior. Nothing switches on the numerator degrees
#'     of freedom. At one degree of freedom this is the square of the moderated t that
#'     \code{prolfqua::ContrastsModerated()} reports for a single contrast, and above
#'     one it is what \code{test_trend()} reports for several coefficients, so one
#'     variance estimator covers every \code{config$test_term} and p-values stay
#'     comparable across a sweep of terms of differing degrees of freedom. The
#'     unmoderated F and its p-value are returned alongside, in
#'     \code{F.value.unmod} and \code{p.value.unmod}, so the effect of moderating can
#'     be seen without refitting; with \code{config$test_moderate=FALSE} they are
#'     equal to the reported ones.
#'   The prior that variance is shrunk toward is one number shared by every feature
#'     unless \code{config$test_trend} is \code{TRUE}, which fits it against the mean
#'     intensity of each feature as a natural spline of up to four degrees of freedom
#'     on the log residual variances. That covariate is the one
#'     \code{limma::eBayes(trend=TRUE)} uses, the row mean of the response over the
#'     observations where the feature was seen, so with the trend on this method and
#'     \code{test_trend()} moderate identically and differ only in which package fits
#'     the models. It is off by default, since a mean-variance trend is a property of
#'     the data worth asking for deliberately, and it needs enough features to fit the
#'     spline: given too few, or a mean intensity with no spread, the prior stays flat
#'     and a warning says so, and the \code{trend} column of \code{hits} records what
#'     was actually done. Note the trend is defined on log intensities, so requesting
#'     it for a response that is not log transformed also warns. It is ignored,
#'     with a warning, when \code{config$test_moderate} is \code{FALSE}: it sets the
#'     prior of a shrinkage that is then not performed.
#'   A feature is dropped, with a warning naming how many and which, when its own
#'     data do not support the test: if removing the tested columns does not reduce
#'     the rank of that feature's fit by \code{design$df_intend}, the hypothesis is
#'     not estimable for it. This happens when missingness leaves a covariate
#'     constant, or a factor level unobserved, among the samples where that feature
#'     was measured. \code{filter_features_by_estimability()} screens for the same
#'     condition earlier in the pipeline, so features reach this point only when
#'     \code{test()} is called without it.
#'   The \code{FDR} column of \code{hits} is \code{stats::p.adjust(..., "BH")} over
#'     the features tested, computed here rather than taken from \code{prolfqua},
#'     which adjusts within each row of its per-term ANOVA table.
#'   Flow is:
#'     \tabular{l}{
#'       1. Reshape data into long format with intensities, feature meta, and sample meta. \cr
#'       2. Build the design with \code{f.design_test_cols()} and add its columns to
#'            the long data as numeric columns. \cr
#'       3. Make \code{prolfqua::AnalysisTableAnnotation} and \code{prolfqua::LFQData}
#'            objects. \cr
#'       4. Make two \code{prolfqua::strategy_lm} objects, for the full design and for
#'            the design without the columns carrying \code{config$test_term}. \cr
#'       5. Build a \code{prolfqua} model from each. \cr
#'       6. Moderate the per-feature error variance across features with
#'            \code{prolfqua::squeezeVarRob()}, toward a flat prior or, with
#'            \code{config$test_trend}, one fitted against mean feature intensity. \cr
#'       7. Return the per-feature F-test comparing the two fits, dropping features
#'            for which the test is not estimable. \cr
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
#'     \code{frm}                   \cr \tab Formula (formula) to be fit. \cr
#'     \code{test_term}             \cr \tab Term (character) to be tested for non-zero coefficient. \cr
#'     \code{contrast}      \cr \tab Weighted sum (character scalar) of coefficients of \code{config$frm} to test instead of \code{config$test_term}; "" for none. \cr
#'     \code{reference_levels}      \cr \tab Named character vector with the reference level of each factor variable in \code{config$frm}. \cr
#'     \code{covariate_types}       \cr \tab Optional; classification of variables in \code{config$frm}, as set by \code{initialize()}. \cr
#'     \code{factor_levels}         \cr \tab Optional; resolved levels of each factor variable, as set by \code{initialize()}. \cr
#'     \code{test_moderate}         \cr \tab Optional logical; whether to shrink the error variance across features. Defaults to \code{TRUE} when absent. \cr
#'     \code{test_trend}            \cr \tab Optional logical; whether the prior of that shrinkage is fitted against mean feature intensity rather than flat. Defaults to \code{FALSE} when absent. Unrelated to \code{config$test_method="trend"}. \cr
#'     \code{normalization_method}  \cr \tab If present and \code{is_log_transformed} unset, used to infer it. \cr
#'   }
#' @param is_log_transformed Logical scalar: whether \code{state$expression} has
#'   been log transformed. Defaults to \code{config$is_log_transformed}, which
#'   \code{h0testr::initialize()} and \code{h0testr::normalize()} maintain;
#'   passing both is an error unless they agree.
#' @return
#'   A list with components:
#'   \tabular{ll}{
#'     \code{hits}  \cr \tab \code{data.frame} of the per-feature F-tests of
#'       \code{config$test_term}: the feature id columns
#'       (\code{config$gene_id_col} and \code{config$feat_id_col}, which are the same
#'       column once \code{combine_features()} has run) and:
#'       \code{c("isSingular", "nrcoef", "factor", "Df", "Sum.Sq", "Mean.Sq",
#'       "F.value", "p.value", "moderated", "trend", "s2.denom", "df.denom",
#'       "df.prior", "s2.prior", "F.value.unmod", "p.value.unmod", "FDR")}. One row per
#'       feature, sorted by
#'       \code{p.value}; \code{isSingular} and \code{nrcoef} describe the full fit,
#'       and \code{factor} is \code{config$test_term}. \code{F.value} and
#'       \code{p.value} are the reported test, moderated unless
#'       \code{config$test_moderate} is \code{FALSE}; \code{s2.denom} and
#'       \code{df.denom} are the error variance and denominator degrees of freedom it
#'       used, \code{df.prior} and \code{s2.prior} the degrees of freedom and the value
#'       of the prior the variance was shrunk toward, \code{trend} whether that prior
#'       was fitted against mean feature intensity, in which case \code{s2.prior}
#'       varies by feature, and \code{F.value.unmod} and \code{p.value.unmod} the same
#'       test without moderation. Features for which the test is
#'       not estimable are absent; see Details. \cr
#'     \code{fit}   \cr \tab Model returned by \code{prolfqua::build_model()} for the
#'       full design. \cr
#'     \code{fit_reduced} \cr \tab Model returned by \code{prolfqua::build_model()}
#'       for the design without the columns carrying \code{config$test_term}. \cr
#'     \code{design} \cr \tab The design used, as returned by
#'       \code{f.design_test_cols()}: \code{X} is the full design matrix, whose
#'       column names are the coefficient names of \code{fit}, \code{cols_test}
#'       indexes the columns carrying the test, and \code{df_intend} is the numerator
#'       degrees of freedom the test is meant to have. \cr
#'   }
#'   Note the granularity of \code{hits} follows the input: called on
#'     peptide-level data, as in the example below, it gives one row per peptide,
#'     while in the \code{h0testr::run()} and \code{h0testr::tune()} pipelines
#'     \code{combine_features()} has already run, so it gives one row per
#'     gene/protein group, like every other test method.
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
#'   frm=~grp+sex,
#'   test_term="grp",
#'   reference_levels=c(grp="ctl", sex="F")
#' )
#' out <- h0testr::initialize(state, config, minimal=TRUE)
#'
#' ## actual test:
#' result <- h0testr::test_prolfqua(out$state, out$config, is_log_transformed=FALSE)
#' head(result$hits)
#'
#' ## an interaction is fine too: testing 'grp' in ~grp*sex is a joint test over grp
#' ##   and grp:sex, which the nested model comparison carries as one F-test. The
#' ##   design columns it covers, and the reference level of each factor, being the
#' ##   level with no column of its own:
#' config$frm <- ~grp*sex
#' out <- h0testr::initialize(state, config, minimal=TRUE)
#' result <- h0testr::test_prolfqua(out$state, out$config, is_log_transformed=FALSE)
#' colnames(result$design$X)
#' colnames(result$design$X)[result$design$cols_test]
#' head(result$hits)
#'
#' ## the same 2 df test with the error variance moderated across features, which is
#' ##   the default, next to the same test without it. The moderated p-value is the
#' ##   one reported; both are in every result:
#' head(result$hits[, c("Df", "p.value", "p.value.unmod", "df.denom", "df.prior")])
#' config$test_moderate <- FALSE
#' out <- h0testr::initialize(state, config, minimal=TRUE)
#' plain <- h0testr::test_prolfqua(out$state, out$config, is_log_transformed=FALSE)
#' head(plain$hits[, c("Df", "p.value", "p.value.unmod", "df.denom", "df.prior")])
#'
#' ## and the prior can be fitted against mean feature intensity instead of being one
#' ##   number, which is the prior test_trend() uses. s2.prior is then a value per
#' ##   feature, and the range of it is the fitted mean-variance trend:
#' config$test_moderate <- TRUE
#' config$test_trend <- TRUE
#' out <- h0testr::initialize(state, config, minimal=TRUE)
#' trended <- h0testr::test_prolfqua(out$state, out$config, is_log_transformed=FALSE)
#' head(trended$hits[, c("Df", "p.value", "moderated", "trend", "s2.prior")])
#' range(trended$hits$s2.prior)
#' range(result$hits$s2.prior)                     ## one value, the flat prior

test_prolfqua <- function(state, config, is_log_transformed=NULL) {

  is_log_transformed <- f.is_log_transformed(is_log_transformed, config,
    "test_prolfqua")
  
  parsed <- f.parse_frm(config$frm, config)

  idvars <- unique(c(config$gene_id_col, config$feat_id_col))
  dat <- data.frame(state$features[, idvars, drop=F], state$expression)
  times <- colnames(state$expression)
  i <- grepl("^[[:digit:]]", times)
  times[i] <- paste0("X", times[i])
  if(any(duplicated(times))) {
    f.err("test_prolfqua: cannot get unique observation names in data.frame", "\n",
      "Use observation names that follow R rules for list element naming;",
      "Should start with letter or underbar followed by letter.", config=config
    )
  }
  if(!all(names(dat) %in% c(idvars, times))) {
    f.err("test_prolfqua: !all(names(dat) %in% c(idvars, times))", "\n",
      "names(dat):", names(dat), "\n",
      "c(idvars, times):", c(idvars, times), "\n",
      config=config
    )
  }

  dat <- stats::reshape(dat, direction="long", varying=times, 
    v.names="intensity", idvar=idvars, timevar="sample", times=times)

  rownames(dat) <- NULL
  noms <- colnames(state$expression)
  names(noms) <- times
  dat$sample <- noms[dat$sample]

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
  
  ## the design matrix, built explicitly so that the reduced model can be formed by
  ##   dropping columns from it. Re-deriving the reduced model as formula text
  ##   instead lets stats::model.matrix() re-code the remaining factors to full rank
  ##   and restore the span that was meant to be removed, leaving nothing to test;
  ##   see f.design_test_cols(), which every other engine here uses for the same
  ##   reason. initialize() has already ordered the levels of each factor covariate
  ##   with the declared reference level first, and they are ordered again here so
  ##   that a direct caller who has not run initialize() also gets the declared
  ##   reference level rather than one re-derived by sorting:

  types <- f.covariate_types(state, config)
  st_fit <- state

  for(trm in sort(unique(parsed$vars))) {
    if(types[[trm]] %in% "factor") {
      st_fit$samples[[trm]] <- f.relevel_covariate(state$samples[[trm]], trm,
        config, "test_prolfqua")
    }
  }

  design <- f.design_test_cols(st_fit, config)
  X <- design$X

  ## the covariates reach the fit as the numeric columns of the design matrix, so
  ##   the fit sees no factors at all and categorical and continuous covariates are
  ##   handled by exactly the same code. Nothing is registered in
  ##   prolfqua::AnalysisTableAnnotation$factors, which holds categorical
  ##   annotations only, and nothing needs to be: prolfqua::LFQData$new() takes
  ##   setup=FALSE by default, so dat reaches prolfqua::build_model() verbatim, and
  ##   prolfqua::strategy_lm() fits a plain stats::lm(). The column names are the
  ##   coefficient names, so colnames(design$X) is the record of which level of each
  ##   factor is the reference:

  cols <- make.names(colnames(X), unique=T)
  taken <- intersect(cols, c(idvars, "sample", "intensity"))

  if(length(taken)) {
    f.err("test_prolfqua: design matrix column name(s)", paste(taken, collapse=", "),
      "collide with the columns test_prolfqua() builds;", "\n",
      "rename the covariate(s) in state$samples and config$frm;", "\n",
      "columns built here:", paste(c(idvars, "sample", "intensity"), collapse=", "),
      config=config)
  }

  idx <- match(dat$sample, samps[[config$obs_col]])
  for(k in seq_along(cols)) dat[[cols[k]]] <- X[idx, k]

  ## the reduced model. For config$test_term it is a subset of the columns of X, and
  ##   is named by subsetting, so that nothing about a term test moves. config$contrast
  ##   constrains the model instead of dropping terms from it, so its reduced design is
  ##   a re-parameterization of X rather than a subset of its columns, and those
  ##   columns go into dat under names of their own:

  if(is.null(design$contrast)) {

    cols_red <- cols[-design$cols_test]

  } else {

    cols_red <- paste0("h0red", seq_len(ncol(design$X_red)))
    taken2 <- intersect(cols_red, c(cols, idvars, "sample", "intensity"))

    if(length(taken2)) {
      f.err("test_prolfqua: column name(s)", paste(taken2, collapse=", "),
        "collide with the names test_prolfqua() gives the columns of the reduced",
        "design for config$contrast;", "\n",
        "rename the covariate(s) in state$samples and config$frm", config=config)
    }

    for(k in seq_along(cols_red)) dat[[cols_red[k]]] <- design$X_red[idx, k]
  }

  meta <- prolfqua::AnalysisTableAnnotation$new()
  meta$workIntensity <- "intensity"
  meta$is_response_transformed <- is_log_transformed
  meta$hierarchy[[config$gene_id_col]] <- config$gene_id_col
  meta$hierarchy[[config$feat_id_col]] <- config$feat_id_col

  obj <- prolfqua::LFQData$new(data=dat, config=meta)

  ## the test is the comparison of the full design against the design with the
  ##   columns carrying config$test_term removed, which is an exact F-test for a
  ##   linear model and is the same hypothesis test_proda() performs by likelihood
  ##   ratio. It replaces reading prolfqua's per-term anova table, which could report
  ##   only one term at a time and needed the terms of config$frm reordered so that
  ##   the tested one came last. Nothing moves for a single-term test: the Type I sum
  ##   of squares for a term entered last is by construction the difference in
  ##   residual sums of squares between these two fits. The intercept, if
  ##   config$frm has one, is already a column of X, hence the 0 + :

  ## config$contrast constrains the model rather than dropping terms from it, so its
  ##   reduced design is a re-parameterization of X and not a subset of its columns;
  ##   those columns therefore go into dat under names of their own, alongside the
  ##   columns of X that the full model is fitted to. For config$test_term the
  ##   reduced design is a subset, and is named by subsetting, so that nothing about
  ##   a term test moves:

  frm_full <- paste("intensity ~ 0 +", paste(cols, collapse=" + "))
  frm_red <- if(length(cols_red)) {
    paste("intensity ~ 0 +", paste(cols_red, collapse=" + "))
  } else {
    "intensity ~ 0"
  }

  f.msg("test_prolfqua:", f.test_label(design, config), "; full:", frm_full, "\n",
    "  reduced:", frm_red, "; numerator df:", design$df_intend, config=config)

  fit_full <- prolfqua::build_model(data=obj$data,
    model_strategy=prolfqua::strategy_lm(frm_full),
    subject_Id=obj$config$hierarchy_keys())

  fit_red <- prolfqua::build_model(data=obj$data,
    model_strategy=prolfqua::strategy_lm(frm_red),
    subject_Id=obj$config$hierarchy_keys())

  ## config$test_trend fits the prior variance of the moderation against mean feature
  ##   intensity instead of shrinking every feature toward one number. That covariate
  ##   is limma::eBayes(trend=TRUE)'s: the row mean of the response handed to the fit,
  ##   over the observations where the feature was seen, which limma::lmFit() records
  ##   as Amean. Taken from state$expression rather than from dat so that it is
  ##   visibly that same quantity; dat is a reshape of the same matrix:

  covariate <- NULL

  if(isTRUE(config$test_trend)) {
    covariate <- rowMeans(state$expression, na.rm=T)
    names(covariate) <- as.character(state$features[[config$feat_id_col]])

    if(!is_log_transformed) {
      f.msg("WARNING: test_prolfqua: config$test_trend is TRUE but the response is",
        "not log transformed, so the prior variance is being fitted against mean",
        "untransformed intensity;", "\n", "the mean-variance trend",
        "limma::eBayes(trend=TRUE) models is a trend in log intensity", config=config)
    }
  }

  tbl <- f.prolfqua_nested_f(fit_full, fit_red, design, config, covariate=covariate)

  f.msg("tested", length(unique(tbl[[config$feat_col]])), "features; found",
    sum(tbl$FDR < 0.05, na.rm=T), "hits", config=config)

  return(list(hits=tbl, fit=fit_full, fit_reduced=fit_red, design=design))
}

#' Hypothesis testing using \code{limma::voom}
#' @description
#'   Tests for differential expression using the \code{limma::voom()} function.
#' @details
#'   The \code{limma::voom()} model is fit to \code{config$frm} and a test
#'     is performed for whether the effect of \code{config$test_term} on
#'     \code{state$expression} is zero.
#'   Coefficients are selected exactly as in \code{h0testr::test_trend()}: the design
#'     matrix columns assigned to \code{config$test_term} and to every term
#'     containing it, so naming a variable that also appears in an interaction gives a
#'     joint test over the interaction as well. One coefficient gives a moderated
#'     t-test with a \code{logFC} column, several give an F-test with an \code{F}
#'     column.
#'   Note \code{limma::voom()} models a count mean-variance relationship, so it is
#'     appropriate for count-like input rather than for already log-transformed
#'     abundances.
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
#'     \code{frm}            \cr \tab Formula (formula) to be fit. \cr
#'     \code{test_term}      \cr \tab Term (character) to be tested for non-zero coefficient. \cr
#'     \code{contrast}      \cr \tab Weighted sum (character scalar) of coefficients of \code{config$frm} to test instead of \code{config$test_term}; "" for none. \cr
#'     \code{reference_levels} \cr \tab Named character vector with the reference level of each factor variable in \code{config$frm} (see examples). \cr
#'   }
#' @param normalize.method Character in 
#'   \code{c("TMM", "TMMwsp", "RLE", "upperquartile", "none")}.
#' @return 
#'   A list with components:
#'   \tabular{ll}{
#'     \code{hits}  \cr \tab \code{data.frame} of results from \code{limma::topTable()}. \cr
#'     \code{fit}   \cr \tab Model returned by \code{limma::eBayes()}. \cr
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
#' config <- h0testr::new_config()    ## defaults
#' config$save_state <- FALSE           ## default is TRUE
#' config$feat_col <- config$feat_id_col <- config$gene_id_col <- "feature_id"
#' config$obs_col <- config$obs_id_col <- config$sample_id_col <- "observation_id"
#' config$frm <- ~condition
#' config$test_term <- "condition"
#' config$reference_levels <- c(condition="placebo")
#'
#' ## set up and check configuration, including covariates:
#' out <- h0testr::initialize(state, config, minimal=TRUE)
#' 
#' result <- h0testr::test_voom(out$state, out$config)
#' head(result$hits)

test_voom <- function(state, config, normalize.method="none") {

  check_config(config)
  f.check_state(state, config)
  
  exprs <- state$expression
  i <- apply(exprs, 1, function(v) any(is.na(v)))
  exprs <- exprs[!i, , drop=F]
  
  ## the design and the columns of it carrying the test; see test_trend() for why
  ##   these come from f.design_test_cols() rather than from coefficient names:

  design <- f.design_test_cols(state, config)

  obj <- limma::voom(exprs, design$X, plot=F, normalize.method=normalize.method)
  fit <- limma::lmFit(obj, design$X)
  lc <- f.limma_contrast_fit(fit, design, config)
  fit <- limma::eBayes(lc$fit, trend=F)

  ## a single coefficient gives a t-test and a logFC column; several give an F-test:

  tbl <- limma::topTable(fit, coef=lc$coef, number=Inf)

  f.msg("test_voom:", f.test_label(design, config), "; design columns:",
    ncol(design$X), "; test columns:", length(design$cols_test), "; df:",
    design$df_intend, config=config)
  f.msg("tested", nrow(exprs), "features", config=config)
  f.msg("found", sum(tbl$adj.P.Val < 0.05, na.rm=T), "hits", config=config)
  
  return(list(hits=tbl, fit=fit))
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
#'  The \code{lmFit()} model is fit to \code{config$frm} and a test is
#'    performed on each \code{config$feat_col} for whether the effect of
#'    \code{config$test_term} on \code{state$expression} is zero.
#'   The coefficients carrying that test are the columns of the design matrix
#'     assigned to \code{config$test_term} and to every term containing it, which is
#'     the same selection used by \code{h0testr::test_lm()} and by
#'     \code{h0testr::filter_features_by_estimability()}. So naming a variable that
#'     also appears in an interaction tests the interaction too: with
#'     \code{config$frm = ~sex * batch} and \code{config$test_term = "sex"}, the test
#'     is a joint 2 degree of freedom test of \code{sexM} and \code{sexM:batchb2}.
#'     Testing a factor with more than two levels is likewise a joint test over its
#'     contrasts.
#'   A test of one coefficient is reported by \code{limma::topTable()} as a moderated
#'     t-test with a \code{logFC} column; a test of several is reported as an F-test
#'     with an \code{F} column, and no \code{logFC}, since several coefficients have
#'     no single fold change. Any covariate type \code{stats::model.matrix()}
#'     accepts works, including continuous covariates, interactions among them, and
#'     formulas with no intercept.
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
#'     \code{frm}            \cr \tab Formula (formula) to be fit. \cr
#'     \code{test_term}      \cr \tab Term (scalar character) to be tested for non-zero coefficient. \cr
#'     \code{contrast}      \cr \tab Weighted sum (character scalar) of coefficients of \code{config$frm} to test instead of \code{config$test_term}; "" for none. \cr
#'     \code{reference_levels} \cr \tab Named character vector with the reference level of each factor variable in \code{config$frm} (see examples). \cr
#'   }
#' @return 
#'   A list with components:
#'   \tabular{ll}{
#'     \code{hits}  \cr \tab \code{data.frame} of results from \code{limma::topTable()}. \cr
#'     \code{fit}   \cr \tab Model returned by \code{limma::eBayes()}. \cr
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
#' config <- h0testr::new_config()    ## defaults
#' config$save_state <- FALSE           ## default is TRUE
#' config$feat_col <- config$feat_id_col <- config$gene_id_col <- "feature_id"
#' config$obs_col <- config$obs_id_col <- config$sample_id_col <- "observation_id"
#' config$frm <- ~condition
#' config$test_term <- "condition"
#' config$reference_levels <- c(condition="placebo")
#'
#' ## set up and check covariates and parameters:
#' out <- h0testr::initialize(state, config, minimal=TRUE)
#' 
#' tbl <- h0testr::test_trend(out$state, out$config)
#' print(tbl)

test_trend <- function(state, config) {

  check_config(config)
  f.check_state(state, config)

  if(!is.matrix(state$expression)) {
    f.err("test_trend: !is.matrix(state$expression)", config=config)
  }
  
  ## the design and the columns of it carrying the test, from the same helper
  ##   test_lm() and filter_features_by_estimability() use, so that all three test
  ##   the hypothesis config$test_term names. Selecting coefficients by name instead
  ##   (from the column names of ~0 + test_term) silently violated marginality: those
  ##   names never include a higher-order term containing the test variable, so
  ##   testing 'sex' in ~sex*batch became a 1 df test of sexM alone rather than a
  ##   joint test of sexM and sexM:batchb2:

  design <- f.design_test_cols(state, config)

  fit <- limma::lmFit(state$expression, design$X)
  lc <- f.limma_contrast_fit(fit, design, config)
  fit <- limma::eBayes(lc$fit, trend=T)

  ## a single coefficient gives a t-test and a logFC column; several give an F-test:

  tbl <- limma::topTable(fit, coef=lc$coef, number=Inf)

  f.msg("test_trend:", f.test_label(design, config), "; design columns:",
    ncol(design$X), "; test columns:", length(design$cols_test), "; df:",
    design$df_intend, config=config)
  f.msg("tested", nrow(state$expression), "features", config=config)
  f.msg("found", sum(tbl$adj.P.Val < 0.05, na.rm=T), "hits", config=config)
  
  return(list(hits=tbl, fit=fit))
}

## The effect size reported in the logfc column of the standardized table, from the
##   coefficients of the design matrix columns carrying the test.
##   One column carrying the test means one coefficient, and that coefficient is the
##   signed contrast the engines that report one already report: the difference between
##   the two levels of a factor, or, for a continuous covariate, the change per unit of
##   it. That is returned unchanged, so nothing an engine reports is redefined here.
##   Several columns have no single contrast to report, which is why every engine
##   leaves logfc empty for a joint test. What is reported instead is the total swing:
##   the range, over the observations, of the fitted contribution of the terms under
##   test. It is the largest difference the tested terms can account for between any
##   two observations, and it specializes usefully in each case: for a factor with more
##   than two levels it is the largest difference between any two levels, including the
##   pair that reference level coding does not name a coefficient for; for a continuous
##   covariate it is the slope times the range of the covariate, that is the total
##   change across the observed range rather than a change per unit, which is what
##   makes it comparable with a factor's contrast; for an interaction the columns of
##   the interaction are included, since they are part of what is being tested.
##   It is unsigned: with more than one coefficient there is no single direction to
##   report. In the one column case it would equal the absolute value of the
##   coefficient, so the two are on the same scale.
##   coefs is a features by length(design$cols_test) matrix, its columns in the order
##   of design$cols_test; a feature with any coefficient missing gets NA, since its
##   contribution is then unknown:

f.logfc_effect <- function(coefs, design, config) {

  x_test <- design$X[, design$cols_test, drop=F]

  if(!is.matrix(coefs) || ncol(coefs) != ncol(x_test)) {
    f.err("f.logfc_effect: coefs does not match the columns under test;",
      "columns under test:", ncol(x_test), "; columns of coefs:",
      if(is.matrix(coefs)) ncol(coefs) else paste("not a matrix:", class(coefs)),
      config=config)
  }

  ## config$contrast has a signed effect size of its own, the weighted sum of the
  ##   coefficients it names, which is the quantity being tested and is on the same
  ##   scale as a single coefficient. The engines built on limma report it themselves,
  ##   limma::contrasts.fit() having made it their one coefficient, so this is reached
  ##   for the engines that compare a constrained design instead:

  if(!is.null(design$contrast)) {
    out <- drop(coefs %*% design$contrast[design$cols_test])
    names(out) <- rownames(coefs)
    return(out)
  }

  if(ncol(x_test) %in% 1) {
    out <- coefs[, 1, drop=T]
    names(out) <- rownames(coefs)
    return(out)
  }

  ## contribution of the tested terms to the fitted value of every observation, one
  ##   column per feature:

  contrib <- x_test %*% t(coefs)
  out <- apply(contrib, 2, function(v) diff(range(v)))
  names(out) <- rownames(coefs)

  return(out)
}

## Coefficients of the design matrix columns carrying the test, one row per feature and
##   one column per column of design$cols_test, taken from the fit that produced the
##   p-values rather than from a second fit of our own. NULL when the engine does not
##   expose them, which f.fill_standard() reports rather than silently leaving the
##   effect size empty. Each engine keeps them somewhere different:

f.test_coefs <- function(result, method, design, config) {

  cols <- colnames(design$X)[design$cols_test]

  ## limma moderates the variance and not the coefficients, so fit$coefficients holds
  ##   the least squares estimates for the design that was handed to limma::lmFit(),
  ##   whose columns are the columns of design$X, under the same names:

  if(method %in% c("trend", "voom", "deqms")) {
    coefs <- result$fit$coefficients
    if(!is.matrix(coefs) || !all(cols %in% colnames(coefs))) return(NULL)
    return(coefs[, cols, drop=F])
  }

  ## test_lm() already reports the coefficients under test, one column of its hit
  ##   table each. They arrive there through data.frame(), which applies
  ##   base::make.names() to the design matrix column names, and the intercept is
  ##   renamed from the X.Intercept. that produces:

  if(method %in% "lm") {
    nom <- ifelse(cols %in% "(Intercept)", "Intercept", make.names(cols))
    if(anyDuplicated(nom) || !all(nom %in% names(result$hits))) return(NULL)
    coefs <- as.matrix(result$hits[, nom, drop=F])
    dimnames(coefs) <- list(as.character(result$hits$feature), cols)
    return(coefs)
  }

  ## msqrob2 fits by robust regression and moderates the variance rather than the
  ##   coefficients, so each gene's fitted StatModel carries the estimates for the
  ##   parameters msqrob2::msqrob() built from config$frm, under the same names as the
  ##   columns of design$X; see f.msqrob_wald(), which reads the same models. A model
  ##   that could not be fit leaves NA:

  if(method %in% "msqrob") {

    dat <- SummarizedExperiment::rowData(result$fit[["genes"]])
    if(!all(c("msqrobModels", config$gene_id_col) %in% names(dat))) return(NULL)

    models <- dat$msqrobModels
    coefs <- matrix(as.numeric(NA), nrow=length(models), ncol=length(cols),
      dimnames=list(as.character(dat[[config$gene_id_col]]), cols))

    for(idx in seq_along(models)) {
      beta <- try(msqrob2::getCoef(models[[idx]]), silent=T)
      if(inherits(beta, "try-error") || is.null(names(beta))) next
      coefs[idx, ] <- beta[match(cols, names(beta))]
    }

    return(coefs)
  }

  ## proDA::proDA() renames the intercept column, so the coefficient matrix knows it
  ##   under the name test_proda() looks it up by:

  if(method %in% "proda") {
    coefs <- stats::coefficients(result$fit)
    nom <- ifelse(cols %in% "(Intercept)", "Intercept", cols)
    if(!is.matrix(coefs) || !all(nom %in% colnames(coefs))) return(NULL)
    coefs <- coefs[, nom, drop=F]
    colnames(coefs) <- cols
    return(coefs)
  }

  ## prolfqua::strategy_lm() fits a plain stats::lm() of the response on the columns
  ##   of design$X, under the make.names() forms of their names that test_prolfqua()
  ##   built, one model per feature in fit$modelDF. A model that could not be fit, or
  ##   that dropped a column as non-estimable for that feature, leaves NA:

  if(method %in% "prolfqua") {

    mdf <- as.data.frame(result$fit$modelDF)
    id_col <- intersect(c(config$feat_col, config$gene_id_col, config$feat_id_col),
      names(mdf))
    if(!length(id_col) || !("linear_model" %in% names(mdf))) return(NULL)

    nom <- make.names(colnames(design$X), unique=T)[design$cols_test]
    coefs <- matrix(as.numeric(NA), nrow=nrow(mdf), ncol=length(cols),
      dimnames=list(as.character(mdf[[id_col[1]]]), cols))

    for(idx in seq_len(nrow(mdf))) {
      beta <- try(stats::coef(mdf$linear_model[[idx]]), silent=T)
      if(inherits(beta, "try-error") || is.null(names(beta))) next
      coefs[idx, ] <- beta[match(nom, names(beta))]
    }

    return(coefs)
  }

  return(NULL)
}

## Average feature expression for the expr column of the standardized table: the mean
##   of the values handed to the test, over the observations where the feature was
##   seen, which is the same quantity limma reports as AveExpr. The engines that take
##   peptide level input report gene level results, so their means are over the
##   peptides of each gene as well as over the observations:

f.feature_means <- function(state, method, config) {

  out <- rowMeans(state$expression, na.rm=T)
  out[is.nan(out)] <- NA                  ## a feature seen in no observation

  if(f.gene_level_method(method)) {
    return(tapply(out, f.gene_ids(state$features, config, "f.feature_means"),
      mean, na.rm=T))
  }

  names(out) <- as.character(state$features[[f.test_id_col(method, config)]])

  return(out)
}

## Fill the two columns of the standardized table that an engine's own result table
##   need not carry. Kept out of the f.format_*() functions so that each of those
##   remains a rename of what its engine returned:

f.fill_standard <- function(tbl2, result, state, design, method, config) {

  ## logfc, only where the engine reported none. Whether it did is a property of the
  ##   test rather than of the feature, so this is all rows or none within a run:

  if(any(is.na(tbl2$logfc))) {

    coefs <- f.test_coefs(result, method, design, config)

    if(is.null(coefs)) {
      f.msg("WARNING: test: method", method, "does not expose the coefficients of",
        "the columns carrying the test, so the logfc column is left empty",
        config=config)
    } else {
      eff <- f.logfc_effect(coefs, design, config)
      i <- is.na(tbl2$logfc) & tbl2$feature %in% names(eff)
      tbl2$logfc[i] <- eff[tbl2$feature[i]]
      f.msg("test: logfc:", if(!is.null(design$contrast)) {
        paste("config$contrast", trimws(config$contrast), ", the weighted sum of",
          length(design$cols_test), "coefficients")
      } else if(length(design$cols_test) %in% 1) {
        paste("coefficient of", colnames(design$X)[design$cols_test])
      } else {
        paste("total swing over", length(design$cols_test), "coefficients (",
          paste(colnames(design$X)[design$cols_test], collapse=", "), ")")
      }, config=config)
    }
  }

  if(all(is.na(tbl2$expr))) {
    means <- f.feature_means(state, method, config)
    i <- tbl2$feature %in% names(means)
    tbl2$expr[i] <- means[tbl2$feature[i]]
  }

  return(tbl2)
}

## helper for test():

f.format_lm <- function(tbl, id_col, config) {

  if(!is.data.frame(tbl)) {
    f.err("f.format_lm: !is.data.frame(tbl); class(tbl): ",
      class(tbl), config=config)
  }

  if(!(id_col %in% names(tbl))) {
    f.err("f.format_lm: id_col:", id_col, "not %in% names(tbl); names(tbl):",
      names(tbl), config=config)
  }

  nom <- c("pval", "p.adj")
  if(!all(nom %in% names(tbl))) {
    f.err("f.format_lm: expected names not %in% names(tbl); names(tbl):", 
      names(tbl), "; expected names:", nom, config=config)
  }
  
  tbl <- data.frame(feature=tbl[[id_col]], expr=as.numeric(NA),
    logfc=as.numeric(NA), stat=as.numeric(NA), lod=as.numeric(NA),
    pval=tbl$pval, adj_pval=tbl$p.adj)
  
  tbl <- tbl[order(tbl$pval, decreasing=F), , drop=F]
  rownames(tbl) <- NULL
  
  return(tbl)
}

## helper for test():

f.format_msqrob <- function(tbl, id_col, config) {

  if(!is.data.frame(tbl)) {
    f.err("f.format_msqrob: !is.data.frame(tbl); class(tbl): ",
      class(tbl), config=config)
  }

  if(!(id_col %in% names(tbl))) {
    f.err("f.format_msqrob: id_col:", id_col, "not %in% names(tbl); names(tbl):",
      names(tbl), config=config)
  }

  ## two shapes, since msqrob2::hypothesisTest() returns a moderated t and a log fold
  ##   change for the single contrast test_msqrob() runs when one design column carries
  ##   the test, and f.msqrob_wald() returns an F statistic and no single coefficient
  ##   when several do. Same split, and the same missing logfc, as f.format_proda():

  if(all(c("logFC", "t", "pval", "adjPval") %in% names(tbl))) {

    tbl <- data.frame(feature=tbl[[id_col]], expr=as.numeric(NA),
      logfc=tbl$logFC, stat=tbl$t, lod=as.numeric(NA),
      pval=tbl$pval, adj_pval=tbl$adjPval)

  } else if(all(c("f_statistic", "pval", "adjPval") %in% names(tbl))) {

    tbl <- data.frame(feature=tbl[[id_col]], expr=as.numeric(NA),
      logfc=as.numeric(NA), stat=tbl$f_statistic, lod=as.numeric(NA),
      pval=tbl$pval, adj_pval=tbl$adjPval)

  } else {
    f.err("f.format_msqrob: expected names not %in% names(tbl); names(tbl):",
      names(tbl), config=config)
  }

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

  ## two shapes, since proDA::test_diff() returns a t statistic and a difference
  ##   for the single contrast test_proda() runs when one design column carries the
  ##   test, and an F statistic with no difference for the likelihood ratio test it
  ##   runs when several do. Same split, and the same missing logfc, as
  ##   f.format_limma() has for limma's F test:

  if(all(c("name", "avg_abundance", "diff", "t_statistic", "pval", "adj_pval") %in%
    names(tbl))) {

    tbl <- data.frame(feature=tbl$name, expr=tbl$avg_abundance,
      logfc=tbl$diff, stat=tbl$t_statistic, lod=as.numeric(NA),
      pval=tbl$pval, adj_pval=tbl$adj_pval)

  } else if(all(c("name", "avg_abundance", "f_statistic", "pval", "adj_pval") %in%
    names(tbl))) {

    tbl <- data.frame(feature=tbl$name, expr=tbl$avg_abundance,
      logfc=as.numeric(NA), stat=tbl$f_statistic, lod=as.numeric(NA),
      pval=tbl$pval, adj_pval=tbl$adj_pval)

  } else {
    f.err("f.format_proda: expected names not %in% names(tbl); names(tbl):",
      names(tbl), config=config)
  }

  tbl <- tbl[order(tbl$pval, decreasing=F), , drop=F]
  rownames(tbl) <- NULL
  
  return(tbl)
}

## helper for test():

f.format_prolfqua <- function(tbl, id_col, config) {

  if(!is.data.frame(tbl)) {
    f.err("f.format_prolfqua: !is.data.frame(tbl); class(tbl): ",
      class(tbl), config=config)
  }

  if(!(id_col %in% names(tbl))) {
    f.err("f.format_prolfqua: id_col:", id_col,
      "not %in% names(tbl); names(tbl):", names(tbl), config=config)
  }

  nom <- c("F.value", "p.value", "FDR")
  if(!all(nom %in% names(tbl))) {
    f.err("f.format_prolfqua: expected names not %in% names(tbl); names(tbl):", 
      names(tbl), config=config)
  }
  
  tbl <- data.frame(feature=tbl[[id_col]], expr=as.numeric(NA),
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
    
  if(all(c("AveExpr", "logFC", "t", "B", "P.Value", "adj.P.Val") %in% names(tbl))) {
  
    tbl <- data.frame(feature=rownames(tbl), expr=tbl$AveExpr, 
      logfc=tbl$logFC, stat=tbl$t, lod=tbl$B, 
      pval=tbl$P.Value, adj_pval=tbl$adj.P.Val)
    
  } else if(all(c("AveExpr", "F", "P.Value", "adj.P.Val") %in% names(tbl))) {
  
    tbl <- data.frame(feature=rownames(tbl), expr=tbl$AveExpr, 
      logfc=as.numeric(NA), stat=tbl$F, lod=as.numeric(NA), 
      pval=tbl$P.Value, adj_pval=tbl$adj.P.Val)
    
  } else {
    f.err("f.format_limma: expected names not %in% names(tbl); names(tbl):", 
      names(tbl), config=config)
  }
  
  tbl <- tbl[order(tbl$pval, decreasing=F), , drop=F]
  rownames(tbl) <- NULL
  
  return(tbl)
}

f.format_limma.old <- function(tbl, config) {
  
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
    c("lm", "trend", "deqms", "msqrob", "proda", "prolfqua", "voom")
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
#'     \code{frm}            \cr \tab Formula (formula) to be fit. \cr
#'     \code{test_term}      \cr \tab Term (character scalar) to be tested for non-zero coefficient. \cr
#'     \code{contrast}      \cr \tab Weighted sum (character scalar) of coefficients of \code{config$frm} to test instead of \code{config$test_term}; "" for none. \cr
#'     \code{reference_levels} \cr \tab Named character vector with the reference level of each factor variable in \code{config$frm} (see examples). \cr
#'     \code{test_method}    \cr \tab Character scalar in \code{c("lm", "trend", "deqms", "msqrob", "proda", "prolfqua", "voom")}. \cr
#'   }
#' @param method Name of test method where 
#'   \code{method \%in\% h0testr::test_methods()}.
#' @param is_log_transformed Logical scalar: whether \code{state$expression} has
#'   been log transformed. Only consulted for \code{method \%in\% c("proda",
#'   "prolfqua")}. Defaults to \code{config$is_log_transformed}, which
#'   \code{h0testr::initialize()} and \code{h0testr::normalize()} maintain;
#'   passing both is an error unless they agree.
#' @param prior_df Prior degrees of freedom for method \code{proda}; 
#'   where \code{2 <= prior_df <= n_features}.
#' @return A list with the following elements: \cr
#'   \tabular{ll}{
#'     \code{original} \cr \tab A \code{data.frame} with results in native format returned by test. \cr
#'     \code{standard} \cr \tab A \code{data.frame} with results in a standardized format. \cr
#'     \code{fit}      \cr \tab Fitted model returned by the selected testing procedure. \cr
#'   } \cr
#'   The \code{standard} \code{data.frame} has the following fields: \cr
#'   \tabular{ll}{
#'     \code{feature}   \cr \tab Name of feature tested. \cr
#'     \code{expr}      \cr \tab Average feature expression. \cr
#'     \code{logfc}     \cr \tab Estimated effect size; see below. \cr
#'     \code{stat}      \cr \tab Value of test statistic. \cr
#'     \code{lod}       \cr \tab Log-odds of differential expression. \cr
#'     \code{pval}      \cr \tab Raw p-value resulting from test. \cr
#'     \code{adj_pval}  \cr \tab Adjusted (for multiple testing) p-value. \cr
#'   }
#'   What \code{logfc} holds depends on how many design matrix columns
#'     \code{config$test_term} resolves to, which is a property of
#'     \code{config$frm} and \code{config$test_term} and so is the same for every
#'     row of one result:
#'   \itemize{
#'     \item One column: the coefficient of that column, signed. For a two level
#'       factor that is the difference between its levels on the scale of
#'       \code{state$expression}, so a log fold change when the input is log
#'       transformed. For a \strong{continuous} covariate it is instead the change
#'       \strong{per unit} of that covariate, which is not a fold change between
#'       groups and whose size depends on the units the covariate is recorded in:
#'       an effect per month is a twelfth of the same effect per year. Thresholds
#'       on \code{abs(logfc)} therefore have to be chosen with the covariate's
#'       units in mind, although rankings and p-values are unaffected.
#'     \item More than one column: a joint test of several coefficients has no
#'       single contrast to report, so what is reported is the total swing, the
#'       range over the observations of the fitted contribution of the terms under
#'       test. This is the largest difference those terms can account for between
#'       any two observations: for a factor of more than two levels, the largest
#'       difference between any two of its levels; for a continuous covariate, the
#'       slope times the range of the covariate, that is the total change across
#'       the observed range. It is unsigned, since several coefficients have no one
#'       direction, and it is on the same scale as the single coefficient case.
#'     \item \code{config$contrast}: the weighted sum of the coefficients it names,
#'       signed, which is the quantity being tested and is on the same scale as the
#'       single coefficient case. A contrast is one degree of freedom however many
#'       coefficients it weights, so unlike a joint test it always has one number to
#'       report. See \code{h0testr::new_config()} for how to write one.
#'   }
#'   \code{feature} holds the value of \code{config$feat_col} identifying the row
#'     of \code{state$expression} the result describes, except for \code{method
#'     \%in\% c("deqms", "msqrob")}, which aggregate internally and report one row
#'     per value of \code{config$gene_id_col} whatever level the input is at. So a
#'     method other than those two reports precursors when handed precursors and
#'     genes when handed the output of \code{h0testr::combine_features()}, which
#'     sets \code{config$feat_col} to \code{config$gene_id_col}. The feature
#'     metadata reported in \code{original} is at the matching level: for the two
#'     gene level methods it is the gene level form of \code{state$features},
#'     built exactly as \code{h0testr::combine_features()} builds it.
#'   \code{expr} is the mean of the values handed to the test, over the
#'     observations where the feature was seen; for \code{method \%in\% c("deqms",
#'     "msqrob")}, which take feature level input and report gene level results, it
#'     is also over the features of each gene.
#'   \code{stat} is the statistic the engine itself reports, with one exception: a
#'     joint test with \code{method="msqrob"} reports an F computed by \code{h0testr}
#'     from the fitted \code{msqrob2} models, since
#'     \code{msqrob2::hypothesisTest()} answers one contrast at a time. See
#'     \code{h0testr::test_msqrob()} for what that statistic is. A
#'     \code{config$contrast} run has no such exception: every engine reports its
#'     own statistic for a contrast, including \code{method="msqrob"}, which
#'     answers it with \code{msqrob2::hypothesisTest()}, and
#'     \code{method="deqms"}, which cannot run the corresponding
#'     \code{config$test_term} at all.
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
#'   reference_levels=c(condition="placebo")
#' )
#' 
#' ## set up and check covariates and parameters:
#' out <- h0testr::initialize(state, config, minimal=TRUE)
#' 
#' out <- h0testr::test(out$state, out$config, method="trend")
#' head(out$original)
#' head(out$standard)
#' summary(out$fit)

test <- function(state, config, method=NULL, 
    is_log_transformed=NULL, prior_df=NULL) {
  
  if(is.null(method) || method %in% "") method <- config$test_method
  if(is.null(method) || method %in% "") {
    f.err("test: method and config$test_method both unset", config=config)
  }
  
  if(is.null(prior_df)) prior_df <- config$test_prior_df
  if(method %in% "proda" && is.null(prior_df)) {
    f.err("test: method %in% 'proda' && is.null(prior_df)", config=config)
  }
  
  ## only these two methods are told the scale; resolved here rather than in
  ##   them so that an unusable combination is caught before the fit:

  if(method %in% c("proda", "prolfqua")) {
    is_log_transformed <- f.is_log_transformed(is_log_transformed, config,
      "test")
  }
  
  f.msg("test: method:", method, "; is_log_transformed:", is_log_transformed,
    "; prior_df:", prior_df, config=config)

  ## every method below tests config$test_term against a reduced model, whether
  ##   by dropping terms or by contrasting coefficients, so none of them has
  ##   anything to report when the two models span the same space; checked once
  ##   here rather than in each method:

  design <- f.design_test_cols(state, config)

  ## the column of state$features that identifies the rows the engine will return, and
  ##   so the ids the standardized table is keyed by. Decided once here rather than
  ##   inside each f.format_*(), which used to read a column of its own choosing:
  ##   f.format_lm() read config$gene_id_col although test_lm() is row-wise on
  ##   state$expression, so a state whose features had not been aggregated came back
  ##   labelled with gene ids and failed the check below. See f.test_id_col():

  test_col <- f.test_id_col(method, config)

  if(method %in% "lm") {
    result <- test_lm(state, config)
    tbl2 <- f.format_lm(result$hits, test_col, config)
  } else if(method %in% "trend") {
    result <- test_trend(state, config)
    tbl2 <- f.format_limma(result$hits, config)
  } else if(method %in% "deqms") {
    result <- test_deqms(state, config)
    tbl2 <- f.format_limma(result$hits, config)
  } else if(method %in% "msqrob") {
    result <- test_msqrob(state, config)
    tbl2 <- f.format_msqrob(result$hits, test_col, config)
  } else if(method %in% "proda") {
    result <- test_proda(state, config,
      is_log_transformed=is_log_transformed, prior_df=prior_df)
    tbl2 <- f.format_proda(result$hits, config)
  } else if(method %in% "prolfqua") {
    result <- test_prolfqua(state, config,
      is_log_transformed=is_log_transformed)
    tbl2 <- f.format_prolfqua(result$hits, test_col, config)
  } else if(method %in% "voom") {
    result <- test_voom(state, config)
    tbl2 <- f.format_limma(result$hits, config)
  } else if(method %in% "none") {
    f.msg("skipping testing: method %in% 'none'", config=config)
    return(NULL)
  } else f.err("test: unexpected method:", method, config=config)

  ## the effect size and the average expression, where the engine's own table does not
  ##   carry them; see f.logfc_effect() for what logfc holds for a joint test:

  tbl2 <- f.fill_standard(tbl2, result, state, design, method, config)

  ## the feature metadata to report alongside the result. For a gene level method that
  ##   is the gene level form of the table, built by the same function combine_features()
  ##   uses, so that a gene row carries the same metadata whichever route produced it:

  feats <- state$features
  if(f.gene_level_method(method)) {
    feats <- f.gene_features(feats, config, "test")$features
  }
  rownames(feats) <- feats[[test_col]]

  if(!all(tbl2$feature %in% rownames(feats))) {
    i <- !(tbl2$feature %in% rownames(feats))
    f.err("test: method", method, "returned", sum(i), "of", nrow(tbl2), "results",
      "whose feature id is not in state$features[[", test_col, "]], so they cannot",
      "be matched back to the feature metadata;", "\n",
      "  first few returned:", utils::head(tbl2$feature[i], 5), "\n",
      "  first few available:", utils::head(rownames(feats), 5), "\n",
      "  config$feat_col:", config$feat_col,
      "; config$gene_id_col:", config$gene_id_col, config=config)
  }

  ## where the feature id sits in result$hits, which every engine answers
  ##   differently: the limma family puts it in the rownames, proDA::test_diff() in a
  ##   'name' column, msqrob2 and prolfqua in the id column named by config. The
  ##   f.format_*() functions each know which, but they return only the standardized
  ##   table, so the key is recovered here by taking the first candidate that
  ##   reproduces every id in it. Guessing instead is not visible in the standardized
  ##   table, only in the original one below, which comes back as a full set of NA
  ##   columns; hence taking a candidate only on a complete match, and erroring
  ##   rather than emitting that table when no candidate gives one:

  o <- NULL
  for(key in list(rownames(result$hits), result$hits[["feature"]],
    result$hits[["name"]], result$hits[[test_col]])) {

    if(is.null(key)) next
    o2 <- match(tbl2$feature, as.character(key))
    if(!any(is.na(o2))) {
      o <- o2
      break
    }
  }

  if(is.null(o)) {
    f.err("test: cannot match the standardized results back to the results",
      method, "returned, so the original results cannot be reported;", "\n",
      "feature ids:", utils::head(tbl2$feature, 5), "\n",
      "columns of the returned results:", names(result$hits), config=config)
  }

  tbl <- cbind(feats[tbl2$feature, , drop=F], result$hits[o, , drop=F])
  rownames(tbl) <- NULL
  
  if((!is.null(config$save_state)) && config$save_state) {
    
    file_out <- paste0(config$dir_out, "/", length(config$run_order) + 3, 
      config$result_mid_out, ".reformat", config$suffix_out)
    f.log("writing reformatted results to", file_out, config=config)
    f.save_tsv(tbl2, file_out, config)

    file_out <- paste0(config$dir_out, "/", length(config$run_order) + 3,
      config$result_mid_out, ".original", config$suffix_out)
    f.log("writing original results to", file_out, config=config)
    f.save_tsv(tbl, file_out, config)
  }
  
  return(list(original=tbl, standard=tbl2, fit=result$fit))
}
