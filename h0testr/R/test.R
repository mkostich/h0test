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

## The variance prior DEqMS::spectraCounteBayes() fits, fitted on the genes that can
##   carry one. Its own row handling is silently wrong otherwise: it fits
##   loess(log(fit$sigma^2) ~ log2(fit$count)), and stats::loess() defaults to
##   na.action=na.omit, so a gene whose residual variance is not finite is dropped from
##   the fit. A gene fitted on as many observations as the design has columns has
##   df.residual 0 and sigma NA, which is what unimputed missing values leave behind, so
##   that is reachable from a documented workflow rather than hypothetical.
##   stats::fitted() then comes back shorter than the fit and DEqMS recycles it against
##   fit$df.residual, which hands every gene at or after the first dropped row another
##   gene's prior variance. The recycling restores the length, so
##   f.deqms_moderated_f()'s length check cannot see it, and nothing else says anything:
##   with one ineligible gene among 200, a gene's sca.postvar came back 1.122234 with
##   the ineligible one in row 1 and 1.126142 with it in row 100, the same gene with the
##   same feature count and the same data both times.
##   So the prior is fitted here on the eligible genes alone and the per-gene results
##   are put back where they belong, which leaves the rest without a moderated statistic
##   rather than with a wrong one. That subset fit is DEqMS's own: fitting it this way
##   agrees to all.equal() with DEqMS on an expression matrix that never held the other
##   genes. limma's `[.MArrayLM`() subsets the components it knows about and leaves
##   fit$count alone, which is the same misalignment in miniature, hence setting that by
##   hand below.
##   Genes with no residual degrees of freedom are the ones limma itself cannot test
##   either: their P.Value comes back NA from limma::topTable() for the same reason.
##   Returns fit with the components DEqMS::spectraCounteBayes() adds, each per-gene one
##   filled in for the eligible genes and NA elsewhere, in the row order of
##   fit$coefficients:

f.deqms_prior <- function(fit, config, who="f.deqms_prior") {

  nom <- rownames(fit$coefficients)
  n <- length(nom)

  ## log2(0) is -Inf and would be dropped by na.omit the same way, so the count is
  ##   checked for being usable rather than merely present; test_deqms() has already
  ##   refused an NA count by here:

  ok <- is.finite(fit$sigma) & is.finite(fit$df.residual) & fit$df.residual > 0 &
    is.finite(fit$count) & fit$count > 0

  if(sum(ok) < 2) {
    f.err(who, ": only", sum(ok), "of", n, "genes have a finite residual variance and",
      "a usable feature count, so there is nothing for",
      "DEqMS::spectraCounteBayes() to fit its variance prior from;", "\n",
      "  genes with no residual degrees of freedom:",
      sum(!is.finite(fit$df.residual) | fit$df.residual < 1), "\n",
      "  to fix, impute first (config$impute_method other than 'none'), screen the",
      "features with h0testr::filter_features_by_estimability() and a",
      "config$df_resid_min of at least 1, or use config$test_method 'trend', which is",
      "the same limma fit without the count based prior", config=config)
  }

  ## the up-front guard in test_deqms() counts every gene, and excluding the ineligible
  ##   ones can leave the survivors with a single count between them, which DEqMS itself
  ##   fails on with an error about a missing value in an if() condition:

  if(length(unique(fit$count[ok])) < 2) {
    f.err(who, ": the", sum(ok), "genes with a finite residual variance all have the",
      "same number of features (", unique(fit$count[ok]), "), so there is no spread",
      "for DEqMS::spectraCounteBayes() to fit its variance prior against;", "\n",
      "  the", n - sum(ok), "gene(s) left out of the prior carried the rest of the",
      "spread;", "\n",
      "  to fix, impute first (config$impute_method other than 'none'), or use",
      "config$test_method 'trend', which is the same limma fit without the count",
      "based prior", config=config)
  }

  if(any(!ok)) {
    f.msg("WARNING:", who, ":", sum(!ok), "of", n, "genes are left out of DEqMS's",
      "variance prior and reported without a moderated statistic;", "\n",
      "  no residual degrees of freedom:",
      sum(!is.finite(fit$df.residual) | fit$df.residual < 1),
      "; unusable feature count:", sum(!is.finite(fit$count) | fit$count < 1), "\n",
      "  limma's own columns of those rows are not NA, limma::eBayes() shrinking such a",
      "gene toward a prior of its own that needs no residual variance from the gene;",
      "\n", "  h0testr::impute(), or h0testr::filter_features_by_estimability() with a",
      "config$df_resid_min of at least 1, leaves every gene with residual degrees of",
      "freedom", config=config)
  }

  sub <- fit[ok, ]
  sub$count <- fit$count[ok]

  ## the whole point of this function, so said rather than assumed:

  if(!identical(rownames(sub$coefficients), nom[ok])) {
    f.err(who, ": subsetting the fit did not preserve the gene order, so the prior",
      "cannot be mapped back;", "\n", "  first few expected:",
      utils::head(nom[ok], 5), "; first few found:",
      utils::head(rownames(sub$coefficients), 5), config=config)
  }

  ## the loess needs those counts spread rather than merely varied, and answers two ways
  ##   when almost every gene sits at one of them: its fitted values come back NaN for
  ##   the genes at the crowded count, and when they come back NaN for all of them the
  ##   mean DEqMS matches its prior degrees of freedom against is NaN too and the search
  ##   for them walks off the end of a vector, which surfaces as an error about a missing
  ##   value in an if() condition. Two of sixty genes at one count and the rest at another
  ##   gave the first, twenty genes of two peptides with two single-peptide genes among
  ##   them lost 18 of 20 the same way, and one to three genes of 19, 24 or 60 away from a
  ##   single count gave the second, while an even split fits either way.
  ##   A gene whose posterior variance came back NaN is reported without a moderated
  ##   statistic, like a gene the prior could not be fitted from, and the count of them is
  ##   said in the log; only DEqMS erroring outright is refused, there being no prior at
  ##   all then, and its own message mentions neither the counts nor what to do:

  tb <- table(fit$count[ok])
  tb <- paste(names(tb), as.integer(tb), sep="x", collapse=" ")

  sub <- try(DEqMS::spectraCounteBayes(sub, fit.method="loess"), silent=TRUE)

  if(inherits(sub, "try-error")) {
    f.err(who, ": DEqMS::spectraCounteBayes() could not fit a variance prior from the",
      sum(ok), "genes with residual degrees of freedom;", "\n",
      "  features per gene (count x genes):", tb, "\n",
      "  its loess needs those counts spread rather than merely varied, and fails this",
      "way when nearly every gene sits at one of them;", "\n",
      "  DEqMS said:", trimws(as.character(sub)), "\n",
      "  to fix, use config$test_method 'trend', which is the same limma fit without",
      "the count based prior, or supply feature level data whose numbers of features",
      "per gene are less lopsided", config=config)
  }

  bad <- !is.finite(as.numeric(sub$sca.postvar))

  if(any(bad)) {
    f.msg("WARNING:", who, ":", sum(bad), "of the", sum(ok), "genes the variance prior",
      "was fitted from came back with a non-finite posterior variance, and are reported",
      "without a moderated statistic;", "\n",
      "  features per gene (count x genes):", tb, "\n",
      "  DEqMS's loess needs those counts spread rather than merely varied, and comes",
      "back NaN for the genes at a crowded count when nearly every gene sits at one",
      "of them;", "\n",
      "  config$test_method 'trend' is the same limma fit without the count based prior",
      config=config)
  }

  f.msg(who, ": DEqMS variance prior fitted from", sum(ok), "of", n, "genes over",
    length(unique(fit$count[ok])), "distinct feature counts; prior df:",
    format(sub$sca.dfprior), config=config)

  fit$fit.method <- sub$fit.method
  fit$model <- sub$model
  fit$sca.dfprior <- sub$sca.dfprior          ## one number for the whole fit

  ## a per-gene quantity of the subset fit, back in the row order of the whole fit:

  f.fill <- function(x) {
    out <- rep(as.numeric(NA), n)
    names(out) <- nom
    out[ok] <- as.numeric(x)
    return(out)
  }

  fit$sca.postvar <- f.fill(sub$sca.postvar)
  fit$sca.priorvar <- f.fill(sub$sca.priorvar)

  ## and the same for the gene by coefficient matrices DEqMS forms from them.
  ##   spectraCounteBayes() takes its coef_col missing here, and a missing argument
  ##   passed to `[` is an empty index, so these carry every coefficient rather than one:

  f.fill_mat <- function(x) {
    out <- matrix(as.numeric(NA), nrow=n, ncol=ncol(x),
      dimnames=list(nom, colnames(x)))
    out[ok, ] <- x
    return(out)
  }

  fit$sca.t <- f.fill_mat(sub$sca.t)
  fit$sca.p <- f.fill_mat(sub$sca.p)

  return(fit)
}

## The moderated test of the design matrix columns carrying the test, from a fit
##   DEqMS::spectraCounteBayes() has moderated. DEqMS's own statistic is one
##   coefficient's moderated t, and its API stops there: spectraCounteBayes() takes a
##   coef_col and DEqMS::outputResult() reports one column. The moderation itself does
##   not: that function fits a variance prior against the number of features behind
##   each gene and returns, per gene, a posterior variance in $sca.postvar and a prior
##   degrees of freedom in $sca.dfprior, neither of which mentions a coefficient. Its
##   coef_col enters only in the last two statements, where it forms sca.t and sca.p
##   from them. So the prior is a variance prior of exactly limma's kind, and a joint
##   test over several coefficients is the same substitution limma makes for its own F:
##   the ordinary F with the per-gene residual variance replaced by the posterior
##   variance and the denominator degrees of freedom raised by those of the prior. See
##   f.moderate_var(), which says the same thing about prolfqua's prior.
##   At one numerator degree of freedom that F is the square of DEqMS's own moderated
##   t, so both cases are computed here rather than only the joint one, and the numbers
##   this reports for a single coefficient are DEqMS's own to within floating point.
##   That case is deliberately computed from fit$stdev.unscaled, which limma fits per
##   gene, rather than from fit$cov.coefficients, which it takes from the complete
##   design: the two agree only when every gene was fitted on every observation, and
##   the single coefficient case has to keep working when they do not. Hence also the
##   guard in test_deqms() for the joint case, which has no per-gene covariance matrix
##   available; limma's own F carries the same caveat and does not refuse, but a
##   reported statistic that is quietly wrong for the genes with missing values is
##   worse than one that is not reported.
##   cols indexes the columns of fit$coefficients carrying the test, which after
##   f.limma_contrast_fit() is the single column of the contrast; the fit that comes
##   back from limma::contrasts.fit() carries $cov.coefficients and $stdev.unscaled for
##   that column, so a contrast needs no separate code here.
##   Returns a data.frame with one row per gene, in the row order of fit$coefficients:
##   $t is the signed moderated t when one column carries the test and NA otherwise,
##   $F the moderated F in both cases, and $p.value the p-value of $F. A gene whose
##   posterior variance could not be fitted, or that is missing one of the tested
##   coefficients, comes back NA rather than dropping out of the table:

f.deqms_moderated_f <- function(fit, cols, config, who="f.deqms_moderated_f") {

  need <- c("coefficients", "stdev.unscaled", "df.residual", "sca.postvar",
    "sca.dfprior")

  if(!all(need %in% names(fit))) {
    f.err(who, ": the fit is missing", paste(setdiff(need, names(fit)), collapse=", "),
      ", so it did not come from DEqMS::spectraCounteBayes();", "\n",
      "  names(fit):", names(fit), config=config)
  }

  nom <- colnames(fit$coefficients)[cols]
  betas <- fit$coefficients[, cols, drop=F]
  post_var <- as.numeric(fit$sca.postvar)
  df_den <- as.numeric(fit$sca.dfprior) + as.numeric(fit$df.residual)

  if(length(post_var) != nrow(betas) || length(df_den) != nrow(betas)) {
    f.err(who, ": the moderated variance does not line up with the coefficients;",
      "genes:", nrow(betas), "; posterior variances:", length(post_var),
      "; denominator degrees of freedom:", length(df_den), config=config)
  }

  tval <- rep(as.numeric(NA), nrow(betas))

  if(length(cols) %in% 1) {

    ## DEqMS's own statistic, recomputed: coefficient over stdev.unscaled times the
    ##   square root of the posterior variance, which is what it divides by:

    se <- fit$stdev.unscaled[, cols] * sqrt(post_var)
    tval <- betas[, 1] / se
    fval <- tval^2
    df_num <- 1L

  } else {

    ## the unscaled covariance of the tested coefficients, indexed by name because
    ##   limma drops the columns of a rank deficient design from this matrix rather
    ##   than keeping them as NA, so position and name need not agree:

    V <- fit$cov.coefficients

    if(is.null(V) || !all(nom %in% rownames(V))) {
      f.err(who, ": the fit does not carry an unscaled covariance for the",
        "coefficient(s)", paste(setdiff(nom, rownames(V)), collapse=", "),
        "under test, so the joint test cannot be formed;", "\n",
        "  coefficients under test:", nom, "; available:", rownames(V),
        config=config)
    }

    V <- V[nom, nom, drop=F]
    df_num <- qr(V)$rank

    if(df_num < length(nom)) {
      f.err(who, ": the", length(nom), "coefficient(s) under test (",
        paste(nom, collapse=", "), ") span only", df_num, "dimension(s) of the",
        "design, so the joint test is not estimable;", "\n",
        "  filter_features_by_estimability() screens features against this same",
        "rank, so a whole design that fails it is a property of config$frm",
        config=config)
    }

    Vinv <- try(solve(V), silent=T)
    if(inherits(Vinv, "try-error")) {
      f.err(who, ": the unscaled covariance of the coefficients under test (",
        paste(nom, collapse=", "), ") could not be inverted", config=config)
    }

    ## rowSums() rather than a loop over genes: the same quadratic form
    ##   t(b) %*% Vinv %*% b for every row of betas, and the matrix is small:

    quad <- rowSums((betas %*% Vinv) * betas)
    fval <- (quad / df_num) / post_var
  }

  pval <- stats::pf(fval, df_num, df_den, lower.tail=F)

  ## a gene whose posterior variance came back non-finite, which is a gene f.deqms_prior()
  ##   left out of the prior because it has no residual degrees of freedom for one to be
  ##   fitted from. Said out loud because the p-value is then missing for that gene and
  ##   nothing else reports it. Feature counts with too little spread are a different
  ##   failure, refused up front by test_deqms() and f.deqms_prior() rather than reaching
  ##   here, DEqMS::spectraCounteBayes() erroring rather than returning NA on them:

  lost <- !is.finite(fval)

  if(any(lost)) {
    f.msg("WARNING:", who, ":", sum(lost), "of", length(fval), "genes have no",
      "moderated statistic, so their p-value is NA;", "\n",
      "non-finite posterior variance:", sum(!is.finite(post_var)),
      "; no residual degrees of freedom:",
      sum(!is.finite(fit$df.residual) | fit$df.residual < 1),
      "; missing a coefficient under test:", sum(!stats::complete.cases(betas)),
      config=config)
  }

  out <- data.frame(t=tval, F=fval, df_num=df_num, df_den=df_den, p.value=pval)
  rownames(out) <- rownames(fit$coefficients)

  return(out)
}

#' Hypothesis testing using the \code{DEqMS} package
#' @description
#'   Tests for differential expression using the 
#'     \code{DEqMS::spectraCounteBayes()} function.
#' @details
#'   The \code{DEqMS::spectraCounteBayes()} model is fit to \code{config$frm}
#'     and a moderated test is performed for whether the effect of
#'     \code{config$test_term} on \code{state$expression} is zero.
#'   The coefficients carrying that test are the columns of the design matrix
#'     assigned to \code{config$test_term} and to every term containing it, which is
#'     the same selection used by \code{h0testr::test_lm()},
#'     \code{h0testr::test_trend()} and
#'     \code{h0testr::filter_features_by_estimability()}. So naming a variable that
#'     also appears in an interaction tests the interaction too: with
#'     \code{config$frm = ~sex * batch} and \code{config$test_term = "sex"}, the test
#'     is a joint 2 degree of freedom test of \code{sexM} and \code{sexM:batchb2}.
#'     Testing a factor with more than two levels is likewise a joint test over its
#'     contrasts.
#'   A test of one coefficient is reported as \code{DEqMS}'s own moderated t, with a
#'     \code{logFC} column; a test of several is reported as a moderated F, and no
#'     \code{logFC}, since several coefficients have no single fold change.
#'     \code{DEqMS} itself reports only the former: \code{DEqMS::spectraCounteBayes()}
#'     moderates the t-statistic of one coefficient and the package has no F-analogue
#'     anywhere. Its moderation does not have that limit. That function fits a variance
#'     prior against the number of features behind each gene and returns a per-gene
#'     posterior variance and a prior degrees of freedom, neither of which mentions a
#'     coefficient, so the joint test is the ordinary F with the residual variance
#'     replaced by the posterior one and the denominator degrees of freedom raised by
#'     those of the prior. At one numerator degree of freedom that F is the square of
#'     \code{DEqMS}'s own moderated t, so the single-coefficient case reports
#'     \code{DEqMS}'s statistic unchanged. Earlier versions refused anything but that
#'     case.
#'   A joint test additionally requires that no value of the aggregated expression
#'     matrix be missing, and it is an error if any is. \code{limma} fits each gene on
#'     the observations that gene has, so \code{fit$stdev.unscaled} is per gene while
#'     \code{fit$cov.coefficients}, which a joint test needs, comes from the complete
#'     design; the joint statistic would then be wrong for exactly the genes with
#'     missing values. Running \code{h0testr::impute()} first, which the documented
#'     workflow does, satisfies this. A single coefficient, which
#'     \code{config$contrast} also reduces to, is unaffected.
#'   \code{config$contrast} tests a weighted sum of coefficients instead of a term, is
#'     one degree of freedom however many coefficients it weights, and reaches this
#'     engine as a single coefficient, \code{limma::contrasts.fit()} having made it
#'     one. It is a different hypothesis rather than a way around a restriction: a
#'     contrast compares the levels it names, holding the other variables of any
#'     higher-order term at their reference level, which
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
#'   Counts that vary but are lopsided are the same shortage in weaker form, and
#'     \code{DEqMS::spectraCounteBayes()} answers it two ways. Its loess comes back
#'     \code{NaN} for the genes at a crowded count, which are then reported without a
#'     moderated statistic and counted in the log; twenty genes of two peptides with two
#'     single-peptide genes among them lost 18 of 20 that way. When it comes back
#'     \code{NaN} for every gene the function fails outright, with an error about a missing
#'     value in an \code{if()} condition, and that is refused with the counts named, one to
#'     three genes of 19, 24 or 60 away from a single count being enough. An even split
#'     fits either way. The same fixes apply.
#'   A gene with no residual degrees of freedom, which is a gene the aggregated matrix
#'     has as many observations of as \code{config$frm} has design columns, is left out
#'     of that prior and reported with \code{NA} in every \code{sca.} column, which
#'     is said in the log. \code{limma}'s own columns of that row are not \code{NA}:
#'     \code{limma::eBayes()} shrinks such a gene toward its own prior, which needs no
#'     residual variance of the gene's own, and reports a moderated \code{P.Value} on the
#'     prior's degrees of freedom alone. \code{DEqMS}'s prior is a loess against the
#'     feature counts and has no value to contribute at a gene it was not fitted from.
#'     This is not a restriction so much as
#'     avoidance of a silent error: \code{DEqMS::spectraCounteBayes()} fits its prior
#'     with \code{stats::loess()}, whose default \code{na.action} drops such a gene, and
#'     then recycles the shortened predictions against every gene, so the prior reaching
#'     each gene depends on where the untestable ones sit in
#'     \code{state$expression}. Running \code{h0testr::impute()} first, which the
#'     documented workflow does, or screening features with
#'     \code{h0testr::filter_features_by_estimability()} and a \code{config$df_resid_min}
#'     of at least 1, leaves every gene with residual degrees of freedom.
#'     \code{stats::p.adjust()} takes its \code{n} from the p-values that are not
#'     \code{NA}, so the excluded genes are left out of the multiplicity correction
#'     rather than counted in it.
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
#'       5. Adjust statistics using \code{DEqMS::spectraCounteBayes()}, fitted on the
#'            genes that have residual degrees of freedom. \cr
#'       6. Form the moderated test of the coefficients carrying the test from the
#'            variance prior that fitted. \cr
#'       7. Generate hit table with \code{limma::topTable()}, with the moderated
#'            statistics appended. \cr
#'     }
#'   \code{DEqMS::outputResult()} built the hit table until the joint test was added,
#'     and is no longer used: it takes a single \code{coef_col}, there being nothing
#'     joint for it to report. Every column it produced is still in the table, with
#'     the joint statistic and its degrees of freedom added.
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
#'     \code{test_trend}    \cr \tab Optional logical; answers the \code{trend} argument when that is not given. Defaults to \code{FALSE} when both are absent. Unrelated to \code{config$test_method="trend"}. \cr
#'   }
#' @param trend Logical scalar. Whether \code{limma::eBayes()} fits its variance prior
#'   against mean gene intensity rather than shrinking every gene toward one number.
#'   Defaults to \code{config$test_trend}, and to \code{FALSE} when that is absent, which
#'   is the fit this function performed before the key reached it. \strong{This moves
#'   less than it appears to}: \code{DEqMS::spectraCounteBayes()} fits its own prior from
#'   \code{fit$sigma}, \code{fit$df.residual} and \code{fit$count}, which
#'   \code{limma::eBayes()} does not alter, and the reported statistic comes from that
#'   prior. So \code{trend} sets \code{P.Value}, \code{t}, \code{B}, \code{s2.prior} and
#'   \code{s2.post} of \code{hits}, and leaves every \code{sca.} column, which is what
#'   \code{h0testr::test()} reports, unchanged. For a trended prior that changes the
#'   answer, use \code{config$test_method="trend"}.
#' @return
#'   A list with components:
#'   \tabular{ll}{
#'     \code{hits}  \cr \tab \code{data.frame} of results; columns:
#'       \code{c("logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B", "gene", "count", "sca.t", "sca.F", "sca.df.num", "sca.df.den", "sca.P.Value", "sca.adj.pval")}.
#'       Initial statistics from \code{limma}. Columns beginning with \code{sca.} come from
#'       \code{DEqMS}'s variance prior, and are the ones reported by \code{h0testr::test()}.
#'       \code{sca.t} is \code{NA} when more than one coefficient carries the test, and the
#'       per-coefficient \code{logFC} column is then replaced by one column per tested
#'       coefficient, as \code{limma::topTable()} reports them. They are also \code{NA}
#'       for a gene with no residual degrees of freedom, which has no prior fitted for
#'       it. \cr
#'     \code{fit}   \cr \tab Model returned by \code{DEqMS::spectraCounteBayes}, whose
#'       \code{sca.} components were fitted from the genes with residual degrees of
#'       freedom and are \code{NA} for the rest; \code{fit$model} is the
#'       \code{stats::loess()} of those genes. \cr
#'   }
#'   \code{logFC} is an effect size on the scale of \code{state$expression}: for a two level
#'     factor, the difference between its levels, so a log fold change when the input is log
#'     transformed; for a \strong{continuous} covariate, the change \strong{per unit} of it,
#'     whose size depends on the units the covariate is recorded in. This is what
#'     \code{h0testr::test()} reports as \code{logfc} when one design matrix column carries the
#'     test; for a joint test there is no single contrast, and \code{h0testr::test()} reports the
#'     total swing instead. See \code{h0testr::test()}.
#' @examples
#' ## setup of expression data: ten peptides per gene, a third of them dropped, and no
#' ##   missing values, so that the example is about the test rather than about missingness:
#' set.seed(101)
#' samps <- h0testr::sim_samples(factors=list(grp=c("ctl", "trt"), sex=c("F", "M")),
#'   n_per_cell=3)
#' sim <- h0testr::sim_design(samps, frm=~grp + sex, test_term="grp", n_genes=100,
#'   peps_per_gene=10, p_drop=0.33, mnar_c0=-Inf, mnar_c1=0, mcar_p=0)
#' state <- sim$state
#' config <- sim$config
#' rm(samps, sim)
#'
#' ## no grp:sex term here, so that this example shows the single-coefficient case,
#' ##   which reports DEqMS's own moderated t. By marginality, testing "grp" in
#' ##   ~grp+sex+grp:sex would be a joint test of grptrt and grptrt:sexM, reported as
#' ##   a moderated F; test "grp:sex" to test the interaction itself:
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
#'
#' ## the same fit with limma's own prior fitted against mean gene intensity instead of
#' ##   flat. Equivalently config$test_trend <- TRUE, which h0testr::test(method="deqms")
#' ##   passes through. limma's P.Value moves and DEqMS's sca.P.Value does not, the
#' ##   count-based prior the latter comes from being fitted from quantities
#' ##   limma::eBayes() leaves alone; see the note on the trend argument:
#' trended <- h0testr::test_deqms(out$state, out$config, trend=TRUE)
#' i <- rownames(result$hits)
#' c(limma=cor(result$hits$P.Value, trended$hits[i, "P.Value"]),
#'   deqms=cor(result$hits$sca.P.Value, trended$hits[i, "sca.P.Value"]))

test_deqms <- function(state, config, trend=NULL) {

  ## NULL rather than FALSE, so that config$test_trend reaches this engine: the key
  ##   describes the prior of a moderation and this fit has one to describe, and a
  ##   caller who never touches config$test_trend still gets the untrended fit that
  ##   was the previous default. See f.is_trend():

  trend <- f.is_trend(trend, config, "test_deqms")

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
  
  ## the design and the columns carrying the test, from the same helper test_lm() and
  ##   filter_features_by_estimability() use, so that all three test the hypothesis
  ##   config$test_term names. This was capped at a single column, because
  ##   DEqMS::spectraCounteBayes() reports one coefficient's moderated t and the package
  ##   has no F-analogue; the moderation itself has no such limit, and
  ##   f.deqms_moderated_f() forms the joint test from what that function returns.
  ##   Selecting by coefficient name instead of by design column would still be wrong,
  ##   for the reason it always was: the name match happens to yield exactly one column
  ##   for a two-level factor or a numeric covariate inside an interaction, so testing
  ##   'sex' in ~sex*batch matched sexM alone and quietly dropped sexM:batchb2:

  design <- f.design_test_cols(out$state, out$config)

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

  ## a joint test needs the covariance of the tested coefficients, and limma supplies
  ##   only one such matrix for the whole fit, taken from the complete design: with a
  ##   missing value it fits each gene on the observations that gene has, so
  ##   fit$stdev.unscaled is per gene but fit$cov.coefficients is not, and the joint
  ##   statistic would be wrong for exactly the genes with missing values. limma's own
  ##   F has the same caveat and does not refuse; refused here rather than reported
  ##   wrong, since h0testr::impute() runs before h0testr::test() in the documented
  ##   workflow and so this is normally already satisfied. A single coefficient, which a
  ##   contrast also reduces to, is unaffected: f.deqms_moderated_f() reads
  ##   fit$stdev.unscaled there. Per-missingness-pattern recomputation is the way to
  ##   lift this if it turns out to bite:

  if(is.null(design$contrast) && length(design$cols_test) > 1 &&
      anyNA(out$state$expression)) {

    f.err("test_deqms: testing config$test_term '", config$test_term, "' in",
      deparse(design$parsed$frm), "is a joint test of", length(design$cols_test),
      "coefficients (",
      paste(colnames(design$X)[design$cols_test], collapse=", "), "), and", "\n",
      "  ", sum(is.na(out$state$expression)), "of",
      length(out$state$expression), "values of the aggregated expression matrix",
      "are missing, which leaves limma no per-gene covariance for the coefficients",
      "under test;", "\n",
      "  to fix, impute first (config$impute_method other than 'none'), or use",
      "config$contrast, or config$test_method 'trend', 'voom', 'lm', 'msqrob',",
      "'msqrob_agg', 'prolfqua', 'prolfqua_lmer' or 'proda'", config=config)
  }

  ## config$contrast reaches this engine by a shorter route: whatever it weights,
  ##   limma::contrasts.fit() leaves one coefficient, and the fit it returns carries
  ##   $cov.coefficients and $stdev.unscaled for that coefficient, so nothing below
  ##   needs to know which of the two kinds of hypothesis is being tested:

  fit <- limma::lmFit(out$state$expression, design$X)

  ## limma::eBayes() refuses a fit in which no gene has a residual degree of freedom, and
  ##   f.deqms_prior() cannot fit a variance prior from fewer than two such genes. Checked
  ##   here so that what to do about it is said: limma's own message for this is "No
  ##   residual degrees of freedom in linear model fits", which does not mention that
  ##   imputing or filtering first is what this engine needs:

  if(sum(fit$df.residual > 0, na.rm=T) < 2) {
    f.err("test_deqms: only", sum(fit$df.residual > 0, na.rm=T), "of",
      length(fit$df.residual), "genes of the aggregated expression matrix have a",
      "residual degree of freedom, so there is nothing to fit a variance prior from;",
      "\n", "  design columns:", ncol(design$X), "; observations:",
      ncol(out$state$expression), "\n",
      "  to fix, impute first (config$impute_method other than 'none'), screen the",
      "features with h0testr::filter_features_by_estimability() and a",
      "config$df_resid_min of at least 1, or drop terms from config$frm", config=config)
  }

  lc <- f.limma_contrast_fit(fit, design, out$config)
  idx <- lc$coef

  ## said out loud because the two fits differ only in the prior and the hit table does
  ##   not record which was used. Worth being explicit about how little this reaches:
  ##   DEqMS::spectraCounteBayes() fits its own prior from fit$sigma, fit$df.residual and
  ##   fit$count, none of which limma::eBayes() alters, and f.deqms_moderated_f() forms
  ##   the reported statistic from that prior (sca.postvar, sca.dfprior). So the trended
  ##   prior moves limma's own columns of hits, P.Value, t, B, s2.prior and s2.post, and
  ##   leaves every sca.* column, and therefore everything h0testr::test() reports for
  ##   this method, exactly where it was. Passed through anyway, because the argument is
  ##   limma's to take and the columns it moves are returned to the caller, but a caller
  ##   after a trended prior that changes the answer wants test_method "trend":

  f.msg("test_deqms: limma::eBayes prior fitted against", if(trend) {
    "mean gene intensity (trend=TRUE)"
  } else "one number for every gene (trend=FALSE)", "\n",
    " this sets limma's own columns of the result; the reported sca.* statistics come",
    "from the prior DEqMS fits against the feature counts, which is unaffected",
    config=config)

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

  ## rather than DEqMS::spectraCounteBayes() directly: a gene with no residual degrees
  ##   of freedom is dropped from the loess that fits the prior and the predictions are
  ##   then recycled against every gene, which misaligns the prior silently. See
  ##   f.deqms_prior(), which fits it on the genes that can carry one:

  fit <- f.deqms_prior(fit, out$config, "test_deqms")

  ## the moderated test, and the table it is reported in. DEqMS::outputResult() built
  ##   this table, and is not used: it takes a single coef_col, since there is nothing
  ##   joint for it to report, and it reads fit$sca.t and fit$sca.p, which
  ##   DEqMS::spectraCounteBayes() forms for one coefficient at a time. The columns it
  ##   produced are all still here, with the joint statistic and its degrees of freedom
  ##   added, so that a caller reading the original table sees what it saw before:

  mod <- f.deqms_moderated_f(fit, idx, out$config, "test_deqms")

  hits <- limma::topTable(fit, coef=idx, number=Inf, sort.by="none")
  mod <- mod[rownames(hits), , drop=F]

  hits$gene <- rownames(hits)
  hits$count <- fit$count[rownames(hits)]
  hits$sca.t <- mod$t
  hits$sca.F <- mod$F
  hits$sca.df.num <- mod$df_num
  hits$sca.df.den <- mod$df_den
  hits$sca.P.Value <- mod$p.value

  ## stats::p.adjust() takes n from the number of p-values that are not NA, so a gene
  ##   whose posterior variance could not be fitted is left out of the correction
  ##   rather than counted in it:

  hits$sca.adj.pval <- stats::p.adjust(hits$sca.P.Value, method="BH")

  hits <- hits[order(hits$sca.P.Value, decreasing=F), , drop=F]

  f.msg("test_deqms:", f.test_label(design, config), "; design columns:",
    ncol(design$X), "; test columns:", length(design$cols_test), "; df:",
    design$df_intend, config=config)
  f.msg("tested", nrow(hits), "genes over", nrow(state$expression), "features",
    config=config)
  f.msg("found", sum(hits$sca.adj.pval < 0.05, na.rm=T), "hits", config=config)

  return(list(hits=hits, fit=fit))
}

#' Hypothesis testing using the \code{msqrob2} package
#' @description
#'   Tests for differential expression using the \code{msqrob2::msqrob()} function, or
#'   with \code{aggregate=TRUE} the \code{msqrob2::msqrobAggregate()} function, which
#'   fits one mixed model per gene over the rows of its features instead.
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
#'   With \code{aggregate=TRUE}, which \code{config$test_method="msqrob_agg"} selects,
#'     steps 3 and 4 are replaced by a single call to
#'     \code{msqrob2::msqrobAggregate()}. That function is not an alternative to
#'     aggregating so much as a different use of it: it aggregates to make the assay
#'     that carries the result, one row per gene, but fits \code{msqrob2::msqrobLmer()}
#'     on the un-aggregated assay with the features of a gene as one group, so the
#'     models are feature level and the results gene level. Nothing reads the aggregated
#'     values, which is why \code{base::colMeans()} is passed as the summary rather than
#'     \code{msqrob2}'s default \code{MsCoreUtils::robustSummary()}: the latter reports a
#'     sample whose every value in a group is exactly \code{0} as \code{NA}, silently,
#'     as \code{h0testr::combine_features()} documents, and there is nothing to gain by
#'     inviting that into a quantity no test uses.
#'   The model \code{aggregate=TRUE} fits is
#'     \code{intensity ~ <design columns> + (1|config$feat_id_col) +
#'     (1|config$obs_col)}, the second random term dropped when
#'     \code{config$test_random_obs} is \code{FALSE}. The random feature effect is what
#'     \code{msqrob2} documents: it absorbs the baseline of each feature as one variance
#'     component rather than as one coefficient per feature, which is what lets a
#'     thinly observed feature contribute without dominating. The random observation
#'     effect is added because the covariates come from \code{state$samples} and so vary
#'     across observations rather than across a gene's features: without a term for the
#'     observation, the features of a gene stand as independent measurements of it and
#'     the test is anti-conservative. Simulated under the null with 300 genes, 5 features
#'     per gene and 12 observations in two groups of six, with a per-observation effect
#'     the size of the residual, the rejection rate at the 5\% level was 0.360 with the
#'     feature effect alone and 0.093 with both; with no per-observation effect present,
#'     0.063 and 0.053. \code{config$test_random_obs=FALSE} is kept so that the
#'     structure \code{msqrob2} documents can be compared, and warns when set.
#'   \strong{The residual 0.093 above is \code{msqrob2}'s own calibration, and it is
#'     worth knowing about}: \code{msqrob2} refers its moderated t to
#'     \code{dfPosterior}, which is \code{limma}'s prior degrees of freedom plus an
#'     effective residual df for the whole feature level fit, and that is much larger
#'     than the degrees of freedom the tested contrast actually has. In the simulation
#'     above the median was 253 where the observations supply 9. The covariance of the
#'     fixed effects does account for both random terms, which is what keeps the
#'     rejection rate near the level rather than at 0.360, but the reference
#'     distribution is too generous and the test is mildly anti-conservative. This is
#'     what \code{msqrob2} reports and it is reported unchanged, as
#'     \code{config$test_method="prolfqua"} reports \code{prolfqua}'s own statistic;
#'     \code{config$test_method="prolfqua_lmer"} fits the same random structure and
#'     takes Satterthwaite degrees of freedom for the contrast instead, giving 0.053 on
#'     the same simulation, and is the better calibrated of the two feature level mixed
#'     model paths.
#'   \code{config$test_ridge=TRUE} penalizes the fixed effects,
#'     \code{msqrob2::msqrobAggregate(ridge=TRUE)}. It is \code{FALSE} by default, which
#'     is \code{msqrob2}'s own default in every one of its entry points, for three
#'     reasons: the penalty biases the coefficients toward zero by construction, so
#'     \code{logFC} stops being comparable to what the other methods report; it renames
#'     every fitted fixed effect \code{ridge<column>}, which this function translates
#'     but which any reader of \code{result$fit} will meet; and it refuses a mean model
#'     with fewer than two non-intercept columns, so \code{config$frm=~grp} for a two
#'     level factor stops working and has to be written \code{~0+grp}.
#'   A gene with fewer than two observed features is fitted at the gene level instead,
#'     by \code{msqrob2::msqrob()} on the aggregated values. Its random feature effect
#'     would be a single unknown confounded with the intercept, so that variance is not
#'     identified and \code{lme4} refuses the fit outright, which \code{msqrob2} records
#'     as a \code{fitError} and would report as a row of \code{NA}s. Such genes are
#'     common, and dropping them would leave a gene table that does not line up with
#'     what another \code{config$test_method} produces, so they get the same fixed model
#'     without the term that could not be estimated, which is the model
#'     \code{config$test_method="msqrob"} fits. The gene level fit is run for every gene
#'     and used for these, so that the variance prior \code{msqrob2} shrinks toward is
#'     estimated across all of them rather than across the handful that needed it;
#'     \code{hits$fit_type} records which route each row took.
#'   \code{hits$nNonZero} is the number of observations in which the gene was measured
#'     at all, meaning in which at least one of its features has a value. It is computed
#'     from \code{state$expression} rather than read from the aggregated
#'     \code{rowData()}, for both paths: \code{QFeatures::aggregateFeatures()} keeps only
#'     the \code{rowData()} columns that are constant within a group, so a per feature
#'     count of observed values survives aggregation only when every feature of a gene
#'     was measured the same number of times and is dropped silently otherwise. On the
#'     inputs where reading it worked the two agree, a per feature count that is constant
#'     within a gene being that gene's count of observations.
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
#'     \code{msqrob2}}, which is one of the two places in this package where the reported
#'     statistic is not the engine's own; the other is the joint case of
#'     \code{h0testr::test_deqms()}, for the same reason and by the same route. The
#'     single column case is untouched and is
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
#'     \code{feat_id_col}      \cr \tab Name of column in \code{state$features} with unique feature ids; must differ from \code{config$gene_id_col} when \code{aggregate=TRUE}. \cr
#'     \code{test_random_obs}  \cr \tab Optional logical; whether the \code{aggregate=TRUE} fit includes a random observation effect alongside the random feature effect. Defaults to \code{TRUE} when absent, which is the calibrated model; see Details. Ignored when \code{aggregate=FALSE}. \cr
#'     \code{test_ridge}       \cr \tab Optional logical; whether the \code{aggregate=TRUE} fit penalizes the fixed effects. Defaults to \code{FALSE} when absent, which is \code{msqrob2}'s own default; see Details. Ignored when \code{aggregate=FALSE}. \cr
#'   }
#' @param maxit Integer scalar >= 1. How many iterations to use for \code{rlm} fitting.
#' @param aggregate Logical scalar. \code{FALSE}, the default, aggregates the features of
#'   a gene and fits the aggregated values, \code{msqrob2::msqrob()}. \code{TRUE} fits
#'   one mixed model per gene over the rows of its features instead,
#'   \code{msqrob2::msqrobAggregate()}, which is what
#'   \code{config$test_method="msqrob_agg"} selects. \code{TRUE} needs feature level
#'   input, so \code{config$feat_id_col} and \code{config$gene_id_col} must name
#'   different columns of \code{state$features}; naming one column as both is refused,
#'   the fit it describes being the one \code{aggregate=FALSE} performs.
#' @return
#'   A list with components:
#'   \tabular{ll}{
#'     \code{hits}  \cr \tab \code{data.frame} of results; columns \code{config$gene_id_col},
#'       \code{c("nNonZero", ".n")}, \code{fit_type} when \code{aggregate=TRUE}, and then
#'       either
#'       \code{c("logFC", "se", "df", "t", "pval", "adjPval")} when one design matrix
#'       column carries the test, or
#'       \code{c("f_statistic", "df1", "df2", "pval", "adjPval")} when several do. \cr
#'     \code{fit}   \cr \tab Model returned by \code{msqrob2::hypothesisTest()}, or by
#'       \code{msqrob2::msqrob()} for a joint test, which does not go through
#'       \code{msqrob2::hypothesisTest()}. With \code{aggregate=TRUE} it is what
#'       \code{msqrob2::msqrobAggregate()} returned, with the same two possibilities on
#'       top of it; the fitted models are in
#'       \code{rowData(fit[["genes"]])$msqrobModels} either way. \cr
#'   }
#'   \code{logFC} is an effect size on the scale of \code{state$expression}: for a two level
#'     factor, the difference between its levels, so a log fold change when the input is log
#'     transformed; for a \strong{continuous} covariate, the change \strong{per unit} of it,
#'     whose size depends on the units the covariate is recorded in. This is what
#'     \code{h0testr::test()} reports as \code{logfc}. The joint test reports no effect size at
#'     all, having no single contrast to report, and \code{h0testr::test()} reports the total
#'     swing instead. See \code{h0testr::test()}.
#' @examples
#' ## setup of expression data: ten peptides per gene, a third of them dropped, and no
#' ##   missing values, so that the example is about the test rather than about missingness:
#' set.seed(101)
#' samps <- h0testr::sim_samples(factors=list(grp=c("ctl", "trt"), sex=c("F", "M")),
#'   n_per_cell=3)
#' sim <- h0testr::sim_design(samps, frm=~grp + sex, test_term="grp", n_genes=100,
#'   peps_per_gene=10, p_drop=0.33, mnar_c0=-Inf, mnar_c1=0, mcar_p=0)
#' state <- sim$state
#' config <- sim$config
#' rm(samps, sim)
#'
#' ## no grp:sex term here, so testing "grp" is the single coefficient grptrt, which
#' ##   msqrob2::hypothesisTest() reports as a moderated t. Adding grp:sex would make
#' ##   it a joint test of grptrt and grptrt:sexM, reported as an F computed here:
#' out <- h0testr::initialize(state, config, minimal=TRUE)
#'
#' ## actual test:
#' result <- h0testr::test_msqrob(out$state, out$config)
#' head(result$hits)
#'
#' ## the same genes tested from their peptides instead, one mixed model each, which
#' ##   test_method="msqrob_agg" selects. A few genes are enough to show the shape:
#' keep <- out$state$features$gene_id %in% unique(out$state$features$gene_id)[1:8]
#' small <- list(expression=out$state$expression[keep, , drop=FALSE],
#'   features=out$state$features[keep, , drop=FALSE], samples=out$state$samples)
#' mixed <- h0testr::test_msqrob(small, out$config, aggregate=TRUE)
#' head(mixed$hits[, c("gene_id", "nNonZero", ".n", "fit_type", "logFC", "df", "pval")])
#'
#' ## fit_type says which genes could not carry a random feature effect and were fitted
#' ##   at the gene level instead, and the formula used is in the log:
#' table(mixed$hits$fit_type)

test_msqrob <- function(state, config, maxit=100, aggregate=FALSE) {

  check_config(config)
  f.check_state(state, config)

  ## the aggregate path models the feature level, so there has to be one. With one
  ##   feature per gene the random feature effect is a single unknown confounded with
  ##   the intercept, its variance is not identified, and the model reduces to the one
  ##   aggregate=FALSE fits, so that is what to use. Same refusal, and for the same
  ##   reason, as test_prolfqua(mixed=TRUE):

  if(aggregate && config$feat_id_col %in% config$gene_id_col) {
    f.err("test_msqrob: the aggregate path needs feature level input, so",
      "config$feat_id_col and config$gene_id_col must name different columns of",
      "state$features, and both are '", config$feat_id_col, "';", "\n",
      "with one feature per gene the random feature effect is not identified and the",
      "model reduces to the one test_method 'msqrob' fits, so use that instead, or",
      "supply un-aggregated data carrying a gene id column", config=config)
  }

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

  ## f.parse_frm()$frm, not config$frm: f.design_test_cols() below selects the
  ##   tested column from a design built on the parsed formula, whose interaction
  ##   labels have their variables sorted, and stats::model.matrix() names an
  ##   interaction column in the order the term is written, so ~sex*batch yields
  ##   sexM:batchb2 here and batchb2:sexM there. Same fit either way, but the
  ##   contrast below is matched to the fit by name:

  parsed <- f.parse_frm(config$frm, config)

  ## one gene level fit of the aggregated values, or one feature level mixed model per
  ##   gene; either way the aggregated assay is called "genes" and the models sit in
  ##   rowData(obj[["genes"]])$msqrobModels, so everything below is shared:

  fit_type <- NULL

  if(aggregate) {

    fitted <- f.msqrob_agg(obj, state, config, parsed$frm, maxit)
    obj <- fitted$obj
    fit_type <- fitted$fit_type

  } else {

    ## aggregate peptides into genes:
    obj <- QFeatures::aggregateFeatures(obj, i="features",
      fcol=config$gene_id_col, na.rm=T, name="genes", fun=base::colMeans)

    obj <- msqrob2::msqrob(object=obj, i="genes", formula=parsed$frm, maxitRob=maxit)
  }

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

  ## config$test_ridge renames every fitted fixed effect, so the design column names
  ##   are translated once here and everything below asks for the fitted names. The
  ##   translation is checked against a fit rather than trusted, since a name msqrob2
  ##   does not have is answered with NAs and not with an error:

  ridge <- aggregate && isTRUE(config$test_ridge)
  cols_fit <- f.msqrob_parms(cols_pick, ridge)

  if(ridge) {

    models <- SummarizedExperiment::rowData(obj[["genes"]])$msqrobModels
    ok <- which(!vapply(models, msqrob2::getFitMethod, character(1)) %in% "fitError")

    if(!length(ok)) {
      f.err("test_msqrob: no gene could be fit, so there is nothing to test",
        config=config)
    }

    nom <- names(msqrob2::getCoef(models[[ok[1]]]))

    if(!all(cols_fit %in% nom)) {
      f.err("test_msqrob: parameter", cols_fit, "carrying the test of",
        f.test_label(design, config), "is not among those msqrob2 fitted with",
        "config$test_ridge TRUE;", "\n", "parameters fit:",
        paste(nom, collapse=", "), config=config)
    }
  }

  ## the number of observations each gene was measured in, and the number of features
  ##   behind it. See f.msqrob_n_obs() for why the first is not read from the
  ##   aggregated rowData():

  dat <- SummarizedExperiment::rowData(obj[["genes"]])
  dat <- as.data.frame(dat[, c(config$gene_id_col, ".n"), drop=F])
  dat$nNonZero <- as.numeric(f.msqrob_n_obs(state, config)[dat[[config$gene_id_col]]])
  dat <- dat[, c(config$gene_id_col, "nNonZero", ".n"), drop=F]

  if(!is.null(fit_type)) {
    dat$fit_type <- as.character(fit_type[dat[[config$gene_id_col]]])
  }

  ## config$contrast is a weighted sum of coefficients, which is one contrast however
  ##   many it weights, so msqrob2::hypothesisTest() answers it directly and the
  ##   statistic reported is the package's own moderated t rather than the joint F
  ##   below. The weight matrix is the one f.contrast_vector() built, handed over as a
  ##   matrix rather than as text for msqrob2::makeContrast() to re-parse, so that
  ##   there is no second reading of config$contrast to disagree with the first:

  if(!is.null(design$contrast)) {

    nom <- make.names(trimws(config$contrast))

    con <- matrix(design$contrast, ncol=1,
      dimnames=list(f.msqrob_parms(colnames(design$X), ridge), nom))

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

  if(length(cols_fit) %in% 1) {

    con <- msqrob2::makeContrast(contrasts=paste0(cols_fit, "=0"),
      parameterNames=c(cols_fit))

    obj <- msqrob2::hypothesisTest(object=obj, i="genes", contrast=con,
      modelColumn="msqrobModels")

    tbl <- SummarizedExperiment::rowData(obj[["genes"]])[[cols_fit]]
    tbl <- cbind(dat[rownames(tbl), , drop=F], tbl)
    tbl <- tbl[order(tbl$adjPval, -abs(tbl$logFC)), ]

  } else {

    tbl <- f.msqrob_wald(obj, cols_fit, config)
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

    ## msqrob2 leaves the covariance of a mixed fit as a Matrix rather than as a base
    ##   matrix, msqrob2:::.noridge_msqrobLmer() not coercing it the way its ridge
    ##   counterpart does, so coerce here before subsetting by name:

    vcv <- try(as.matrix(vcv), silent=T)
    if(inherits(vcv, "try-error")) next
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

## helper for test_msqrob(): the names under which msqrob2 fitted the columns of the
##   design. msqrob2 penalizes the fixed effects, ridge=TRUE, by refitting them as a
##   random effect over a grouping factor it calls "ridge" whose levels are the design
##   columns, and lme4 names a random effect estimate by pasting the name of the
##   grouping factor in front of the level, so each column comes back as
##   "ridge<column>". The intercept is not penalized and keeps its own name; a design
##   that suppresses the intercept still gets one, unpenalized, which is simply not
##   asked about. Checked against msqrob2 1.8.0 and pinned by the tests rather than
##   assumed. With ridge=FALSE, which is msqrob2's default throughout and this
##   package's, the fitted names are the design column names unchanged:

f.msqrob_parms <- function(cols, ridge) {
  if(!isTRUE(ridge)) return(cols)
  return(ifelse(cols %in% "(Intercept)", cols, paste0("ridge", cols)))
}

## helper for test_msqrob(): the number of observations in which each gene was measured
##   at all, meaning in which at least one of its features has a value. Computed here
##   rather than read from the aggregated rowData() because
##   QFeatures::aggregateFeatures() reduces rowData() by keeping only the columns that
##   are constant within a group, so a per feature count of observed values survives
##   aggregation only when every feature of a gene was measured the same number of
##   times, and is silently dropped otherwise. Reading it was a latent error for any
##   data with feature specific missingness, which is essentially all of it. On the
##   inputs where the old lookup worked the two agree: a per feature count that is
##   constant within a gene is that gene's count of observations:

f.msqrob_n_obs <- function(state, config, who="test_msqrob") {
  genes <- f.gene_ids(state$features, config, who)
  ok <- !is.na(state$expression)
  n <- tapply(1:nrow(ok), genes, function(idxs) {
    sum(apply(ok[idxs, , drop=F], 2, any))
  })
  return(n)
}

## helper for test_msqrob(): the aggregate=TRUE fit. msqrob2::msqrobAggregate() is not
##   an alternative to aggregating so much as a different use of it: it aggregates to
##   make the assay that carries the result, one row per gene, but fits
##   msqrob2::msqrobLmer() on the un-aggregated assay with the features of a gene as one
##   group, so the models are feature level and the results gene level. The aggregated
##   values take no part in any fit, which is why base::colMeans() is passed as the
##   summary rather than msqrob2's default MsCoreUtils::robustSummary(): the latter
##   summarizes a sample whose every value in a group is exactly 0 as NA, silently, see
##   f.combine_features_robust_summary(), and there is nothing to be gained by inviting
##   that in a quantity nothing reads.
##   The random part is the feature effect msqrob2 documents plus, unless
##   config$test_random_obs is FALSE, a random observation effect. Both are needed for
##   the same reason as in f.prolfqua_mixed(): the covariates come from state$samples
##   and so vary across observations rather than across a gene's features, and without a
##   term for the observation the features of a gene stand as independent measurements of
##   it. Simulated under the null at 300 genes, 5 features per gene and 12 observations,
##   with a per-observation effect the size of the residual, rejection at the 5% level
##   was 0.360 with the feature effect alone and 0.093 with both; with no such effect
##   present, 0.063 and 0.053. See test_msqrob() for what the residual 0.093 is.
##   Returns the object with the fitted models in rowData(obj[["genes"]])$msqrobModels
##   and the route each gene took:

f.msqrob_agg <- function(obj, state, config, frm, maxit) {

  ridge <- isTRUE(config$test_ridge)
  random_obs <- is.null(config$test_random_obs) || isTRUE(config$test_random_obs)

  ran <- paste0("(1|`", config$feat_id_col, "`)")

  if(random_obs) {
    ran <- paste(ran, "+", paste0("(1|`", config$obs_col, "`)"))
  } else {
    f.msg("WARNING: test_msqrob: config$test_random_obs is FALSE, so the model carries",
      "a random feature effect but no random observation effect, which is the structure",
      "msqrob2 documents;", "\n", "the covariates vary across observations rather than",
      "across a gene's features, so without that term the features of a gene stand as",
      "independent measurements of it and the test is anti-conservative; see",
      "test_msqrob() for the simulated rejection rates", config=config)
  }

  rhs <- paste(deparse(frm[[length(frm)]]), collapse=" ")
  frm_mix <- stats::as.formula(paste("~", rhs, "+", ran))

  f.msg("test_msqrob: fitting one mixed model per gene over the rows of its features:",
    paste(deparse(frm_mix), collapse=" "), if(ridge) "with the fixed effects penalized",
    config=config)

  ## the aggregated assay is named "genes" for both paths, so that everything
  ##   downstream of the fit reads the same place:

  obj <- f.quiet_fits(msqrob2::msqrobAggregate(object=obj, i="features",
    fcol=config$gene_id_col, name="genes", formula=frm_mix, ridge=ridge, robust=T,
    maxitRob=maxit, modelColumnName="msqrobModels",
    aggregateFun=function(x, ...) base::colMeans(x, na.rm=T)),
    "mixed", config, "test_msqrob")

  dat <- SummarizedExperiment::rowData(obj[["genes"]])
  models <- dat$msqrobModels
  ids <- as.character(dat[[config$gene_id_col]])
  fit_type <- vapply(models, msqrob2::getFitMethod, character(1))

  ## genes with fewer than two observed features are fitted at the gene level instead:
  ##   their random feature effect would be a single unknown confounded with the
  ##   intercept, so its variance is not identified and lme4 refuses the fit, which
  ##   msqrob2 records as a fitError and would report as a row of NAs. Such genes are
  ##   common, and dropping them would leave a gene table that does not line up with
  ##   what another test_method produces, so they get the same fixed model without the
  ##   term that could not be estimated, which is what test_method="msqrob" fits.
  ##   Fitted for every gene and used for these, rather than fitted for these alone,
  ##   so that the variance prior msqrob2 shrinks toward is estimated across all of
  ##   them; hits$fit_type records which route each row took:

  n_feat <- tapply(rowSums(!is.na(state$expression)) > 0,
    f.gene_ids(state$features, config, "test_msqrob"), sum)
  thin <- names(n_feat)[n_feat < 2]
  need <- which(fit_type %in% "fitError" & ids %in% thin)

  if(length(need)) {

    f.msg("test_msqrob:", length(need), "gene(s) have fewer than two observed",
      "features, so the random feature effect is not identified and lme4 refuses the",
      "fit; fitting those at the gene level instead, which is the same fixed model",
      "without the term that could not be estimated;", "\n", "affected genes:",
      paste(utils::head(ids[need], 5), collapse=", "),
      if(length(need) > 5) paste("and", length(need) - 5, "more") else "",
      config=config)

    obj2 <- try(f.quiet_fits(msqrob2::msqrob(object=obj, i="genes", formula=frm,
      ridge=ridge, robust=T, maxitRob=maxit, modelColumnName="msqrobModelsGene"),
      "gene level", config, "test_msqrob"), silent=T)

    if(inherits(obj2, "try-error")) {
      f.msg("WARNING: test_msqrob: the gene level fit for those genes failed, so they",
        "are reported as NA;", "\n", "  ",
        trimws(conditionMessage(attr(obj2, "condition"))), config=config)
    } else {
      m2 <- SummarizedExperiment::rowData(obj2[["genes"]])$msqrobModelsGene
      for(idx in need) {
        models[[idx]] <- m2[[idx]]
        fit_type[idx] <- msqrob2::getFitMethod(m2[[idx]])
      }
      SummarizedExperiment::rowData(obj[["genes"]])[["msqrobModels"]] <- models
    }
  }

  return(list(obj=obj, fit_type=stats::setNames(fit_type, ids)))
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
#'   \code{diff} is an effect size on the scale of \code{state$expression}: for a two level
#'     factor, the difference between its levels, so a log fold change when the input is log
#'     transformed; for a \strong{continuous} covariate, the change \strong{per unit} of it,
#'     whose size depends on the units the covariate is recorded in. This is what
#'     \code{h0testr::test()} reports as \code{logfc}. The likelihood ratio test reports no
#'     effect size at all, having no single contrast to report, and \code{h0testr::test()}
#'     reports the total swing instead. See \code{h0testr::test()}.
#' @examples
#' ## setup of expression data: ten peptides per gene, a third of them dropped, and no
#' ##   missing values, so that the example is about the test rather than about missingness:
#' set.seed(101)
#' samps <- h0testr::sim_samples(factors=list(grp=c("ctl", "trt"), sex=c("F", "M")),
#'   n_per_cell=3)
#' sim <- h0testr::sim_design(samps, frm=~grp + sex, test_term="grp", n_genes=100,
#'   peps_per_gene=10, p_drop=0.33, mnar_c0=-Inf, mnar_c1=0, mcar_p=0)
#' state <- sim$state
#' config <- sim$config
#' rm(samps, sim)
#'
#' ## no grp:sex term here, so testing "grp" is the single coefficient grptrt, and a
#' ##   fold change is reported. Adding one would, by marginality, make it a joint
#' ##   test of grptrt and grptrt:sexM, run as a likelihood ratio test reporting an F
#' ##   statistic and no fold change:
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
        "use test_method 'lm', 'trend', 'voom', 'deqms', 'msqrob', 'msqrob_agg',",
        "'prolfqua' or 'prolfqua_lmer' to test against zero",
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

## The rows of L for f.wald_f(): the hypothesis as a matrix over the coefficients of
##   the design, which is the identity rows of the tested columns for config$test_term
##   and the single row of weights for config$contrast. One shape for both, so that
##   nothing downstream branches on which was asked for. Named by the make.names()
##   forms of the design matrix column names, those being the names the fits carry:

f.test_L <- function(design) {

  cols <- make.names(colnames(design$X), unique=T)

  L <- if(is.null(design$contrast)) {
    diag(length(cols))[design$cols_test, , drop=F]
  } else {
    matrix(design$contrast, nrow=1)
  }

  colnames(L) <- cols

  return(L)
}

## The Wald F-test of L %*% beta == 0 from one fitted model. For a mixed model that is
##   what lmerTest::contest() computes, an F whose denominator degrees of freedom come
##   from the Satterthwaite approximation, so that a covariate which varies across
##   observations is tested against between-observation variation rather than against
##   the feature by observation residual; for a least squares fit it is the ordinary
##   F, which is also the nested model comparison f.nested_f() performs. Both are
##   needed because the mixed path falls back to stats::lm() for a gene with a single
##   feature, and both come back in one shape so that the two kinds of row sit in one
##   table.
##   Differencing the residual sums of squares of two mixed fits, which is what
##   f.prolfqua_nested_f() does for the least squares path, is not available here: the
##   variance components are re-estimated for each fit, so sigma^2 times the residual
##   degrees of freedom is not a residual sum of squares that decomposes, and the
##   difference of two of them is not an F. Hence a Wald test of one fit rather than a
##   comparison of two.
##   NULL when the fit cannot support the test, which the caller reports as a dropped
##   gene: a coefficient that lme4 dropped for rank deficiency or that least squares
##   left aliased is absent or NA, and the hypothesis is then not estimable for that
##   gene. A singular fit, meaning a variance component estimated at zero, is not such
##   a case: the F is still the right test, and the model has simply collapsed toward
##   the one without that effect:

f.wald_f <- function(fit, L) {

  if(inherits(fit, "merMod")) {

    beta <- try(lme4::fixef(fit), silent=T)
    if(inherits(beta, "try-error") || is.null(names(beta))) return(NULL)
    if(!all(colnames(L) %in% names(beta))) return(NULL)

    tst <- try(lmerTest::contest(fit, L[, names(beta), drop=F], joint=T), silent=T)
    if(inherits(tst, "try-error") || !nrow(tst)) return(NULL)

    return(data.frame(Df=tst[["NumDF"]], Sum.Sq=tst[["Sum Sq"]],
      Mean.Sq=tst[["Mean Sq"]], F.value=tst[["F value"]], p.value=tst[["Pr(>F)"]],
      df.denom=tst[["DenDF"]]))
  }

  beta <- try(stats::coef(fit), silent=T)
  if(inherits(beta, "try-error") || is.null(names(beta))) return(NULL)
  if(!all(colnames(L) %in% names(beta))) return(NULL)

  nom <- colnames(L)
  if(any(is.na(beta[nom]))) return(NULL)

  ## the Wald statistic is the extra sum of squares divided by the error variance, so
  ##   multiplying it back by that variance recovers the sum of squares the nested
  ##   comparison would report, and at one degree of freedom the whole thing is the
  ##   square of the t-statistic:

  V <- try(stats::vcov(fit)[nom, nom, drop=F], silent=T)
  if(inherits(V, "try-error")) return(NULL)

  Lb <- L %*% beta[nom]
  w <- try(drop(t(Lb) %*% solve(L %*% V %*% t(L), Lb)), silent=T)
  if(inherits(w, "try-error") || !is.finite(w)) return(NULL)

  df_num <- qr(L)$rank
  df_den <- stats::df.residual(fit)
  if(!is.finite(df_den) || df_den <= 0 || df_num < 1) return(NULL)

  s2 <- stats::sigma(fit)^2
  fval <- w / df_num

  return(data.frame(Df=df_num, Sum.Sq=w * s2, Mean.Sq=w * s2 / df_num,
    F.value=fval, p.value=stats::pf(fval, df_num, df_den, lower.tail=F),
    df.denom=df_den))
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
##   FALSE, against a flat prior unless trend is TRUE and a covariate to
##   fit the prior against is supplied, keyed by config$feat_id_col so that it can be
##   subset to the features that survive the drop below.
##   Features whose data do not support the requested test are dropped explicitly
##   and reported: a rank shortfall means the reduced fit lost fewer columns than
##   the test needs, which happens when missingness leaves a covariate aliased for
##   that feature. Worth being loud about, because both silent alternatives are
##   worse: stats::anova() omits the term for such a feature so it used to vanish
##   from the results unannounced, and the comparison below would otherwise report
##   it with a zero numerator degrees of freedom and an NaN F:

f.prolfqua_nested_f <- function(fit_full, fit_red, design, config, covariate=NULL,
    trend=NULL) {

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

  ## the caller resolves this, config$test_trend answering when it passes nothing, so
  ##   that a direct call behaves as before; see f.is_trend():

  trend <- f.is_trend(trend, config, "f.prolfqua_nested_f")

  if(moderate) {

    ## subset after the drop above, so that features the test cannot be run for do
    ##   not contribute to the prior either. A name that is not in the covariate
    ##   gives NA, which f.moderate_var() reports and falls back from:

    cov <- NULL

    if(trend) {
      if(is.null(covariate)) {
        f.err("f.prolfqua_nested_f: the trended prior was asked for but no covariate",
          "was supplied to fit the prior variance against", config=config)
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
      "unmoderated F-test", if(trend) paste("and ignoring the trended prior, which",
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

## The fixed effects of a fit from either strategy, under the names test_prolfqua()
##   gave the design matrix columns. stats::coef() of a mixed model returns the per
##   group coefficients as a list rather than the fixed effects, hence the branch:

f.fixed_coefs <- function(fit) {
  if(inherits(fit, "merMod")) return(lme4::fixef(fit))
  return(stats::coef(fit))
}

## Fit a collection of models with what they have to say about themselves collected
##   rather than emitted. A mixed model reports non-convergence and a fit at the
##   boundary of the parameter space per model, which over a few thousand genes is a few
##   thousand lines saying something already recorded per gene in hits$isSingular, so
##   the count and the first message go to the log and the rest is muffled. Warnings are
##   counted rather than passed through because they arrive from inside
##   prolfqua::build_model()'s dplyr::mutate(), which relabels them as its own, and from
##   inside msqrob2::msqrobLmer()'s BiocParallel::bplapply(), which does the same.
##   Shared by test_prolfqua(mixed=TRUE) and test_msqrob(aggregate=TRUE), hence who:

f.quiet_fits <- function(expr, what, config, who="test_prolfqua") {

  n <- 0L
  first <- NULL

  out <- withCallingHandlers(suppressMessages(expr),
    warning=function(w) {
      n <<- n + 1L
      if(is.null(first)) first <<- conditionMessage(w)
      invokeRestart("muffleWarning")
    })

  if(n) {
    f.msg(paste0("WARNING: ", who, ":"), n, "warning(s) raised while fitting the",
      what, "models; the first was:", "\n", "  ", trimws(first), config=config)
  }

  return(out)
}

## The mixed path of test_prolfqua(): one model per gene over the long table of its
##   features, giving gene level inference from feature level data without aggregating.
##   Takes the LFQData object, the design and the make.names() forms of its column
##   names that test_prolfqua() has already put into the long table, so that the fixed
##   part of the model is the same design every other engine here is handed.
##   The random effects are a random intercept per feature and a random intercept per
##   observation.
##   The feature effect absorbs the feature baselines, which differ by orders of
##   magnitude among the precursors of one protein. Random rather than a fixed effect
##   per feature: it costs one variance component instead of one coefficient per
##   feature, and shrinks a thinly observed feature toward the gene mean rather than
##   estimating it from its own few observations, which is what makes the fit unbiased
##   under feature specific missingness. Aggregating the features of a gene first
##   instead compares a different subset of them in each observation.
##   The observation effect is the error stratum the fixed effects belong to. Every
##   covariate of config$frm comes from state$samples, so it varies across observations
##   and not within one; without this term the features of a gene are asserted to be
##   independent measurements of that observation and the Satterthwaite denominator
##   degrees of freedom come out near the number of feature by observation rows rather
##   than near the number of observations. Simulated under the null at 12 observations,
##   5 features per gene and a per observation effect the size of the residual, that
##   rejects 26% of the time at the 5% level, against 4.7% with the term present. Where
##   the per observation variance really is zero the term costs a little conservatism,
##   3.0% against 5.3%, which is the safe direction. config$test_random_obs=FALSE drops
##   it, giving the feature only structure prolfqua documents, for comparison.
##   Not included is a random slope: letting each feature respond differently to the
##   tested covariate makes the fixed effect an average over a distribution of feature
##   specific responses, which is a different hypothesis, and it is rarely identifiable
##   at proteomics sample sizes:

f.prolfqua_mixed <- function(obj, design, cols, config, trend=NULL) {

  dat <- obj$data

  ran <- paste0("(1|`", config$feat_id_col, "`)")
  random_obs <- is.null(config$test_random_obs) || isTRUE(config$test_random_obs)

  if(random_obs) {
    ran <- paste(ran, "+ (1|sample)")
  } else {
    f.msg("WARNING: test_prolfqua: config$test_random_obs is FALSE, so the fit carries",
      "no random observation effect and the features of a gene are treated as",
      "independent measurements of it;", "\n", "the test is anti-conservative to the",
      "extent that the features of a gene co-vary within an observation, which under",
      "simulation was a rejection rate of 26% at the 5% level", config=config)
  }

  ## config$test_moderate and config$test_trend describe a shrinkage of the per feature
  ##   error variance across features, which this path does not perform: its denominator
  ##   is a Satterthwaite combination of variance components rather than one residual
  ##   variance, so there is nothing single to shrink. Said out loud rather than left
  ##   implicit, because both keys are consulted by the least squares path and their
  ##   being ignored here is otherwise invisible. A note rather than a warning:
  ##   new_config() sets test_moderate=TRUE, so nothing is wrong with a default config
  ##   arriving here, and flagging every default run as a warning would be noise:

  ignored <- character(0)
  if("test_moderate" %in% names(config) && isTRUE(config$test_moderate)) {
    ignored <- c(ignored, "config$test_moderate")
  }
  ## named by what is actually set rather than by which of the two the caller used: the
  ##   value arrives resolved, so a TRUE that came from config$test_trend is
  ##   indistinguishable here from one passed to test_prolfqua(), and naming the key when
  ##   the key is what says TRUE points at the thing to change either way:

  if(f.is_trend(trend, config, "f.prolfqua_mixed")) {
    ignored <- c(ignored, if(isTRUE(config$test_trend)) "config$test_trend" else "trend")
  }

  if(length(ignored)) {
    f.msg("NOTE: test_prolfqua:", paste(ignored, collapse=" and "),
      if(length(ignored) > 1) "are" else "is", "TRUE, but the mixed path does not",
      "moderate the error variance across genes and does not consult",
      if(length(ignored) > 1) "them" else "it", ";", "\n",
      "the denominator of a Satterthwaite F is a combination of variance components",
      "rather than one residual variance, so there is no single quantity to shrink;",
      "the reported test is unmoderated and hits$moderated is FALSE", config=config)
  }

  ## genes with a single observed feature are fitted by least squares: lme4 refuses a
  ##   grouping factor with one level, correctly, the variance of that effect having no
  ##   between-feature contrast to be estimated from. Counted over the observed values,
  ##   since a feature missing everywhere is dropped by the fit and so contributes no
  ##   level:

  genes <- as.character(dat[[config$gene_id_col]])
  ok <- !is.na(dat$intensity)

  n_feat <- tapply(as.character(dat[[config$feat_id_col]])[ok], genes[ok],
    function(v) length(unique(v)))

  multi <- names(n_feat)[n_feat > 1]
  single <- names(n_feat)[n_feat %in% 1]
  gone <- setdiff(unique(genes), names(n_feat))

  if(length(gone)) {
    f.msg("WARNING: test_prolfqua: dropping", length(gone), "of",
      length(unique(genes)), "genes with no observed value in any feature;", "\n",
      "first few:", paste(utils::head(gone, 5), collapse=", "), config=config)
  }

  if(!length(multi) && !length(single)) {
    f.err("test_prolfqua: no gene has an observed value", config=config)
  }

  frm_fix <- paste("intensity ~ 0 +", paste(cols, collapse=" + "))
  frm_mix <- paste(frm_fix, "+", ran)

  f.msg("test_prolfqua:", f.test_label(design, config), "; mixed:", frm_mix, "\n",
    "  genes with several observed features:", length(multi), "; with one:",
    length(single), "; numerator df:", design$df_intend, config=config)

  mods <- list()

  if(length(multi)) {
    mods[[length(mods) + 1]] <- list(fit_type="lmer",
      model=f.quiet_fits(prolfqua::build_model(
        data=dat[genes %in% multi, , drop=F],
        model_strategy=prolfqua::strategy_lmer(frm_mix),
        subject_Id=config$gene_id_col), "lmer", config))
  }

  if(length(single)) {
    f.msg("test_prolfqua: fitting", length(single), "gene(s) with a single observed",
      "feature by least squares instead:", frm_fix, config=config)
    mods[[length(mods) + 1]] <- list(fit_type="lm",
      model=f.quiet_fits(prolfqua::build_model(
        data=dat[genes %in% single, , drop=F],
        model_strategy=prolfqua::strategy_lm(frm_fix),
        subject_Id=config$gene_id_col), "lm", config))
  }

  out <- f.prolfqua_mixed_f(mods, design, config, "test_prolfqua")

  f.msg("tested", nrow(out$hits), "genes; found",
    sum(out$hits$FDR < 0.05, na.rm=T), "hits", config=config)

  ## the two collections of fits are returned separately and named for what they are,
  ##   rather than as one object, because they were built by different strategies and a
  ##   caller reading a fit needs to know which. Either can be NULL, when no gene took
  ##   that route. The coefficients are assembled where the fits are, since they come
  ##   from both:

  types <- vapply(mods, function(m) m$fit_type, character(1))

  return(list(hits=out$hits, coefs=out$coefs,
    fit=if("lmer" %in% types) mods[[which(types %in% "lmer")[1]]]$model else NULL,
    fit_lm=if("lm" %in% types) mods[[which(types %in% "lm")[1]]]$model else NULL,
    fit_reduced=NULL, design=design))
}

## Per-gene Wald F-tests from the mixed path of test_prolfqua(), where one model is
##   fitted per gene over its features rather than one model per feature. Takes the
##   models as a list of prolfqua Model objects each labelled with how it was fitted,
##   because a gene with a single feature is fitted by least squares: its random
##   feature effect would be a single unknown confounded with the intercept, so its
##   variance is not identified and lme4 refuses the fit outright. That leaves the
##   model those genes get identical to the mixed model minus a term that could not
##   have been estimated, and the fit_type column of the result records which route
##   each row took.
##   Genes whose data do not support the test are dropped and reported, on the same
##   terms as f.prolfqua_nested_f() drops features: a fit that failed, or one where the
##   test does not come out at the intended degrees of freedom, which happens when
##   missingness leaves a tested column aliased for that gene.
##   The error variance is not moderated across genes here. The denominator of a
##   Satterthwaite F is a combination of variance components rather than one residual
##   variance, so there is no single quantity to shrink; s2.denom reports the effective
##   denominator the F actually used, being Mean.Sq divided by F.value, and the columns
##   describing a moderation are NA:

f.prolfqua_mixed_f <- function(mods, design, config, caller="f.prolfqua_mixed_f") {

  L <- f.test_L(design)
  cols <- colnames(design$X)[design$cols_test]
  nom <- make.names(colnames(design$X), unique=T)[design$cols_test]

  rows <- list()
  betas <- list()
  ids_out <- character(0)
  bad <- character(0)
  n_in <- 0
  n_by_type <- integer(0)

  for(mod in mods) {

    mdf <- as.data.frame(mod$model$modelDF)
    if(!nrow(mdf)) next

    ids <- as.character(mdf[[config$gene_id_col]])
    n_in <- n_in + nrow(mdf)
    n_by_type[mod$fit_type] <- nrow(mdf)

    for(idx in seq_len(nrow(mdf))) {

      fit <- mdf$linear_model[[idx]]
      stat <- if(is.null(fit) || is.character(fit)) NULL else f.wald_f(fit, L)

      if(is.null(stat) || !(stat$Df %in% design$df_intend)) {
        bad <- c(bad, ids[idx])
        next
      }

      beta <- try(f.fixed_coefs(fit), silent=T)
      if(inherits(beta, "try-error")) beta <- stats::setNames(numeric(0), character(0))

      rows[[length(rows) + 1]] <- data.frame(
        isSingular=mdf$isSingular[idx], nrcoef=mdf$nrcoef[idx],
        factor=config$test_term, stat, fit_type=mod$fit_type, stringsAsFactors=F)

      betas[[length(betas) + 1]] <- beta[match(nom, names(beta))]
      ids_out <- c(ids_out, ids[idx])
    }
  }

  if(length(bad)) {
    f.msg("WARNING:", caller, ": dropping", length(bad), "of", n_in,
      "genes whose data do not support the test of",
      if(is.null(design$contrast)) {
        paste0("config$test_term '", config$test_term, "';")
      } else {
        paste0("config$contrast '", trimws(config$contrast), "';")
      }, "\n", "the model could not be fitted for these genes, or the test does not",
      "come out at", design$df_intend, "degree(s) of freedom for them, so the",
      "hypothesis is not estimable;", "\n", "first few:",
      paste(utils::head(bad, 5), collapse=", "), config=config)
  }

  if(!length(rows)) {
    f.err(caller, ": no gene supports the test of",
      if(is.null(design$contrast)) {
        paste0("config$test_term '", config$test_term, "'")
      } else {
        paste0("config$contrast '", trimws(config$contrast), "'")
      }, ";", "genes considered:", n_in, "; design columns carrying the test:",
      length(design$cols_test), "; degrees of freedom intended:", design$df_intend,
      config=config)
  }

  out <- do.call(rbind, rows)
  out[[config$gene_id_col]] <- ids_out

  ## the columns of the least squares path that a Wald F on a mixed fit has no value
  ##   for, carried as NA rather than dropped so that a result table reads the same
  ##   whichever prolfqua path produced it:

  out$moderated <- FALSE
  out$trend <- FALSE
  out$s2.denom <- out$Mean.Sq / out$F.value
  out$df.prior <- as.numeric(NA)
  out$s2.prior <- as.numeric(NA)
  out$F.value.unmod <- as.numeric(NA)
  out$p.value.unmod <- as.numeric(NA)
  out$FDR <- stats::p.adjust(out$p.value, method="BH")

  keep <- c(config$gene_id_col, "isSingular", "nrcoef", "factor", "Df", "Sum.Sq",
    "Mean.Sq", "F.value", "p.value", "fit_type", "moderated", "trend", "s2.denom",
    "df.denom", "df.prior", "s2.prior", "F.value.unmod", "p.value.unmod", "FDR")
  out <- out[, keep, drop=F]

  coefs <- do.call(rbind, betas)
  dimnames(coefs) <- list(ids_out, cols)

  f.msg(caller, ": tested", nrow(out), "of", n_in, "genes; fitted",
    paste(paste0(n_by_type, " by ", names(n_by_type)), collapse=", "),
    "; singular fits:", sum(out$isSingular %in% TRUE), config=config)

  out <- out[order(out$p.value, decreasing=F), , drop=F]
  rownames(out) <- NULL

  return(list(hits=out, coefs=coefs))
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
#'     unless the \code{trend} argument, or \code{config$test_trend} when it is not
#'     given, is \code{TRUE}, which fits it against the mean
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
#'   With \code{mixed=TRUE}, which \code{test_method="prolfqua_lmer"} selects, the
#'     engine changes from \code{prolfqua::strategy_lm()} to
#'     \code{prolfqua::strategy_lmer()} and one model is fitted per gene over the rows
#'     of all of its features, rather than one model per feature. The result is gene
#'     level inference from feature level data with no aggregation step: the fixed part
#'     of the model is the same design matrix, whose columns come from
#'     \code{state$samples} and so are identical across the features of a gene, which
#'     makes its coefficients gene level quantities. \code{hits} then has one row per
#'     gene, keyed by \code{config$gene_id_col}, and \code{h0testr::run()} and
#'     \code{h0testr::tune()} do not call \code{combine_features()} before it, as they
#'     do not for \code{test_method} \code{"deqms"} and \code{"msqrob"}. Handing it
#'     data that have already been aggregated is an error rather than a silent
#'     degenerate fit; see below.
#'   The model is \code{intensity ~ 0 + <design columns> + (1|<feat_id_col>) +
#'     (1|<obs_col>)}. The random feature intercept absorbs the feature baselines,
#'     which differ by orders of magnitude among the precursors of one protein. It is
#'     random rather than a fixed effect per feature because it then costs one variance
#'     component instead of one coefficient per feature, and because a thinly observed
#'     feature is shrunk toward the gene mean rather than estimated from its own few
#'     observations, which is what makes the fit unbiased under feature specific
#'     missingness: aggregating first instead compares a different subset of a gene's
#'     features in each observation.
#'   The random observation intercept is the error stratum the fixed effects belong to.
#'     Every covariate of \code{config$frm} comes from \code{state$samples}, so it
#'     varies across observations and not within one, and without this term the
#'     features of a gene are asserted to be independent measurements of that
#'     observation. The Satterthwaite denominator degrees of freedom then come out near
#'     the number of feature by observation rows instead of near the number of
#'     observations, and the test rejects far more often than it should: simulated under
#'     the null at 12 observations, 5 features per gene and a per observation effect the
#'     size of the residual, 26 percent at the 5 percent level, against 4.7 percent with
#'     the term present. Where the per observation variance really is zero the term
#'     costs a little conservatism, 3.0 percent against 5.3 percent, which is the safe
#'     direction. \code{config$test_random_obs=FALSE} drops it, giving the feature only
#'     structure \code{prolfqua} documents, for comparison; a warning says so. A random
#'     slope is not offered: letting each feature respond differently to the tested
#'     covariate makes the fixed effect an average over a distribution of feature
#'     specific responses, which is a different hypothesis.
#'   The mixed test is a Wald F of one fit rather than a comparison of two, computed by
#'     \code{lmerTest::contest()} from the same coefficients every other engine here
#'     selects, so \code{config$test_term} still covers every term containing it and
#'     \code{config$contrast} is still one degree of freedom. Its denominator degrees of
#'     freedom are the Satterthwaite approximation, reported in \code{df.denom}, and are
#'     fractional in general. Differencing the residual sums of squares of a full and a
#'     reduced fit, which is what the least squares path does, is not available: the
#'     variance components are re-estimated for each fit, so the difference is not an F.
#'   The error variance is not moderated across genes on this path, since the
#'     denominator of a Satterthwaite F is a combination of variance components rather
#'     than one residual variance and there is no single quantity to shrink.
#'     \code{hits$moderated} is \code{FALSE}, \code{df.prior}, \code{s2.prior},
#'     \code{F.value.unmod} and \code{p.value.unmod} are \code{NA}, and
#'     \code{s2.denom} reports the effective denominator variance the F did use.
#'     \code{config$test_moderate} and the trended prior are not consulted, and
#'     a note says so when either is present and \code{TRUE}. A note rather than a
#'     warning because \code{h0testr::new_config()} sets \code{test_moderate=TRUE}, so
#'     it applies to a default configuration and describes nothing wrong.
#'   A gene with a single observed feature is fitted by \code{stats::lm()} instead: its
#'     random feature effect would be a single unknown confounded with the intercept, so
#'     that variance is not identified and \code{lme4} refuses the fit outright, and
#'     what is left is the mixed model without a term that could not have been
#'     estimated. Such genes are reported rather than dropped, since they are common and
#'     dropping them would leave a gene table that does not line up with the one another
#'     \code{test_method} produces; the \code{fit_type} column of \code{hits} records
#'     which route each row took. A gene with no observed value at all, and a gene for
#'     which the test does not come out at the intended degrees of freedom, are dropped
#'     with a warning naming how many and the first few.
#'   A fit at the boundary of the parameter space, meaning a variance component
#'     estimated at zero, is kept and flagged in \code{hits$isSingular}: the F is still
#'     the right test and the model has simply collapsed toward the one without that
#'     effect. Messages and warnings from the fits are counted and summarized in the log
#'     rather than emitted one per gene.
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
#'            \code{trend}, one fitted against mean feature intensity. \cr
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
#'     \code{test_trend}            \cr \tab Optional logical; whether the prior of that shrinkage is fitted against mean feature intensity rather than flat. Answers when the \code{trend} argument is not given; defaults to \code{FALSE} when both are absent. Unrelated to \code{config$test_method="trend"}. \cr
#'     \code{test_random_obs}       \cr \tab Optional logical; whether the \code{mixed=TRUE} fit includes a random observation effect alongside the random feature effect. Defaults to \code{TRUE} when absent, which is the calibrated model; see Details. Ignored when \code{mixed=FALSE}. \cr
#'     \code{feat_id_col}           \cr \tab Name of column in \code{state$features} with unique feature ids; must differ from \code{config$gene_id_col} when \code{mixed=TRUE}. \cr
#'     \code{normalization_method}  \cr \tab If present and \code{is_log_transformed} unset, used to infer it. \cr
#'   }
#' @param is_log_transformed Logical scalar: whether \code{state$expression} has
#'   been log transformed. Defaults to \code{config$is_log_transformed}, which
#'   \code{h0testr::initialize()} and \code{h0testr::normalize()} maintain;
#'   passing both is an error unless they agree.
#' @param mixed Logical scalar: whether to fit one mixed model per gene over the rows
#'   of its features, with the feature and the observation as random effects, instead
#'   of one least squares model per feature. \code{TRUE} is what
#'   \code{config$test_method="prolfqua_lmer"} selects, gives one result row per gene
#'   from feature level input, and requires \code{config$feat_id_col} and
#'   \code{config$gene_id_col} to name different columns. See Details.
#' @param trend Logical scalar: whether the variance prior of the moderation is fitted
#'   against mean feature intensity, which is \code{limma::eBayes(trend=TRUE)}'s
#'   covariate, rather than being flat. Defaults to \code{config$test_trend}, and to
#'   \code{FALSE} when that is absent; an argument that disagrees with the configuration
#'   wins. Not consulted when \code{mixed=TRUE}, whose denominator is a combination of
#'   variance components rather than one residual variance, and a \code{TRUE} there is
#'   reported as a \code{NOTE}.
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
#'   With \code{mixed=TRUE} the components differ: \code{hits} has one row per gene,
#'     keyed by \code{config$gene_id_col} alone, with \code{df.denom} holding the
#'     Satterthwaite denominator degrees of freedom, an extra \code{fit_type} column
#'     recording whether that gene was fitted by \code{lmer} or fell back to \code{lm},
#'     and the columns describing a moderation \code{NA}; \code{coefs} is a matrix of
#'     the tested coefficients, one row per gene, assembled from both collections of
#'     fits; \code{fit} holds the mixed models and \code{fit_lm} the least squares
#'     models for single feature genes, either of which is \code{NULL} when no gene took
#'     that route; and \code{fit_reduced} is \code{NULL}, no reduced model being fitted.
#'   Note the granularity of \code{hits} follows the input: called with
#'     \code{mixed=FALSE} on peptide-level data, as in the example below, it gives one
#'     row per peptide, while in the \code{h0testr::run()} and \code{h0testr::tune()}
#'     pipelines \code{combine_features()} has already run, so it gives one row per
#'     gene/protein group, like every other test method. With \code{mixed=TRUE} it is
#'     gene level whatever the input, that being the point of the mixed model.
#'   No effect size is reported here: the test is an F test of the columns carrying
#'     \code{config$test_term}, whatever their number, so there is no coefficient to report
#'     alongside it. \code{h0testr::test()} fills its \code{logfc} column from the fitted
#'     coefficients instead, taken from \code{fit$modelDF} or, with \code{mixed=TRUE}, from
#'     \code{coefs}: the signed coefficient when one column carries the test, which for a
#'     \strong{continuous} covariate is a change \strong{per unit} of it rather than a fold
#'     change between groups, and the total swing when several do. See \code{h0testr::test()}.
#' @examples
#' ## setup of expression data: ten peptides per gene, a third of them dropped, and no
#' ##   missing values, so that the example is about the test rather than about missingness:
#' set.seed(101)
#' samps <- h0testr::sim_samples(factors=list(grp=c("ctl", "trt"), sex=c("F", "M")),
#'   n_per_cell=3)
#' sim <- h0testr::sim_design(samps, frm=~grp + sex, test_term="grp", n_genes=100,
#'   peps_per_gene=10, p_drop=0.33, mnar_c0=-Inf, mnar_c1=0, mcar_p=0)
#' state <- sim$state
#' config <- sim$config
#' rm(samps, sim)
#'
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
#'
#' ## the mixed path, which test_method="prolfqua_lmer" selects: one model per gene
#' ##   over the rows of its peptides, so the result is gene level although the input
#' ##   is not. A few genes are enough to show the shape:
#' config$test_moderate <- NULL
#' config$test_trend <- NULL
#' keep <- state$features$gene_id %in% unique(state$features$gene_id)[1:8]
#' small <- list(expression=state$expression[keep, , drop=FALSE],
#'   features=state$features[keep, , drop=FALSE], samples=state$samples)
#' out <- h0testr::initialize(small, config, minimal=TRUE)
#' mix <- h0testr::test_prolfqua(out$state, out$config, is_log_transformed=FALSE,
#'   mixed=TRUE)
#' nrow(mix$hits)                                  ## one row per gene, not per peptide
#' head(mix$hits[, c("gene_id", "Df", "F.value", "p.value", "df.denom", "fit_type")])
#'
#' ## df.denom is the Satterthwaite denominator, which the random observation effect
#' ##   pulls down toward the number of observations; dropping that effect gives
#' ##   prolfqua's documented peptide-only structure and a much larger denominator,
#' ##   which is why it is not the default:
#' out$config$test_random_obs <- FALSE
#' peponly <- h0testr::test_prolfqua(out$state, out$config,
#'   is_log_transformed=FALSE, mixed=TRUE)
#' range(mix$hits$df.denom)
#' range(peponly$hits$df.denom)

test_prolfqua <- function(state, config, is_log_transformed=NULL, mixed=FALSE,
    trend=NULL) {

  is_log_transformed <- f.is_log_transformed(is_log_transformed, config,
    "test_prolfqua")

  ## the argument overrides config$test_trend, which is what it answers from when not
  ##   given; see f.is_trend(). Resolved here and threaded down rather than re-read
  ##   further in, so that one value describes the whole call:

  trend <- f.is_trend(trend, config, "test_prolfqua")

  parsed <- f.parse_frm(config$frm, config)

  ## the mixed path models the features of a gene instead of aggregating them, so it
  ##   needs both levels present and distinct. With one feature per gene the random
  ##   feature effect is a single unknown confounded with the intercept, so its variance
  ##   is not identified and lme4 refuses the fit; what remains is the model the least
  ##   squares path already fits:

  if(mixed && config$feat_id_col %in% config$gene_id_col) {
    f.err("test_prolfqua: the mixed path needs feature level input, so",
      "config$feat_id_col and config$gene_id_col must name different columns of",
      "state$features, and both are '", config$feat_id_col, "';", "\n",
      "with one feature per gene the random feature effect is not identified and the",
      "model reduces to the one test_method 'prolfqua' fits, so use that instead, or",
      "supply un-aggregated data carrying a gene id column", config=config)
  }

  idvars <- unique(c(config$gene_id_col, config$feat_id_col))

  ## the ids the mixed path groups by are the ones test() reports its rows under, so
  ##   they are resolved here the same way, giving a feature with no gene id of its own
  ##   an id rather than pooling every such feature into one group:

  if(mixed) {
    state$features[[config$gene_id_col]] <- f.gene_ids(state$features, config,
      "test_prolfqua")
  }
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

  cols_red <- character(0)

  if(mixed) {

    ## the mixed path tests one fit by Wald rather than comparing two, so it forms no
    ##   reduced design; see f.wald_f() for why differencing two mixed fits is not
    ##   available:

  } else if(is.null(design$contrast)) {

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

  ## the mixed path takes over here: one model per gene over the rows of its features,
  ##   with the feature and the observation entering as random effects, rather than one
  ##   model per feature. Everything above is shared, the long table and the design
  ##   being the same either way:

  if(mixed) return(f.prolfqua_mixed(obj, design, cols, config, trend=trend))

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

  if(trend) {
    covariate <- rowMeans(state$expression, na.rm=T)
    names(covariate) <- as.character(state$features[[config$feat_id_col]])

    if(!is_log_transformed) {
      f.msg("WARNING: test_prolfqua: the trended prior was asked for but the response",
        "is not log transformed, so the prior variance is being fitted against mean",
        "untransformed intensity;", "\n", "the mean-variance trend",
        "limma::eBayes(trend=TRUE) models is a trend in log intensity", config=config)
    }
  }

  tbl <- f.prolfqua_nested_f(fit_full, fit_red, design, config, covariate=covariate,
    trend=trend)

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
#'   \code{logFC} is an effect size on the scale of \code{state$expression}: for a two level
#'     factor, the difference between its levels, so a log fold change when the input is log
#'     transformed; for a \strong{continuous} covariate, the change \strong{per unit} of it,
#'     whose size depends on the units the covariate is recorded in. This is what
#'     \code{h0testr::test()} reports as \code{logfc}. A joint test has no \code{logFC}
#'     column at all, having no single fold change to report, and \code{h0testr::test()}
#'     reports the total swing instead. See \code{h0testr::test()}.
#' @examples
#' set.seed(101)
#' ## no missing values: mnar_c0=-Inf, mnar_c1=0, mcar_p=0
#' samps <- h0testr::sim_samples(factors=list(condition=c("placebo", "drug")),
#'   n_per_cell=6)
#' sim <- h0testr::sim_design(samps, frm=~condition, test_term="condition", n_genes=25,
#'   n_genes_signif=5, effects=2, mnar_c0=-Inf, mnar_c1=0, mcar_p=0)
#' state <- sim$state
#' state$expression <- log2(state$expression + 1)
#' 
#' config <- sim$config   ## frm, test_term, the id columns and reference_levels, set
#'                        ##   to match what was simulated; save_state is FALSE
#' config$is_log_transformed <- TRUE   ## normalize() would; logged just above
#' rm(samps, sim)
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
#'   \code{logFC} is an effect size on the scale of \code{state$expression}: for a two level
#'     factor, the difference between its levels, so a log fold change when the input is log
#'     transformed; for a \strong{continuous} covariate, the change \strong{per unit} of it,
#'     whose size depends on the units the covariate is recorded in. This is what
#'     \code{h0testr::test()} reports as \code{logfc}. A joint test has no \code{logFC}
#'     column at all, having no single fold change to report, and \code{h0testr::test()}
#'     reports the total swing instead. See \code{h0testr::test()}.
#' @examples
#' set.seed(101)
#' ## no missing values: mnar_c0=-Inf, mnar_c1=0, mcar_p=0
#' samps <- h0testr::sim_samples(factors=list(condition=c("placebo", "drug")),
#'   n_per_cell=6)
#' sim <- h0testr::sim_design(samps, frm=~condition, test_term="condition", n_genes=25,
#'   n_genes_signif=5, effects=2, mnar_c0=-Inf, mnar_c1=0, mcar_p=0)
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

  if(method %in% c("msqrob", "msqrob_agg")) {

    dat <- SummarizedExperiment::rowData(result$fit[["genes"]])
    if(!all(c("msqrobModels", config$gene_id_col) %in% names(dat))) return(NULL)

    ## a mixed fit carries the estimates of the random effects alongside the fixed ones
    ##   and, with config$test_ridge, under renamed fixed names; see f.msqrob_parms():

    nom <- f.msqrob_parms(cols,
      method %in% "msqrob_agg" && isTRUE(config$test_ridge))

    models <- dat$msqrobModels
    coefs <- matrix(as.numeric(NA), nrow=length(models), ncol=length(cols),
      dimnames=list(as.character(dat[[config$gene_id_col]]), cols))

    for(idx in seq_along(models)) {
      beta <- try(msqrob2::getCoef(models[[idx]]), silent=T)
      if(inherits(beta, "try-error") || is.null(names(beta))) next
      coefs[idx, ] <- beta[match(nom, names(beta))]
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

  ## the mixed path fits one model per gene and falls back to stats::lm() for a gene
  ##   with a single feature, so its coefficients come from two collections of fits and
  ##   are assembled where those are, in f.prolfqua_mixed_f(), rather than re-walked
  ##   here:

  if(method %in% "prolfqua_lmer") {
    coefs <- result$coefs
    if(!is.matrix(coefs) || !all(cols %in% colnames(coefs))) return(NULL)
    return(coefs[, cols, drop=F])
  }

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

## helper for test(). Separate from f.format_limma() because DEqMS's table carries two
##   sets of statistics: the limma columns of the fit that was moderated, and the sca.*
##   columns DEqMS's own count-based prior produced. It was formatted by
##   f.format_limma(), whose first branch matches on the limma columns, all of which are
##   present; so the standardized table received limma's t and P.Value and DEqMS's
##   contribution reached only the original table. That made test_method "deqms" report
##   limma::eBayes(trend=FALSE) on the aggregated matrix in every column anything
##   downstream reads. The sca.* columns are the statistic the method name promises, so
##   they are what is reported here; limma's own numbers remain available under
##   test_method "trend" and in the original table.
##   Two shapes, as f.format_proda() has: a signed moderated t and a fold change when
##   one coefficient carries the test, a moderated F and no fold change when several do,
##   which f.fill_standard() then fills with the total swing. Which one it is comes from
##   the numerator degrees of freedom rather than from which columns are present, since
##   both shapes carry the same columns.
##   The lod column is left empty in both, as it is for every engine except the two that
##   report limma's own statistics. limma's B is a posterior log-odds computed from
##   limma's prior, and f.format_limma() reported it here; carrying it alongside a
##   p-value from DEqMS's prior would mix the two moderations, which is the thing being
##   fixed. It is still in the original table:

f.format_deqms <- function(tbl, config) {

  if(!is.data.frame(tbl)) {
    f.err("f.format_deqms: !is.data.frame(tbl); class(tbl): ",
      class(tbl), config=config)
  }

  nom <- c("gene", "AveExpr", "sca.t", "sca.F", "sca.df.num", "sca.P.Value",
    "sca.adj.pval")

  if(!all(nom %in% names(tbl))) {
    f.err("f.format_deqms: expected names not %in% names(tbl); missing:",
      paste(setdiff(nom, names(tbl)), collapse=", "), "; names(tbl):",
      names(tbl), config=config)
  }

  df_num <- unique(tbl$sca.df.num)

  if(length(df_num) != 1) {
    f.err("f.format_deqms: the table reports", length(df_num), "different numerator",
      "degrees of freedom (", paste(utils::head(df_num, 5), collapse=", "),
      "), but the hypothesis is a property of config$frm and config$test_term and so",
      "is the same for every gene", config=config)
  }

  if(df_num %in% 1) {

    if(!("logFC" %in% names(tbl))) {
      f.err("f.format_deqms: one coefficient carries the test but the table has no",
        "logFC column; names(tbl):", names(tbl), config=config)
    }

    tbl <- data.frame(feature=tbl$gene, expr=tbl$AveExpr,
      logfc=tbl$logFC, stat=tbl$sca.t, lod=as.numeric(NA),
      pval=tbl$sca.P.Value, adj_pval=tbl$sca.adj.pval)

  } else {

    tbl <- data.frame(feature=tbl$gene, expr=tbl$AveExpr,
      logfc=as.numeric(NA), stat=tbl$sca.F, lod=as.numeric(NA),
      pval=tbl$sca.P.Value, adj_pval=tbl$sca.adj.pval)
  }

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
    c("lm", "trend", "deqms", "msqrob", "msqrob_agg", "proda", "prolfqua",
      "prolfqua_lmer", "voom")
  )
}

## helper for test(): the settings only some engines read, so that one set on a run whose
##   method does not consult it can be said to have had no effect. Reported rather than
##   refused, and only when the value differs from the one new_config() ships: tune() sets
##   one config and varies the method, so most combinations of a sweep carry a setting
##   meant for one of the others and flagging those would be noise. Said from here rather
##   than from check_config(), which every step of the workflow calls and which would
##   repeat it once per step; this is the same place, and once per run, as the notes
##   test_prolfqua() and test_msqrob() make about the settings they do read. Keyed on the
##   method test() resolved rather than on config$test_method, test(method=) overriding
##   it. Returns the settings noted, so that a caller can check without reading the log:

f.note_ignored_settings <- function(method, config) {

  knobs <- list(
    test_random_obs=list(default=TRUE, methods=c("prolfqua_lmer", "msqrob_agg")),
    test_ridge=list(default=FALSE, methods="msqrob_agg")
  )

  ignored <- character(0)

  for(nom in names(knobs)) {
    if(!(nom %in% names(config)) || !(length(config[[nom]]) %in% 1)) next
    if(is.na(config[[nom]]) || config[[nom]] %in% knobs[[nom]]$default) next
    if(method %in% knobs[[nom]]$methods) next
    ignored <- c(ignored, paste0("config$", nom, "=", config[[nom]],
      " (read by test_method ", paste(knobs[[nom]]$methods, collapse=", "), ")"))
  }

  if(length(ignored)) {
    f.msg("NOTE: test: method", method, "does not consult",
      paste(ignored, collapse="; "), "\n",
      " so", if(length(ignored) > 1) "those settings have" else "that setting has",
      "no effect on this run", config=config)
  }

  return(ignored)
}

## The test methods whose variance prior can be fitted against mean feature intensity,
##   which is what limma::eBayes(trend=TRUE) does and what config$test_trend, or test()'s
##   trend argument, asks for. test_method "trend" is not here although it always trends:
##   these are the methods a TRUE can be handed to, and that one takes no such argument
##   because trending is its definition. The rest cannot: lm and msqrob do not moderate
##   against a covariate at all, proda fits its own prior, voom puts the mean-variance
##   relationship into precision weights so a trended prior on top of it would count the
##   same thing twice, and the mixed paths have a Satterthwaite denominator rather than
##   one residual variance to shrink. See f.note_trend() for what is said instead:

f.trend_methods <- function() {
  return(c("deqms", "prolfqua"))
}

## helper for test(): a request to trend that the resolved method cannot honor. A warning
##   rather than an error, per the same reasoning as f.note_ignored_settings(): tune() sets
##   one config and varies the method, so a sweep would otherwise die on the first method
##   that does not trend. Keyed on the resolved value rather than on config$test_trend, so
##   that an explicit test(trend=TRUE) is caught too, and naming the source, since which of
##   the two to change differs. test_method "trend" gets the opposite message: it always
##   trends, so TRUE is already what it does and only an explicit FALSE is a request it
##   cannot honor. Returns TRUE when something was said, so that a test need not read the
##   log:

f.note_trend <- function(method, trend, given, config) {

  src <- if(given) "the trend argument" else "config$test_trend"

  if(method %in% "trend") {
    if(given && !trend) {
      f.msg("WARNING: test: trend is FALSE but test_method 'trend' is",
        "limma::eBayes(trend=TRUE), so the prior is fitted against mean gene",
        "intensity regardless;", "\n",
        " for a flat prior use test_method 'deqms', whose limma prior this argument",
        "does set, or 'lm', which does not moderate at all", config=config)
      return(TRUE)
    }
    return(FALSE)
  }

  if(!trend || method %in% f.trend_methods()) return(FALSE)

  f.msg("WARNING: test:", src, "is TRUE, but test_method", method, "does not fit",
    "its variance prior against mean feature intensity, so the run is unaffected;",
    "\n", " the methods that do are", paste(f.trend_methods(), collapse=", "),
    "and 'trend', which always does", config=config)

  return(TRUE)
}

#' Hypothesis testing
#' @description
#'   Wrapper for various hypothesis testing methods.
#' @details
#'   Tests for differential expression using method specified in config. 
#'   See invididual \code{test_*} methods for more details. 
#'   The \code{method} setting meanings are: 
#'   \tabular{ll}{
#'     \code{lm}     \cr \tab Use \code{stats::lm()} on each feature. \cr
#'     \code{trend}  \cr \tab Use \code{limma::eBayes(trend=TRUE)}. \cr
#'     \code{deqms}  \cr \tab Use \code{DEqMS::spectraCounteBayes()}. \cr
#'     \code{msqrob} \cr \tab Use \code{msqrob2::msqrob()}. \cr
#'     \code{msqrob_agg} \cr \tab Use \code{msqrob2::msqrobAggregate()}: one mixed model per gene over the rows of its features, with the feature and the observation as random effects. \cr
#'     \code{proda}  \cr \tab Use \code{proDA::proDA()}. \cr
#'     \code{prolfqua} \cr \tab Use \code{prolfqua::strategy_lm()} on each feature. \cr
#'     \code{prolfqua_lmer} \cr \tab Use \code{prolfqua::strategy_lmer()}: one mixed model per gene over the rows of its features, with the feature and the observation as random effects. \cr
#'     \code{voom}   \cr \tab Use \code{limma::voom()}. \cr
#'   }
#'   Most methods are row-wise on \code{state$expression} and return one row per row of
#'     it, so their granularity follows the input. \code{"deqms"}, \code{"msqrob"},
#'     \code{"msqrob_agg"} and \code{"prolfqua_lmer"} return one row per gene whatever
#'     the input is: the first two aggregate internally, and the last two model the
#'     features of a gene instead of aggregating them, \code{"msqrob_agg"} aggregating
#'     only to carry the result. \code{h0testr::tune()} therefore does not call
#'     \code{combine_features()} before those four, and \code{h0testr::run()} should be
#'     given a \code{config$run_order} without it; \code{"prolfqua_lmer"} and
#'     \code{"msqrob_agg"} need feature level input and refuse a state whose features
#'     are already aggregated.
#'   The two feature level mixed model paths fit the same random structure by different
#'     engines. \code{"prolfqua_lmer"} takes Satterthwaite degrees of freedom for the
#'     tested contrast and is the better calibrated of the two; \code{"msqrob_agg"}
#'     reports \code{msqrob2}'s moderated t against \code{dfPosterior}, which is
#'     generous, and is mildly anti-conservative as a result. See
#'     \code{h0testr::test_msqrob()} for the simulated rejection rates.
#'   Several settings are read by some engines and not by others, so one set on a run
#'     whose method does not consult it does nothing at all:
#'     \code{config$test_random_obs} is read by \code{"prolfqua_lmer"} and
#'     \code{"msqrob_agg"}, and \code{config$test_ridge} by \code{"msqrob_agg"}. Such a
#'     setting is reported as a \code{NOTE} in the log rather than refused, and only when
#'     its value differs from the one \code{h0testr::new_config()} ships:
#'     \code{h0testr::tune()} sets one configuration and varies the method, so most
#'     combinations of a sweep carry a setting meant for one of the others. Said here,
#'     once per run, rather than in \code{h0testr::check_config()}, which every step of
#'     the workflow calls.
#'   The trended variance prior is the same kind of setting, but is resolved here for
#'     every method rather than read by each: the \code{trend} argument answers when
#'     given and \code{config$test_trend} otherwise, and the resolved value is passed to
#'     \code{"prolfqua"}, where it fits the prior against mean feature intensity, and to
#'     \code{"deqms"}, where it sets \code{limma::eBayes(trend=)} but does not change what
#'     is reported: \code{DEqMS} fits its own prior from quantities
#'     \code{limma::eBayes()} leaves alone, so only the \code{limma} columns of
#'     \code{original} move. \code{"trend"} always trends, that being its
#'     definition, so \code{TRUE} is silent there and only an explicit \code{FALSE} draws
#'     a remark. The remaining methods cannot trend: \code{"lm"} and \code{"msqrob"} do
#'     not moderate against a covariate, \code{"proda"} fits its own prior, \code{"voom"}
#'     carries the mean-variance relationship in its precision weights so a trended prior
#'     would count it twice, and the mixed paths have a Satterthwaite denominator rather
#'     than one residual variance to shrink. A \code{TRUE} reaching one of those is a
#'     \code{WARNING} in the log, not a refusal, for the same \code{h0testr::tune()}
#'     reason.
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
#'     \code{test_method}    \cr \tab Character scalar in \code{h0testr::test_methods()}. \cr
#'   }
#' @param method Name of test method where
#'   \code{method \%in\% h0testr::test_methods()}.
#' @param is_log_transformed Logical scalar: whether \code{state$expression} has
#'   been log transformed. Only consulted for \code{method \%in\% c("proda",
#'   "prolfqua", "prolfqua_lmer")}. Defaults to \code{config$is_log_transformed}, which
#'   \code{h0testr::initialize()} and \code{h0testr::normalize()} maintain;
#'   passing both is an error unless they agree.
#' @param prior_df Prior degrees of freedom for method \code{proda};
#'   where \code{2 <= prior_df <= n_features}.
#' @param trend Logical scalar: whether the variance prior of the moderation is fitted
#'   against mean feature intensity rather than being flat. Honored by \code{method
#'   \%in\% c("deqms", "prolfqua")}; \code{"trend"} always trends and the rest cannot,
#'   which is a \code{WARNING} rather than an error. Defaults to
#'   \code{config$test_trend}, and to \code{FALSE} when that is absent too; unlike
#'   \code{is_log_transformed}, an argument that disagrees with the configuration simply
#'   wins, this being a preference rather than a fact about the data. Unrelated to
#'   \code{method="trend"}, which names a different limma fit.
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
#'     of \code{state$expression} the result describes, except for the gene level
#'     methods, \code{method \%in\% c("deqms", "msqrob", "msqrob_agg",
#'     "prolfqua_lmer")}, which report one row per value of
#'     \code{config$gene_id_col} whatever level the input is at. So a
#'     method other than those reports precursors when handed precursors and
#'     genes when handed the output of \code{h0testr::combine_features()}, which
#'     sets \code{config$feat_col} to \code{config$gene_id_col}. The feature
#'     metadata reported in \code{original} is at the matching level: for the gene
#'     level methods it is the gene level form of \code{state$features},
#'     built exactly as \code{h0testr::combine_features()} builds it.
#'   \code{expr} is the mean of the values handed to the test, over the
#'     observations where the feature was seen; for the gene level methods, which
#'     take feature level input and report gene level results, it
#'     is also over the features of each gene.
#'   \code{stat} is the statistic the engine itself reports, with two exceptions, both
#'     for a joint test and both for the same reason, that the engine's API answers one
#'     coefficient at a time while its fit supports the joint test.
#'     \code{method \%in\% c("msqrob", "msqrob_agg")} reports an F computed by
#'     \code{h0testr} from the fitted \code{msqrob2} models, since
#'     \code{msqrob2::hypothesisTest()} answers one contrast at a time; see
#'     \code{h0testr::test_msqrob()}. \code{method="deqms"} reports an F computed from
#'     the per-gene variance \code{DEqMS::spectraCounteBayes()} fits, since
#'     \code{DEqMS} moderates one coefficient's t-statistic and has no F-analogue; see
#'     \code{h0testr::test_deqms()}. Both reduce to the engine's own statistic at one
#'     numerator degree of freedom.
#'   A \code{config$contrast} run has no such exception: every engine reports its
#'     own statistic for a contrast, the two \code{msqrob} methods with
#'     \code{msqrob2::hypothesisTest()} and \code{"deqms"} with the single coefficient
#'     \code{limma::contrasts.fit()} leaves it.
#' @examples
#' set.seed(101)
#' ## no missing values: mnar_c0=-Inf, mnar_c1=0, mcar_p=0
#' samps <- h0testr::sim_samples(factors=list(condition=c("placebo", "drug")),
#'   n_per_cell=6)
#' sim <- h0testr::sim_design(samps, frm=~condition, test_term="condition", n_genes=25,
#'   n_genes_signif=5, effects=2, mnar_c0=-Inf, mnar_c1=0, mcar_p=0)
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
#' 
#' out <- h0testr::test(out$state, out$config, method="trend")
#' head(out$original)
#' head(out$standard)
#' summary(out$fit)

test <- function(state, config, method=NULL,
    is_log_transformed=NULL, prior_df=NULL, trend=NULL) {

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

  if(method %in% c("proda", "prolfqua", "prolfqua_lmer")) {
    is_log_transformed <- f.is_log_transformed(is_log_transformed, config,
      "test")
  }
  
  ## resolved for every method rather than only the two that take it, so that the one
  ##   place the value is settled is also the place a method that cannot honor it says
  ##   so; the argument overrides config$test_trend. See f.is_trend():

  trend_given <- !(is.null(trend) || (is.character(trend) && all(trend %in% "")))
  trend <- f.is_trend(trend, config, "test")

  f.msg("test: method:", method, "; is_log_transformed:", is_log_transformed,
    "; prior_df:", prior_df, "; trend:", trend, config=config)

  ## a setting the resolved method does not read, said once here rather than once per
  ##   step by check_config(); see f.note_ignored_settings() and f.note_trend():

  f.note_ignored_settings(method, config)
  f.note_trend(method, trend, trend_given, config)

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
    result <- test_deqms(state, config, trend=trend)
    tbl2 <- f.format_deqms(result$hits, config)
  } else if(method %in% "msqrob") {
    result <- test_msqrob(state, config)
    tbl2 <- f.format_msqrob(result$hits, test_col, config)
  } else if(method %in% "msqrob_agg") {
    result <- test_msqrob(state, config, aggregate=TRUE)
    tbl2 <- f.format_msqrob(result$hits, test_col, config)
  } else if(method %in% "proda") {
    result <- test_proda(state, config,
      is_log_transformed=is_log_transformed, prior_df=prior_df)
    tbl2 <- f.format_proda(result$hits, config)
  } else if(method %in% "prolfqua") {
    result <- test_prolfqua(state, config,
      is_log_transformed=is_log_transformed, trend=trend)
    tbl2 <- f.format_prolfqua(result$hits, test_col, config)
  } else if(method %in% "prolfqua_lmer") {
    result <- test_prolfqua(state, config,
      is_log_transformed=is_log_transformed, mixed=TRUE)
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
