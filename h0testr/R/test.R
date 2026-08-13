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
##   features by the rank of these same two matrices:

f.test_lm_feat <- function(y, X, cols_test, cols_report) {

  i <- !is.na(y)
  yy <- y[i]
  xf <- X[i, , drop=F]
  xr <- X[i, -cols_test, drop=F]

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

  f.msg("test_lm: test_term:", config$test_term, "; fdr.method:", fdr.method,
    "; design columns:", ncol(X), "; test columns:", length(cols_test),
    "; df:", design$df_intend, config=config)

  hits <- t(apply(state$expression, 1, f.test_lm_feat, X, cols_test,
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
#'   If the number of features per gene/protein-group is the same for all
#'     features, returns same result as \code{h0testr::test_trend()}.
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
  
  tbl <- table(state$features[, config$gene_id_col], useNA="ifany")
  tbl <- as.data.frame(tbl)
  counts <- tbl$Freq
  names(counts) <- tbl$Var1
  
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
  idx <- design$cols_test

  fit <- limma::lmFit(out$state$expression, design$X)
  fit <- limma::eBayes(fit, trend=trend)
  fit$count <- counts[rownames(fit$coefficients)]
  
  if(length(unique(fit$count)) >= 2) {
    fit  <- DEqMS::spectraCounteBayes(fit, fit.method="loess")
    hits <- DEqMS::outputResult(fit, coef_col=idx)
  } else {
    f.msg("WARNING: length(unique(fit$count)) < 2;",
      "falling back to h0testr::test_trend(); unique(fit$count):",
      unique(fit$count),
      config=config
    )
    hits <- test_trend(state, config)$hits
  }

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
#'       5. Calculate test statistics with \code{msqrob2::hypothesisTest()}. \cr
#'     }
#'   \code{msqrob2::hypothesisTest()} returns one table per contrast rather than a
#'     joint test over several, so \code{config$test_term} must resolve to a single
#'     design matrix column; it is an error if it does not. That rules out testing a
#'     factor with more than two levels, and also testing a variable that appears in
#'     an interaction, since by marginality the test then covers every term containing
#'     the variable. Use \code{h0testr::test_lm()}, \code{h0testr::test_trend()} or
#'     \code{h0testr::test_voom()} for those tests.
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
#'     \code{reference_levels} \cr \tab Named character vector with the reference level of each factor variable in \code{config$frm}. \cr
#'     \code{covariate_types}  \cr \tab Optional; classification of variables in \code{config$frm}, as set by \code{initialize()}. \cr
#'     \code{factor_levels}    \cr \tab Optional; resolved levels of each factor variable, as set by \code{initialize()}. \cr
#'   }
#' @param maxit Integer scalar >= 1. How many iterations to use for \code{rlm} fitting.
#' @return
#'   A list with components:
#'   \tabular{ll}{
#'     \code{hits}  \cr \tab \code{data.frame} of results; columns \code{config$gene_id_col} and: 
#'       \code{c("nNonZero .n", "logFC", "se", "df", "t", "pval", "adjPval")}. \cr
#'     \code{fit}   \cr \tab Model returned by \code{msqrob2::hypothesisTest()}. \cr
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
  feats <- state$features
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

  ## the design and the single column carrying the test. msqrob2::hypothesisTest()
  ##   returns one table per contrast rather than a joint test over several, so only a
  ##   test that resolves to one coefficient can be run; f.design_test_cols_max()
  ##   errors otherwise. See test_deqms() for what selecting by coefficient name used
  ##   to hide:

  design <- f.design_test_cols_max(state, config, "test_msqrob", max_cols=1)
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
      "config$test_term", config$test_term, "is not among the parameters",
      "msqrob2::msqrob() fit;", "\n", "parameters fit:",
      paste(parms, collapse=", "), config=config)
  }

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
  
  return(list(hits=tbl, fit=obj))
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

  if(length(design$cols_test) %in% 1) {

    ## one design column carries the test, so the Wald test on that coefficient is
    ##   the test config$test_term names. See test_deqms() for what selecting by
    ##   coefficient name used to hide. proDA::result_names() quotes coefficient
    ##   names containing ':':

    col_pick <- cols_pick
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

    x_red <- design$X[, -design$cols_test, drop=F]

    ## proDA::test_diff() requires a full rank reduced model, and reports a rank
    ##   deficient one as colinear covariates, which says nothing about which part
    ##   of config$frm is responsible:

    rank_red <- f.design_rank(x_red)

    if(rank_red < ncol(x_red)) {
      f.err("test_proda: the reduced model formed by dropping config$test_term '",
        config$test_term, "' is rank deficient, so the likelihood ratio test",
        "against it is not defined;", "\n", "config$frm:",
        deparse(design$parsed$frm), "; reduced model columns:", colnames(x_red),
        "; rank:", rank_red, config=config)
    }

    f.msg("test_proda: config$test_term", config$test_term, "carries",
      length(design$cols_test), "design columns (",
      paste(cols_pick, collapse=", "), "), so testing by likelihood ratio against",
      "the reduced model rather than by a single contrast", config=config)

    tbl <- proDA::test_diff(fit, reduced_model=x_red, sort_by="pval")
  }

  return(list(hits=tbl, fit=fit))
}

#' Hypothesis testing using the \code{prolfq} package
#' @description
#'   Tests for differential expression using the \code{prolfq::build_model()} function.
#' @details
#'   Uses the \code{prolfqua::build_model()} function. Returned results sorted by
#'     p-value.
#'   Covariates classified as \code{"factor"} in \code{config$covariate_types} are
#'     registered with \code{prolfqua::AnalysisTableAnnotation$factors}, which holds
#'     categorical annotations only; continuous (numeric) covariates are passed
#'     through unregistered. Both are supported: \code{prolfqua::build_model()} fits
#'     a plain \code{stats::lm()} to the data.frame as given, so a numeric covariate
#'     is read as continuous and contributes a single one degree of freedom row to
#'     the reported ANOVA table.
#'   The level ordering set by \code{initialize()} is preserved, so the
#'     coefficients of the returned \code{fit} are relative to the reference
#'     level declared in \code{config$reference_levels}. The returned
#'     \code{hits} are term-level ANOVA F-tests, which do not depend on the
#'     choice of reference level.
#'   Results are one row per term of the model rather than one per coefficient, so
#'     a factor with more than two levels is a single multi-df F-test over its
#'     contrasts. What that table cannot express is a joint test over several
#'     terms, which is what the marginality rule makes \code{config$test_term}
#'     mean when the named variable also appears in an interaction: testing
#'     \code{"grp"} in \code{~grp*sex} covers both \code{grp} and \code{grp:sex}.
#'     That is an error here rather than a silently narrower test; use
#'     \code{test_method} \code{"lm"}, \code{"trend"} or \code{"voom"} for such a
#'     \code{test_term}. \code{h0testr::tune()} skips those combinations instead of
#'     stopping.
#'   The terms of \code{config$frm} are reordered so that the tested term comes
#'     last. \code{stats::anova()} decomposes a model sequentially (Type I sums of
#'     squares), so a term is adjusted only for the terms entered before it, and
#'     the p-value for \code{config$test_term} would otherwise depend on where in
#'     \code{config$frm} it was written. Reordering makes the reported test the one
#'     adjusted for every other term; the fit itself is unchanged.
#'   The false discovery rate \code{prolfqua} reports is computed within each term,
#'     so the \code{FDR} column of \code{hits} is adjusted over the features tested
#'     and not over the other terms of the model.
#'   Flow is:
#'     \tabular{l}{
#'       1. Reshape data into long format with intensities, feature meta, and sample meta. \cr
#'       2. Make and populate \code{prolfqua::AnalysisTableAnnotation} object. \cr
#'       3. Make \code{prolfqua::LFQData} object from data and
#'            \code{prolfqua::AnalysisTableAnnotation} object. \cr
#'       4. Make \code{prolfqua::strategy_lm} object from \code{config$frm}, with the
#'            tested term moved last. \cr
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
#'     \code{frm}                   \cr \tab Formula (formula) to be fit. \cr
#'     \code{test_term}             \cr \tab Term (character) to be tested for non-zero coefficient. \cr
#'     \code{reference_levels}      \cr \tab Named character vector with the reference level of each factor variable in \code{config$frm}. \cr
#'     \code{covariate_types}       \cr \tab Optional; classification of variables in \code{config$frm}, as set by \code{initialize()}. \cr
#'     \code{factor_levels}         \cr \tab Optional; resolved levels of each factor variable, as set by \code{initialize()}. \cr
#'     \code{normalization_method}  \cr \tab If present and \code{is_log_transformed} unset, used to infer it. \cr
#'   }
#' @param is_log_transformed Logical scalar: whether \code{state$expression} has
#'   been log transformed. Defaults to \code{config$is_log_transformed}, which
#'   \code{h0testr::initialize()} and \code{h0testr::normalize()} maintain;
#'   passing both is an error unless they agree.
#' @return
#'   A list with components:
#'   \tabular{ll}{
#'     \code{hits}  \cr \tab \code{data.frame} of the \code{prolfqua::build_model()}
#'       ANOVA rows for \code{config$test_term}: the feature id columns
#'       (\code{config$gene_id_col} and \code{config$feat_id_col}, which are the same
#'       column once \code{combine_features()} has run) and:
#'       \code{c("isSingular", "nrcoef", "factor", "Df", "Sum.Sq", "Mean.Sq",
#'       "F.value", "p.value", "FDR")}. One row per feature. \cr
#'     \code{fit}   \cr \tab Model returned by \code{prolfqua::build_model()}. \cr
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
#'   ## no grp:sex term here: this function reports one term of the model at a
#'   ##   time, so with ~grp+sex+grp:sex testing 'grp' would be a joint test over
#'   ##   grp and grp:sex, which it refuses; use 'lm', 'trend' or 'voom' for that:
#'   frm=~grp+sex,
#'   test_term="grp",
#'   reference_levels=c(grp="ctl", sex="F")
#' )
#' out <- h0testr::initialize(state, config, minimal=TRUE)
#'
#' ## actual test:
#' result <- h0testr::test_prolfqua(out$state, out$config, is_log_transformed=FALSE)
#' head(result$hits)

test_prolfqua <- function(state, config, is_log_transformed=NULL) {

  is_log_transformed <- f.is_log_transformed(is_log_transformed, config,
    "test_prolfqua")
  
  ## prolfqua reports the rows of an anova table, one per term of the model, so a
  ##   factor with several levels is already a correct multi-df F-test here; what
  ##   the table cannot express is a joint test over several terms.
  ##   f.normalize_terms() gives exactly the terms the test of config$test_term
  ##   covers, by the marginality rule: testing a variable tests every term
  ##   containing it, so testing 'grp' in ~grp*sex covers both grp and grp:sex.
  ##   Checked here, before the reshaping below, so that an impossible test fails
  ##   before any work is done, and refused rather than quietly reported as the
  ##   test of one of the terms:

  parsed <- f.parse_frm(config$frm, config)
  norm <- f.normalize_terms(config)
  drops <- norm$drop_terms

  if(length(drops) > f.test_max_terms("prolfqua")) {
    f.err("test_prolfqua: testing config$test_term '", config$test_term, "' in",
      deparse(parsed$frm), "is a joint test over", length(drops), "terms (",
      paste(drops, collapse=", "), "), but test_prolfqua() reports one term at a",
      "time;", "\n", "use test_method 'lm', 'trend' or 'voom' for this test_term,",
      "or name a term of config$frm that no higher-order term contains",
      config=config)
  }

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
  
  trms <- sort(unique(parsed$vars))

  meta <- prolfqua::AnalysisTableAnnotation$new()
  meta$workIntensity <- "intensity"
  meta$is_response_transformed <- is_log_transformed
  meta$hierarchy[[config$gene_id_col]] <- config$gene_id_col
  meta$hierarchy[[config$feat_id_col]] <- config$feat_id_col

  ## prolfqua::AnalysisTableAnnotation$factors holds categorical annotations only,
  ##   so only factor covariates are registered there; declaration in
  ##   config$reference_levels is not the test, since logical covariates and
  ##   covariates already stored as factors are categorical without being declared.
  ##   A continuous covariate is passed through unregistered rather than refused:
  ##   the coercion that would break it, as.character() over $factors inside
  ##   prolfqua::setup_analysis(), never runs here, since
  ##   prolfqua::LFQData$new() takes setup=FALSE by default and so keeps dat
  ##   verbatim. Nothing downstream reads $factors either: build_model() is handed
  ##   obj$data with subject_Id from $hierarchy, and strategy_lm() fits a plain
  ##   stats::lm(), which reads a numeric column as continuous:

  types <- f.covariate_types(state, config)

  for(trm in trms) {

    if(types[[trm]] %in% "factor") {

      meta$factors[[trm]] <- trm

      ## initialize() has already ordered the levels, with the declared reference
      ##   level first; coercing to character here would leave the downstream fit
      ##   to re-derive the reference level by sorting, silently changing the
      ##   meaning of the reported coefficients:

      dat[, trm] <- f.relevel_covariate(samps[dat$sample, trm, drop=T], trm, config,
        "test_prolfqua")

    } else {

      ## a continuous term is one term of one degree of freedom in the anova table
      ##   prolfqua reports, so it is within the one term at a time bound
      ##   f.test_max_terms("prolfqua") sets, whether it is the tested term or an
      ##   adjustment:

      dat[, trm] <- samps[dat$sample, trm, drop=T]
    }
  }
  obj <- prolfqua::LFQData$new(data=dat, config=meta)
  
  ## prolfqua::build_model() gets its p-values from stats::anova(), which
  ##   decomposes the model sequentially (Type I sums of squares), so the row for a
  ##   term is adjusted only for the terms entered before it: with ~grp+sex the grp
  ##   p-value ignores sex entirely, and ~sex+grp and ~grp+sex give different
  ##   p-values for the same test on the same data. Putting the tested term last
  ##   makes its row the test adjusted for every other term, which is the
  ##   hypothesis config$test_term names. Only the order of the decomposition
  ##   changes; the fit, the span of the design and the residuals are the same.
  ##   Safe because the check above leaves exactly one term to move, and
  ##   f.test_term_drops() has already refused a term that a higher-order term
  ##   contains:

  lbls <- c(setdiff(parsed$labels, drops), drops)
  rhs <- paste(c(if(parsed$intercept %in% 0) "0", lbls), collapse=" + ")
  frm <- paste("intensity ~", rhs)
  strgy <- prolfqua::strategy_lm(frm)

  f.msg("test_prolfqua: test_term:", config$test_term, "; model:", frm,
    config=config)

  model <- prolfqua::build_model(data=obj$data, model_strategy=strgy,
    subject_Id=obj$config$hierarchy_keys())

  tbl <- as.data.frame(model$get_anova())

  ## match on canonicalized term labels: f.canon_label() has sorted the variables
  ##   of config$test_term alphabetically, but prolfqua's 'factor' column keeps the
  ##   order they appear in config$frm, so ~sex*grp with test_term 'grp:sex'
  ##   matched nothing and returned an empty table without complaint:

  i <- f.canon_label(tbl$factor) %in% norm$test_term

  if(!any(i)) {
    f.err("test_prolfqua: no anova rows for config$test_term '", config$test_term,
      "'; terms in the anova table:", paste(unique(tbl$factor), collapse=", "),
      config=config)
  }

  f.msg("tested", length(unique(tbl[[config$feat_col]])), "features; found",
    sum(tbl$FDR[i] < 0.05, na.rm=T), "hits", config=config)

  return(list(hits=tbl[i, , drop=F], fit=model))
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
  fit <- limma::eBayes(fit, trend=F)

  ## a single coefficient gives a t-test and a logFC column; several give an F-test:

  tbl <- limma::topTable(fit, coef=design$cols_test, number=Inf)

  f.msg("test_voom: test_term:", config$test_term, "; design columns:",
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
  fit <- limma::eBayes(fit, trend=T)

  ## a single coefficient gives a t-test and a logFC column; several give an F-test:

  tbl <- limma::topTable(fit, coef=design$cols_test, number=Inf)

  f.msg("test_trend: test_term:", config$test_term, "; design columns:",
    ncol(design$X), "; test columns:", length(design$cols_test), "; df:",
    design$df_intend, config=config)
  f.msg("tested", nrow(state$expression), "features", config=config)
  f.msg("found", sum(tbl$adj.P.Val < 0.05, na.rm=T), "hits", config=config)
  
  return(list(hits=tbl, fit=fit))
}

## helper for test():

f.format_lm <- function(tbl, config) {
  
  if(!is.data.frame(tbl)) {
    f.err("f.format_lm: !is.data.frame(tbl); class(tbl): ", 
      class(tbl), config=config)
  }
  
  if(!(config$gene_id_col %in% names(tbl))) {
    f.err("f.format_lm: expected names not %in% names(tbl); names(tbl):", 
      names(tbl), config=config)
  }
  
  nom <- c("pval", "p.adj")
  if(!all(nom %in% names(tbl))) {
    f.err("f.format_lm: expected names not %in% names(tbl); names(tbl):", 
      names(tbl), "; expected names:", nom, config=config)
  }
  
  tbl <- data.frame(feature=tbl[[config$gene_id_col]], expr=as.numeric(NA), 
    logfc=as.numeric(NA), stat=as.numeric(NA), lod=as.numeric(NA), 
    pval=tbl$pval, adj_pval=tbl$p.adj)
  
  tbl <- tbl[order(tbl$pval, decreasing=F), , drop=F]
  rownames(tbl) <- NULL
  
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

  f.design_test_cols(state, config)

  if(method %in% "lm") {
    result <- test_lm(state, config)
    tbl2 <- f.format_lm(result$hits, config)
  } else if(method %in% "trend") {
    result <- test_trend(state, config)
    tbl2 <- f.format_limma(result$hits, config)
  } else if(method %in% "deqms") {
    result <- test_deqms(state, config)
    tbl2 <- f.format_limma(result$hits, config)
  } else if(method %in% "msqrob") {
    result <- test_msqrob(state, config)
    tbl2 <- f.format_msqrob(result$hits, config)
  } else if(method %in% "proda") {
    result <- test_proda(state, config, 
      is_log_transformed=is_log_transformed, prior_df=prior_df)
    tbl2 <- f.format_proda(result$hits, config)
  } else if(method %in% "prolfqua") {
    result <- test_prolfqua(state, config, 
      is_log_transformed=is_log_transformed)
    tbl2 <- f.format_prolfqua(result$hits, config)
  } else if(method %in% "voom") {
    result <- test_voom(state, config)
    tbl2 <- f.format_limma(result$hits, config)
  } else if(method %in% "none") {
    f.msg("skipping testing: method %in% 'none'", config=config)
    return(NULL)
  } else f.err("test: unexpected method:", method, config=config)
  
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
