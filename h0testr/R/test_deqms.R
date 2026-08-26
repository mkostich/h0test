## Number of features behind each gene:

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

## Variance prior DEqMS::spectraCounteBayes() fits, fitted on genes that can
##   carry one. Its own row handling is silently wrong otherwise: 

f.deqms_prior <- function(fit, config, who="f.deqms_prior") {

  nom <- rownames(fit$coefficients)
  n <- length(nom)

  ## log2(0) is -Inf and would be dropped by na.omit the same way:
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

  ## whole point of this function:
  if(!identical(rownames(sub$coefficients), nom[ok])) {
    f.err(who, ": subsetting the fit did not preserve the gene order, so the prior",
      "cannot be mapped back;", "\n", "  first few expected:",
      utils::head(nom[ok], 5), "; first few found:",
      utils::head(rownames(sub$coefficients), 5), config=config)
  }

  ## loess needs those counts spread rather than merely varied:
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

  ## per-gene quantity of subset fit, back in row order of whole fit:

  f.fill <- function(x) {
    out <- rep(as.numeric(NA), n)
    names(out) <- nom
    out[ok] <- as.numeric(x)
    return(out)
  }

  fit$sca.postvar <- f.fill(sub$sca.postvar)
  fit$sca.priorvar <- f.fill(sub$sca.priorvar)

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

## Moderated test of design matrix columns carrying test, from a fit
##   DEqMS::spectraCounteBayes() has moderated:

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

    ## DEqMS's own statistic, recomputed: coefficient over stdev.unscaled times 
    ##   square root of posterior variance:

    se <- fit$stdev.unscaled[, cols] * sqrt(post_var)
    tval <- betas[, 1] / se
    fval <- tval^2
    df_num <- 1L

  } else {

    ## unscaled covariance of tested coefficients, indexed by name:
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

    ## rowSums() rather than loop over genes: same quadratic form
    ##   t(b) %*% Vinv %*% b for every row of betas, and matrix is small:

    quad <- rowSums((betas %*% Vinv) * betas)
    fval <- (quad / df_num) / post_var
  }

  pval <- stats::pf(fval, df_num, df_den, lower.tail=F)

  ## gene whose posterior variance came back non-finite, which is a gene f.deqms_prior()
  ##   left out of prior because it has no residual df:

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
#'   A test of one coefficient is reported as \code{DEqMS}'s own moderated t, with a
#'     \code{logFC} column; a test of several is reported as a moderated F, and no
#'     \code{logFC}, since several coefficients have no single fold change.
#'     \code{DEqMS} itself reports only the former: \code{DEqMS::spectraCounteBayes()}
#'     moderates the t-statistic of one coefficient and the package has no F-analogue
#'     anywhere. The joint test is the ordinary F with the residual variance
#'     replaced by the posterior one and the denominator degrees of freedom raised by
#'     those of the prior. At one numerator degree of freedom that F is the square of
#'     \code{DEqMS}'s own moderated t.
#'   \code{config$contrast} tests a weighted sum of coefficients instead of a term, is
#'     one degree of freedom however many coefficients it weights, and reaches this
#'     engine as a single coefficient. 
#'   A gene with no residual degrees of freedom, which is a gene the aggregated matrix
#'     has as many observations of as \code{config$frm} has design columns, is left out
#'     of that prior and reported with \code{NA} in every \code{sca.} column.
#'   Returns gene-level hypothesis testing results based on peptide/precursor-level input.
#'   Flow is:
#'     \tabular{l}{
#'       1. for each gene, count number of associated peptides. \cr
#'       2. Fit linear model to \code{config$frm} using \code{limma::lmFit()}. \cr
#'       3. Calculate statistics using \code{limma::eBayes()} on fitted model. \cr
#'       4. Append peptide counts to model returned by \code{limma::eBayes()}. \cr
#'       5. Adjust statistics using \code{DEqMS::spectraCounteBayes()}, fitted on the
#'            genes that have residual degrees of freedom. \cr
#'       6. Form moderated test of the coefficients carrying the test from the
#'            variance prior that fitted. \cr
#'       7. Generate hit table with \code{limma::topTable()}, with the moderated
#'            statistics appended. \cr
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
#'     \code{test_trend}    \cr \tab Optional logical; answers the \code{trend} argument when that is not given. Defaults to \code{FALSE} when both are absent. Unrelated to \code{config$test_method="trend"}. \cr
#'   }
#' @param trend Logical scalar. Whether \code{limma::eBayes()} fits its variance prior
#'   against mean gene intensity rather than shrinking every gene toward one number.
#'   Defaults to \code{config$test_trend}, and to \code{FALSE} when that is absent. 
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
#' ##   which reports DEqMS's own moderated t:
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
#' ## same fit with limma's prior fitted against mean gene intensity instead of
#' ##   flat. Equivalently config$test_trend <- TRUE, which h0testr::test(method="deqms")
#' ##   passes through. limma's P.Value moves and DEqMS's sca.P.Value does not, the
#' ##   count-based prior the latter comes from being fitted from quantities
#' ##   limma::eBayes() leaves alone; see the note on the trend argument:
#' trended <- h0testr::test_deqms(out$state, out$config, trend=TRUE)
#' i <- rownames(result$hits)
#' c(limma=cor(result$hits$P.Value, trended$hits[i, "P.Value"]),
#'   deqms=cor(result$hits$sca.P.Value, trended$hits[i, "sca.P.Value"]))

test_deqms <- function(state, config, trend=NULL) {

  ## NULL rather than FALSE, so that config$test_trend reaches this engine:

  trend <- f.is_trend(trend, config, "test_deqms")

  check_config(config)
  f.check_state(state, config)

  state <- f.relevel_state_covariates(state, config, caller="test_deqms")
  save_state <- config$save_state
  config$save_state <- FALSE
  out <- combine_features(state, config, method="medianPolish", rescale=FALSE)
  config$save_state <- save_state
  design <- f.design_test_cols(out$state, out$config)
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

  fit <- limma::lmFit(out$state$expression, design$X)

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

  f.msg("test_deqms: limma::eBayes prior fitted against", if(trend) {
    "mean gene intensity (trend=TRUE)"
  } else "one number for every gene (trend=FALSE)", "\n",
    " this sets limma's own columns of the result; the reported sca.* statistics come",
    "from the prior DEqMS fits against the feature counts, which is unaffected",
    config=config)

  fit <- limma::eBayes(lc$fit, trend=trend)
  fit$count <- counts[rownames(fit$coefficients)]

  if(any(is.na(fit$count))) {
    f.err("test_deqms: no feature count for", sum(is.na(fit$count)), "of",
      length(fit$count), "genes of the aggregated expression matrix;", "\n",
      "  first few:", utils::head(rownames(fit$coefficients)[is.na(fit$count)], 5),
      config=config)
  }

  fit <- f.deqms_prior(fit, out$config, "test_deqms")
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

