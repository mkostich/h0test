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

  ## before combine_features(), so that the design f.design_test_cols() builds below
  ##   from the aggregated state carries the reference level config declares rather
  ##   than the alphabetically first one; see test_lm() for what this prevents. A
  ##   no-op for a run that came through test():

  state <- f.relevel_state_covariates(state, config, caller="test_deqms")

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

