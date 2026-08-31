## Apply config$contrast to a limma fit:
##   limma::contrasts.fit() re-expresses the fit in terms of contrasts, 
##   leaving one coefficient per contrast:

f.limma_contrast_fit <- function(fit, design, config) {

  if(is.null(design$contrast)) return(list(fit=fit, coef=design$cols_test))

  con <- matrix(design$contrast, ncol=1,
    dimnames=list(colnames(design$X), trimws(config$contrast)))

  return(list(fit=limma::contrasts.fit(fit, contrasts=con), coef=1L))
}

## F statistic and p-value for a nested pair of linear models, given 
##   residual sums of squares and residual degrees of freedom of each. Exact for a
##   linear model, and equal to the Wald F-test of the same hypothesis:

f.nested_f <- function(rss_red, rss_full, df_red, df_full, s2_err, df_err) {

  df_num <- df_red - df_full
  ss <- rss_red - rss_full
  ms <- ss / df_num
  fval <- ms / s2_err

  return(data.frame(Df=df_num, Sum.Sq=ss, Mean.Sq=ms, F.value=fval,
    p.value=stats::pf(fval, df_num, df_err, lower.tail=F)))
}

## Rows of L for f.wald_f(): hypothesis as a matrix over the coefficients of
##   the design, which is the identity rows of the tested columns for config$test_term
##   and the single row of weights for config$contrast:

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

## The Wald F-test of L %*% beta == 0 from one fitted model. For a mixed model this is
##   an F whose denominator degrees of freedom come from a Satterthwaite approximation:
##   a covariate is tested against between-observation variation rather than
##   the feature by observation residual; for a least squares fit it is the ordinary
##   F, which is also the nested model comparison f.nested_f() performs. Both are
##   needed because the mixed path falls back to stats::lm() for a gene with a single
##   feature.
##   Differencing the residual sums of squares of two mixed fits, which is what
##   f.prolfqua_nested_f() does for the least squares path, is not available: 
##   variance components re-estimated for each fit, so sigma^2 times the residual
##   degrees of freedom is not a decomposable residual sum of squares, and the
##   difference is not an F. Hence a Wald test of one fit rather comparing two.
##   NULL when the fit cannot support the test:

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

  ## Wald statistic is extra sum of squares divided by error variance:

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

## Error variance and its degrees of freedom, shrunk across features toward a
##   common prior. 
##   Moderation does not depend on how many df are being tested. A
##   moderated F is the ordinary F with the per-feature residual variance replaced by
##   the posterior variance and the denominator df raised by those of
##   the prior, which at one numerator df is the square of the moderated t. 
##   A covariate makes prior a function of that covariate instead of one number
##   shared by every feature, fitted as a natural spline of up to four degrees of
##   freedom on the log residual variances. Similar to mean-variance trend of 
##   limma::eBayes(trend=TRUE):

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

  ## a spline needs more features than it has degrees of freedom;
  ##   a covariate without spread is a poor guide:

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

  sv <- limma::squeezeVar(s2, df=df_resid, covariate=covariate, robust=FALSE)

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

## prior variance as one line of a log: a number when flat, a range when trend:

f.prior_txt <- function(x) {

  if(length(unique(x)) > 1) {
    return(paste("range", signif(min(x), 4), "to", signif(max(x), 4)))
  }

  return(as.character(signif(x[1], 4)))
}

