## Per-feature F-tests comparing prolfqua models of a full and a reduced design;
##   error variance moderated across features unless config$test_moderate is
##   FALSE, against a flat prior unless trend is TRUE and a covariate to
##   fit the prior against is supplied (specified by config$feat_id_col):

f.prolfqua_nested_f <- function(fit_full, fit_red, design, config, covariate=NULL,
    trend=NULL) {

  idvars <- unique(c(config$gene_id_col, config$feat_id_col))

  a <- as.data.frame(fit_full$model_df)
  a <- a[, c(idvars, "isSingular", "nr_coef", "df.residual", "sigma"), drop=F]
  names(a)[names(a) %in% "nr_coef"] <- "nrcoef"
  b <- as.data.frame(fit_red$model_df)
  b <- b[, c(idvars, "df.residual", "sigma"), drop=F]
  names(b)[names(b) %in% "df.residual"] <- "df_red"
  names(b)[names(b) %in% "sigma"] <- "sigma_red"

  tbl <- merge(a, b, by=idvars, all=F, sort=F)
  n_in <- nrow(tbl)

  ## sigma is NaN with no residual df left:

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

  ## prolfqua reports stats::sigma():
  rss_red <- tbl$sigma_red^2 * tbl$df_red
  rss_full <- tbl$sigma^2 * tbl$df.residual

  ord <- f.nested_f(rss_red=rss_red, rss_full=rss_full, df_red=tbl$df_red,
    df_full=tbl$df.residual, s2_err=tbl$sigma^2, df_err=tbl$df.residual)

  ## moderated statistic reported:
  moderate <- is.null(config$test_moderate) || isTRUE(config$test_moderate)

  trend <- f.is_trend(trend, config, "f.prolfqua_nested_f")

  if(moderate) {

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

## fixed effects of a fit from either strategy, under the names test_prolfqua()
##   give design matrix columns. stats::coef() of a mixed model returns per
##   group coefficients as list (not fixed effects):

f.fixed_coefs <- function(fit) {
  if(inherits(fit, "merMod")) return(lme4::fixef(fit))
  return(stats::coef(fit))
}

## fit collection of models:

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

## mixed path of test_prolfqua(): one model per gene over the long table of its
##   features, giving gene level inference from feature level data without aggregating.
##   Random effects are a random intercept per feature and a random intercept per observation.
##   Feature effect absorbs feature baselines, which differ by orders of magnitude 
##   among the precursors of one protein. Random rather than fixed effect
##   per feature: costs one variance component instead of one coefficient per
##   feature; shrinks toward the gene mean:

f.prolfqua_mixed <- function(obj, design, cols, config, trend=NULL) {

  dat <- obj$data_long()

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

  ## config$test_moderate and config$test_trend describe shrinkage of per feature
  ##   error variance across features:

  ignored <- character(0)
  if("test_moderate" %in% names(config) && isTRUE(config$test_moderate)) {
    ignored <- c(ignored, "config$test_moderate")
  }
  ## named by what is actually set rather than by which the caller used: 

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

  ## genes with single observed feature are fitted by least squares:
  
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
        subject_id=config$gene_id_col), "lmer", config))
  }

  if(length(single)) {
    f.msg("test_prolfqua: fitting", length(single), "gene(s) with a single observed",
      "feature by least squares instead:", frm_fix, config=config)
    mods[[length(mods) + 1]] <- list(fit_type="lm",
      model=f.quiet_fits(prolfqua::build_model(
        data=dat[genes %in% single, , drop=F],
        model_strategy=prolfqua::strategy_lm(frm_fix),
        subject_id=config$gene_id_col), "lm", config))
  }

  out <- f.prolfqua_mixed_f(mods, design, config, "test_prolfqua")

  f.msg("tested", nrow(out$hits), "genes; found",
    sum(out$hits$FDR < 0.05, na.rm=T), "hits", config=config)

  ## two collections of fits returned separately:

  types <- vapply(mods, function(m) m$fit_type, character(1))

  return(list(hits=out$hits, coefs=out$coefs,
    fit=if("lmer" %in% types) mods[[which(types %in% "lmer")[1]]]$model else NULL,
    fit_lm=if("lm" %in% types) mods[[which(types %in% "lm")[1]]]$model else NULL,
    fit_reduced=NULL, design=design))
}

## per-gene Wald F-tests from mixed path of test_prolfqua(), where one model is
##   fitted per gene over its features rather than one model per feature:

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

    mdf <- as.data.frame(mod$model$model_df)
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
        isSingular=mdf$isSingular[idx], nrcoef=mdf$nr_coef[idx],
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
#'   Tests for differential expression using the \code{prolfqua::build_model()} function.
#' @details
#'   Uses the \code{prolfqua::build_model()} function. Returns results sorted by p-value.
#'   Test is a comparison of nested models: \code{prolfqua::build_model()} fits
#'     full design and a reduced model. The reported F-test is of whether all
#'     columns removed in reduced model are zero. This is exact for a linear model, is 
#'     same hypothesis \code{test_proda()} tests by likelihood ratio, and equals Wald F-test 
#'     of same contrast. Any \code{config$test_term} is therefore testable, including
#'     one that the marginality rule spreads over several terms: testing \code{"grp"}
#'     in \code{~grp*sex} covers \code{grp} and \code{grp:sex} together, as it does
#'     under \code{test_method} \code{"lm"}, \code{"trend"} and \code{"voom"}.
#'   The error variance of each feature is shrunk toward a common prior across
#'     features before the F is formed, unless \code{config$test_moderate} is
#'     \code{FALSE}. The shrinkage is \code{limma::squeezeVar(robust=FALSE)}.
#'   The prior that variance is shrunk toward is one number shared by every feature
#'     unless the \code{trend} argument, or \code{config$test_trend} when it is not
#'     given, is \code{TRUE}, which fits it against the mean
#'     intensity of each feature as a natural spline of up to four degrees of freedom
#'     on the log residual variances. 
#'   A feature is dropped, with a warning, when its own data do not support the test.
#'   The \code{FDR} column of \code{hits} is \code{stats::p.adjust(..., "BH")} over
#'     the features tested, computed here rather than taken from \code{prolfqua}.
#'   With \code{mixed=TRUE}, which \code{test_method="prolfqua_lmer"} selects, the
#'     engine changes from \code{prolfqua::strategy_lm()} to
#'     \code{prolfqua::strategy_lmer()} and one model is fitted per gene over the rows
#'     of all of its features, rather than one model per feature. 
#'   The model is \code{intensity ~ 0 + <design columns> + (1|<feat_id_col>) +
#'     (1|<obs_col>)}. The random feature intercept absorbs the feature baselines,
#'     which differ by orders of magnitude among the precursors of one protein. 
#'   The random observation intercept is the error stratum the fixed effects belong to.
#'     Every covariate of \code{config$frm} comes from \code{state$samples}.
#'   The mixed test is a Wald F of one fit rather than a comparison of two.
#'   The error variance is not moderated across genes on this path, since the
#'     denominator of a Satterthwaite F is a combination of variance components rather
#'     than one residual variance and there is no single quantity to shrink.
#'   A gene with a single observed feature is fitted by \code{stats::lm()} instead: its
#'     random feature effect would be a single unknown confounded with the intercept, so
#'     that variance is not identified and \code{lme4} refuses the fit outright.
#'   A fit at the boundary of the parameter space, meaning a variance component
#'     estimated at zero, is kept and flagged in \code{hits$isSingular}.
#'   Flow is:
#'     \tabular{l}{
#'       1. Reshape data into long format with intensities, feature meta, and sample meta. \cr
#'       2. Build the design with \code{f.design_test_cols()} and add its columns to
#'            the long data as numeric columns. \cr
#'       3. Make \code{prolfqua::AnalysisConfiguration} and \code{prolfqua::LFQData}
#'            objects. \cr
#'       4. Make two \code{prolfqua::strategy_lm} objects, for the full design and for
#'            the design without the columns carrying \code{config$test_term}. \cr
#'       5. Build a \code{prolfqua} model from each. \cr
#'       6. Moderate the per-feature error variance across features with
#'            \code{limma::squeezeVar()}, toward a flat prior or, with
#'            \code{trend}, one fitted against mean feature intensity. \cr
#'       7. Return the per-feature F-test comparing the two fits, dropping features
#'            for which the test is not estimable. \cr
#'     }
#'   See documentation for \code{h0testr::new_config()} 
#'     for more detailed description of configuration parameters. 
#' @param state List with elements like those returned by \code{read_data()}:
#'   \tabular{ll}{
#'     \code{expression} \tab Numeric matrix with non-negative expression values. \cr
#'     \code{features}   \tab A data.frame with feature meta-data for rows of expression. \cr
#'     \code{samples}    \tab A data.frame with observation meta-data for columns of expression. \cr
#'   }
#' @param config List with configuration values. Uses the following keys:
#'   \tabular{ll}{
#'     \code{gene_id_col}           \tab Name of column in \code{state$features} with unique gene/protein group ids. \cr
#'     \code{feat_col}              \tab Name of column in \code{state$features} corresponding to \code{rownames(state$expression)}. \cr
#'     \code{obs_col}               \tab Name of column in \code{state$samples} corresponding to \code{colnames(state$expression)}. \cr
#'     \code{frm}                   \tab Formula (formula) to be fit. \cr
#'     \code{test_term}             \tab Term (character) to be tested for non-zero coefficient. \cr
#'     \code{contrast}      \tab Weighted sum (character scalar) of coefficients of \code{config$frm} to test instead of \code{config$test_term}; "" for none. \cr
#'     \code{reference_levels}      \tab Named character vector with the reference level of each factor variable in \code{config$frm}. \cr
#'     \code{covariate_types}       \tab Optional; classification of variables in \code{config$frm}, as set by \code{init_state()}. \cr
#'     \code{factor_levels}         \tab Optional; resolved levels of each factor variable, as set by \code{init_state()}. \cr
#'     \code{test_moderate}         \tab Optional logical; whether to shrink the error variance across features. Defaults to \code{TRUE} when absent. \cr
#'     \code{test_trend}            \tab Optional logical; whether the prior of that shrinkage is fitted against mean feature intensity rather than flat. Answers when the \code{trend} argument is not given; defaults to \code{FALSE} when both are absent. Unrelated to \code{config$test_method="trend"}. \cr
#'     \code{test_random_obs}       \tab Optional logical; whether the \code{mixed=TRUE} fit includes a random observation effect alongside the random feature effect. Defaults to \code{TRUE} when absent, which is the calibrated model; see Details. Ignored when \code{mixed=FALSE}. \cr
#'     \code{feat_id_col}           \tab Name of column in \code{state$features} with unique feature ids; must differ from \code{config$gene_id_col} when \code{mixed=TRUE}. \cr
#'     \code{normalization_method}  \tab If present and \code{is_log_transformed} unset, used to infer it. \cr
#'   }
#' @param is_log_transformed Logical scalar: whether \code{state$expression} has
#'   been log transformed. Defaults to \code{config$is_log_transformed}, which
#'   \code{h0testr::init_state()} and \code{h0testr::normalize()} maintain;
#'   passing both is an error unless they agree.
#' @param mixed Logical scalar: whether to fit one mixed model per gene over the rows
#'   of its features, with the feature and the observation as random effects, instead
#'   of one least squares model per feature. See Details.
#' @param trend Logical scalar: whether the variance prior of the moderation is fitted
#'   against mean feature intensity, which is \code{limma::eBayes(trend=TRUE)}'s
#'   covariate, rather than being flat. Defaults to \code{config$test_trend}, and to
#'   \code{FALSE} when that is absent.
#' @return
#'   A list with components:
#'   \tabular{ll}{
#'     \code{hits}  \tab \code{data.frame} of the per-feature F-tests of
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
#'     \code{fit}   \tab Model returned by \code{prolfqua::build_model()} for the
#'       full design. \cr
#'     \code{fit_reduced} \tab Model returned by \code{prolfqua::build_model()}
#'       for the design without the columns carrying \code{config$test_term}. \cr
#'     \code{design} \tab The design used, as returned by
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
#'     alongside it. \code{h0testr::test_h0()} fills its \code{logfc} column from the fitted
#'     coefficients instead, taken from \code{fit$model_df} or, with \code{mixed=TRUE}, from
#'     \code{coefs}: the signed coefficient when one column carries the test, which for a
#'     \strong{continuous} covariate is a change \strong{per unit} of it rather than a fold
#'     change between groups, and the total swing when several do. See \code{h0testr::test_h0()}.
#' @examples
#' pkgs <- c("prolfqua", "lme4", "lmerTest")
#' if(all(vapply(pkgs, requireNamespace, logical(1), quietly=TRUE))) {
#'   ## setup of expression data: ten peptides per gene, a third of them dropped, and no
#'   ##   missing values, so that the example is about the test rather than about missingness:
#'   set.seed(101)
#'   samps <- h0testr::sim_samples(factors=list(grp=c("ctl", "trt"), sex=c("F", "M")),
#'     n_per_cell=3)
#'   sim <- h0testr::sim_design(samps, frm=~grp + sex, test_term="grp", n_genes=20,
#'     peps_per_gene=10, p_drop=0.33, mnar_c0=-Inf, mnar_c1=0, mcar_p=0)
#'   state <- sim$state
#'   config <- sim$config
#'   rm(samps, sim)
#'
#'   out <- h0testr::init_state(state, config, minimal=TRUE)
#'
#'   ## actual test:
#'   result <- h0testr::test_prolfqua(out$state, out$config, is_log_transformed=FALSE)
#'   head(result$hits)
#'
#'   ## an interaction is fine too:
#'   config$frm <- ~grp*sex
#'   out <- h0testr::init_state(state, config, minimal=TRUE)
#'   result <- h0testr::test_prolfqua(out$state, out$config, is_log_transformed=FALSE)
#'   colnames(result$design$X)
#'   colnames(result$design$X)[result$design$cols_test]
#'   head(result$hits)
#'
#'   ## the same 2 df test with the error variance moderated across features:
#'   head(result$hits[, c("Df", "p.value", "p.value.unmod", "df.denom", "df.prior")])
#'   config$test_moderate <- FALSE
#'   out <- h0testr::init_state(state, config, minimal=TRUE)
#'   plain <- h0testr::test_prolfqua(out$state, out$config, is_log_transformed=FALSE)
#'   head(plain$hits[, c("Df", "p.value", "p.value.unmod", "df.denom", "df.prior")])
#'
#'   ## the prior can be fitted against mean feature intensity instead of being one
#'   ##   number, which is the prior test_trend() uses:
#'   config$test_moderate <- TRUE
#'   config$test_trend <- TRUE
#'   out <- h0testr::init_state(state, config, minimal=TRUE)
#'   trended <- h0testr::test_prolfqua(out$state, out$config, is_log_transformed=FALSE)
#'   head(trended$hits[, c("Df", "p.value", "moderated", "trend", "s2.prior")])
#'   range(trended$hits$s2.prior)
#'   range(result$hits$s2.prior)                     ## one value, the flat prior
#'
#'   ## mixed path, which test_method="prolfqua_lmer" selects: 
#'   config$test_moderate <- NULL
#'   config$test_trend <- NULL
#'   keep <- state$features$gene_id %in% unique(state$features$gene_id)[1:8]
#'   small <- list(expression=state$expression[keep, , drop=FALSE],
#'     features=state$features[keep, , drop=FALSE], samples=state$samples)
#'   out <- h0testr::init_state(small, config, minimal=TRUE)
#'   mix <- h0testr::test_prolfqua(out$state, out$config, is_log_transformed=FALSE,
#'     mixed=TRUE)
#'   nrow(mix$hits)                                  ## one row per gene, not per peptide
#'   head(mix$hits[, c("gene_id", "Df", "F.value", "p.value", "df.denom", "fit_type")])
#'
#'   ## df.denom is the Satterthwaite denominator, which the random observation effect
#'   ##   pulls down toward the number of observations:
#'   out$config$test_random_obs <- FALSE
#'   peponly <- h0testr::test_prolfqua(out$state, out$config,
#'     is_log_transformed=FALSE, mixed=TRUE)
#'   range(mix$hits$df.denom)
#'   range(peponly$hits$df.denom)
#' }
#' @export

test_prolfqua <- function(state, config, is_log_transformed=NULL, mixed=FALSE,
    trend=NULL) {

  f.need_pkgs(c("prolfqua", if(mixed) c("lme4", "lmerTest")),
    "test_prolfqua", config)

  ## mixed path models the features of a gene instead of aggregating them;
  ##   needs both levels:

  if(mixed && config$feat_id_col %in% config$gene_id_col) {
    f.err("test_prolfqua: the mixed path needs feature level input, so",
      "config$feat_id_col and config$gene_id_col must name different columns of",
      "state$features, and both are '", config$feat_id_col, "';", "\n",
      "with one feature per gene the random feature effect is not identified and the",
      "model reduces to the one test_method 'prolfqua' fits, so use that instead, or",
      "supply un-aggregated data carrying a gene id column", config=config)
  }

  check_config(config)
  f.check_state(state, config)

  is_log_transformed <- f.is_log_transformed(is_log_transformed, config,
    "test_prolfqua")

  ## argument overrides config$test_trend:
  trend <- f.is_trend(trend, config, "test_prolfqua")
  parsed <- f.parse_frm(config$frm, config)
  idvars <- unique(c(config$gene_id_col, config$feat_id_col))

  ## ids mixed path groups:
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
  
  ## the design matrix, built explicitly so that the reduced model can be 
  ##   formed by dropping columns from it:

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

  ## covariates reach fit as numeric columns of design matrix:

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

  ## Reduced model. For config$test_term, a subset of the columns of X:

  cols_red <- character(0)

  if(mixed) {

    ## mixed path tests one fit by Wald rather than comparing two:

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

  meta <- prolfqua::AnalysisConfiguration$new()
  meta$work_intensity <- "intensity"
  meta$is_response_transformed <- is_log_transformed
  meta$sample_name <- "sample"
  meta$hierarchy[[config$gene_id_col]] <- config$gene_id_col
  meta$hierarchy[[config$feat_id_col]] <- config$feat_id_col

  obj <- prolfqua::LFQData$new(data=dat, config=meta)

  ## mixed path takes over here: one model per gene over the rows of its features,
  ##   with feature and the observation entering as random effects:

  if(mixed) return(f.prolfqua_mixed(obj, design, cols, config, trend=trend))

  ## comparison of the full design against design with 
  ##   columns carrying config$test_term removed:
  ## config$contrast constrains model rather than dropping terms from it:

  frm_full <- paste("intensity ~ 0 +", paste(cols, collapse=" + "))
  frm_red <- if(length(cols_red)) {
    paste("intensity ~ 0 +", paste(cols_red, collapse=" + "))
  } else {
    "intensity ~ 0"
  }

  f.msg("test_prolfqua:", f.test_label(design, config), "; full:", frm_full, "\n",
    "  reduced:", frm_red, "; numerator df:", design$df_intend, config=config)

  fit_full <- prolfqua::build_model(data=obj$data_long(),
    model_strategy=prolfqua::strategy_lm(frm_full),
    subject_id=obj$hierarchy_keys())

  fit_red <- prolfqua::build_model(data=obj$data_long(),
    model_strategy=prolfqua::strategy_lm(frm_red),
    subject_id=obj$hierarchy_keys())

  ## config$test_trend fits prior variance of moderation against mean feature
  ##   intensity instead of shrinking every feature toward one number:

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

  ## nrow(), not length(unique(tbl[[config$feat_col]])): f.prolfqua_nested_f() keys
  ##   its table by unique(c(config$gene_id_col, config$feat_id_col)):

  f.msg("tested", nrow(tbl), "features; found",
    sum(tbl$FDR < 0.05, na.rm=T), "hits", config=config)

  return(list(hits=tbl, fit=fit_full, fit_reduced=fit_red, design=design))
}

