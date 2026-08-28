#' Hypothesis testing using the \code{msqrob2} package
#' @description
#'   Tests for differential expression using the \code{msqrob2::msqrob()} function, or
#'   with \code{aggregate=TRUE} the \code{msqrob2::msqrobAggregate()} function.
#' @details
#'   Returns gene-level hypothesis testing results based on peptide/precursor-level input.  
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
#'     \code{msqrob2::msqrobAggregate()}. That function aggregates to make the assay
#'     that carries the result, one row per gene, but fits \code{msqrob2::msqrobLmer()}
#'     on the un-aggregated assay with the features of a gene as one group, so 
#'     models are feature level and results gene level.
#'   The model \code{aggregate=TRUE} fits is
#'     \code{intensity ~ <design columns> + (1|config$feat_id_col) +
#'     (1|config$obs_col)}, the second random term dropped when
#'     \code{config$test_random_obs} is \code{FALSE}.
#'   \code{config$test_ridge=TRUE} penalizes the fixed effects,
#'     \code{msqrob2::msqrobAggregate(ridge=TRUE)}. Is \code{FALSE} by default.
#'   A gene with fewer than two observed features is fitted at the gene level instead,
#'     by \code{msqrob2::msqrob()} on aggregated values.
#'   \code{hits$nNonZero} is the number of observations in which the gene was measured
#'     at all, meaning in which at least one of its features has a value. 
#'   \code{msqrob2::hypothesisTest()} returns one table per contrast rather than a
#'     joint test over several, so it can answer only a \code{config$test_term} that
#'     resolves to a single design matrix column. But each fitted
#'     \code{msqrob2} \code{StatModel} carries the coefficients, the unscaled
#'     covariance of the design, and the moderated variance and posterior degrees of
#'     freedom. For a test spanning several columns the joint statistic is computed
#'     from \code{W}, the quadratic form of the coefficients under test in the inverse 
#'     of their covariance (unscaled covariance of the columns times the posterior variance). 
#'     Reported as \code{F = W / length(cols)} on \code{length(cols)} with \code{dfPosterior}
#'     degrees of freedom.
#'   \strong{That joint statistic is computed by \code{h0testr}, not returned by
#'     \code{msqrob2}}. With one column the
#'     statistic reduces to the square of the moderated t that \code{msqrob2} reports,
#'     which is why it is referred to an F rather than a chi-square.
#'   A gene whose model could not be fit is reported as \code{NA}.
#'   The reported contrast is the tested coefficient under the factor level ordering
#'     \code{config$reference_levels} declares and \code{h0testr::init_state()}
#'     resolves, so for a two level factor \code{logFC} is the non-reference level
#'     minus the reference level. \code{msqrob2::msqrob()} builds its own design, so
#'     that ordering is carried through to it explicitly.
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
#'     \code{covariate_types}  \cr \tab Optional; classification of variables in \code{config$frm}, as set by \code{init_state()}. \cr
#'     \code{factor_levels}    \cr \tab Optional; resolved levels of each factor variable, as set by \code{init_state()}. \cr
#'     \code{feat_id_col}      \cr \tab Name of column in \code{state$features} with unique feature ids; must differ from \code{config$gene_id_col} when \code{aggregate=TRUE}. \cr
#'     \code{test_random_obs}  \cr \tab Optional logical; whether the \code{aggregate=TRUE} fit includes a random observation effect alongside the random feature effect. Defaults to \code{TRUE} when absent, which is the calibrated model; see Details. Ignored when \code{aggregate=FALSE}. \cr
#'     \code{test_ridge}       \cr \tab Optional logical; whether the \code{aggregate=TRUE} fit penalizes the fixed effects. Defaults to \code{FALSE} when absent, which is \code{msqrob2}'s own default; see Details. Ignored when \code{aggregate=FALSE}. \cr
#'   }
#' @param maxit Integer scalar >= 1. How many iterations to use for \code{rlm} fitting.
#' @param aggregate Logical scalar. \code{FALSE}, the default, aggregates the features of
#'   a gene and fits the aggregated values, \code{msqrob2::msqrob()}. \code{TRUE} fits
#'   one mixed model per gene over the rows of its features instead,
#'   \code{msqrob2::msqrobAggregate()}, which is what
#'   \code{config$test_method="msqrob_agg"} selects.
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
#'     whose size depends on the units the covariate is recorded in. The joint test reports no 
#'     effect size, since no single contrast to report, and \code{h0testr::test_h0()} reports the 
#'     total swing instead. See \code{h0testr::test_h0()}.
#' @examples
#' set.seed(101)
#' samps <- h0testr::sim_samples(factors=list(grp=c("ctl", "trt"), sex=c("F", "M")),
#'   n_per_cell=3)
#' sim <- h0testr::sim_design(samps, frm=~grp + sex, test_term="grp", n_genes=100,
#'   peps_per_gene=10, p_drop=0.33, mnar_c0=-Inf, mnar_c1=0, mcar_p=0)
#' state <- sim$state
#' config <- sim$config
#' rm(samps, sim)
#'
#' out <- h0testr::init_state(state, config, minimal=TRUE)
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

  ## the aggregate path models the feature level. With one feature per gene the random 
  ##   feature effect is a single unknown confounded with the intercept, its variance 
  ##   is not identified, and the model reduces to the one aggregate=FALSE fits:

  if(aggregate && config$feat_id_col %in% config$gene_id_col) {
    f.err("test_msqrob: the aggregate path needs feature level input, so",
      "config$feat_id_col and config$gene_id_col must name different columns of",
      "state$features, and both are '", config$feat_id_col, "';", "\n",
      "with one feature per gene the random feature effect is not identified and the",
      "model reduces to the one test_method 'msqrob' fits, so use that instead, or",
      "supply un-aggregated data carrying a gene id column", config=config)
  }

  check_config(config)
  f.check_state(state, config)
  state <- f.relevel_state_covariates(state, config, caller="test_msqrob")
  exprs <- as.data.frame(state$expression)

  n_non0 <- apply(exprs, 1, function(v) sum(!is.na(v)))
  exprs <- cbind(fnames=rownames(exprs), exprs)
  rownames(exprs) <- NULL
  obj <- QFeatures::readQFeatures(table=exprs, ecol=2:ncol(exprs), 
    fnames="fnames", name="features")  

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
  ##   with the level ordering config resolves: 

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
  ##   test column from a design built on parsed formula, whose interaction
  ##   labels have their variables sorted:

  parsed <- f.parse_frm(config$frm, config)

  ## one gene level fit of aggregated values, or one feature level mixed model per gene; 
  ##   the models sit in rowData(obj[["genes"]])$msqrobModels, so everything below is shared:

  fit_type <- NULL

  if(aggregate) {

    fitted <- f.msqrob_agg(obj, state, config, parsed$frm, maxit)
    obj <- fitted$obj
    fit_type <- fitted$fit_type

  } else {

    obj <- QFeatures::aggregateFeatures(obj, i="features",
      fcol=config$gene_id_col, na.rm=T, name="genes", fun=base::colMeans)

    obj <- msqrob2::msqrob(object=obj, i="genes", formula=parsed$frm, maxitRob=maxit)
  }

  ## design and columns carrying the test. msqrob2::hypothesisTest() returns
  ##   one table per contrast rather than joint test over several, so it can answer
  ##   only a test that resolves to one coefficient:

  design <- f.design_test_cols(state, config)
  cols_pick <- colnames(design$X)[design$cols_test]

  ## msqrob2::hypothesisTest() matches the contrast to the fit by parameter name;
  ##   name fit does not have yields a column of NAs instead of an error:

  parms <- colnames(stats::model.matrix(parsed$frm,
    data=as.data.frame(SummarizedExperiment::colData(obj))))

  if(!all(cols_pick %in% parms)) {
    f.err("test_msqrob: design matrix column", cols_pick, "carrying the test of",
      f.test_label(design, config), "is not among the parameters",
      "msqrob2::msqrob() fit;", "\n", "parameters fit:",
      paste(parms, collapse=", "), config=config)
  }

  ## config$test_ridge renames every fitted fixed effect:

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

  ## number of observations each gene measured in, and number of features behind it:

  dat <- SummarizedExperiment::rowData(obj[["genes"]])
  dat <- as.data.frame(dat[, c(config$gene_id_col, ".n"), drop=F])
  dat$nNonZero <- as.numeric(f.msqrob_n_obs(state, config)[dat[[config$gene_id_col]]])
  dat <- dat[, c(config$gene_id_col, "nNonZero", ".n"), drop=F]

  if(!is.null(fit_type)) {
    dat$fit_type <- as.character(fit_type[dat[[config$gene_id_col]]])
  }

  ## config$contrast is a weighted sum of coefficients, which is one contrast;
  ##   msqrob2::hypothesisTest() moderated t rather than the joint F below:

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

  ## one column: the package's own moderated t, unchanged. 
  ##   more than one: joint Wald test over the same fitted models:

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

## helper for test_msqrob(): joint test over several coefficients that
##   msqrob2::hypothesisTest() cannot express. 

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

    ## msqrob2 leaves covariance of a mixed fit as a Matrix rather than as a base matrix:

    vcv <- try(as.matrix(vcv), silent=T)
    if(inherits(vcv, "try-error")) next
    if(!is.matrix(vcv) || !all(cols %in% rownames(vcv))) next
    if(length(sig) != 1 || is.na(sig) || !is.finite(sig) || sig <= 0) next
    if(length(dfp) != 1 || is.na(dfp) || dfp <= 0) next    ## dfp may be Inf

    ## solve() rather than a pseudo-inverse: a singular covariance over the tested
    ##   columns means the data cannot separate them for this gene,:

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

## helper for test_msqrob(): names under which msqrob2 fitted the
##   design. msqrob2 penalizes the fixed effects, ridge=TRUE, by refitting them as a
##   random effect over a grouping factor it calls "ridge" whose levels are the design
##   columns, and lme4 names a random effect estimate by pasting the name of the
##   grouping factor in front of the level:

f.msqrob_parms <- function(cols, ridge) {
  if(!isTRUE(ridge)) return(cols)
  return(ifelse(cols %in% "(Intercept)", cols, paste0("ridge", cols)))
}

## helper for test_msqrob(): the number of observations in which each gene was measured
##   at all, meaning in which at least one of its features has a value:

f.msqrob_n_obs <- function(state, config, who="test_msqrob") {
  genes <- f.gene_ids(state$features, config, who)
  ok <- !is.na(state$expression)
  n <- tapply(1:nrow(ok), genes, function(idxs) {
    sum(apply(ok[idxs, , drop=F], 2, any))
  })
  return(n)
}

## helper for test_msqrob(): aggregate=TRUE fit. msqrob2::msqrobAggregate() is not
##   an alternative to aggregating so much as a different use of it:

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

  ## aggregated assay is named "genes" for both paths:

  obj <- f.quiet_fits(msqrob2::msqrobAggregate(object=obj, i="features",
    fcol=config$gene_id_col, name="genes", formula=frm_mix, ridge=ridge, robust=T,
    maxitRob=maxit, modelColumnName="msqrobModels",
    aggregateFun=function(x, ...) base::colMeans(x, na.rm=T)),
    "mixed", config, "test_msqrob")

  dat <- SummarizedExperiment::rowData(obj[["genes"]])
  models <- dat$msqrobModels
  ids <- as.character(dat[[config$gene_id_col]])
  fit_type <- vapply(models, msqrob2::getFitMethod, character(1))

  ## genes with <2 observed features are fitted at gene level instead:

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

