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
#'   Reachable only by calling this function directly: \code{h0testr::test()} passes no
#'   value for it and no \code{config} key carries one, so a workflow run gets the default.
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

  ## the aggregate path models the feature level, so there has to be one. With one
  ##   feature per gene the random feature effect is a single unknown confounded with
  ##   the intercept, its variance is not identified, and the model reduces to the one
  ##   aggregate=FALSE fits, so that is what to use. Ahead of the checks below for the
  ##   reason given in test_prolfqua(), whose mixed path refuses the same thing in the
  ##   same words: this reads config alone, and it is the refusal that names the
  ##   method to use instead. The usual way to arrive here is config$run_order naming
  ##   "combine_features", which sets config$feat_id_col to config$gene_id_col, there being
  ##   one row per gene afterwards and nothing below it; so the fix is to drop that step
  ##   from run_order rather than to rename an id column. h0testr::tune() already skips it
  ##   for the methods that take feature level input:

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

  ## the colData() loop below already relevels the covariates msqrob2::msqrob() fits
  ##   from, but the state handed to f.design_test_cols() further down was not, so on
  ##   a direct call the two designs disagreed about which level is the reference and
  ##   the cols_pick %in% parms check turned that into an error: a call that should
  ##   have worked failed instead of returning the wrong sign, which is the better of
  ##   the two failures but still not the right one. Both designs are built from a
  ##   releveled state now. A no-op for a run that came through test(); see test_lm():

  state <- f.relevel_state_covariates(state, config, caller="test_msqrob")

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
  ##   with the level ordering config resolves: a character column there is
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

