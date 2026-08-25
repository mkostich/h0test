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
#'   Reachable only by calling this function directly: \code{h0testr::test()} passes no
#'   value for it and no \code{config} key carries one, so a workflow run gets the default.
#'   Unlike \code{prior_df}, which \code{config$test_prior_df} answers.
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

  check_config(config)
  f.check_state(state, config)

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

  ## proDA::proDA() builds its own design from col_data and the formula, so the
  ##   covariates have to reach it with the level ordering config resolves, exactly
  ##   as they reach f.design_test_cols() below; otherwise the two designs differ in
  ##   which level is the reference, which the column name check further down
  ##   reports rather than letting it through. A no-op for a run that came through
  ##   test(); this function is exported, so it does not rely on that:

  state <- f.relevel_state_covariates(state, config, caller="test_proda")

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

