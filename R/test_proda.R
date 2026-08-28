#' Hypothesis testing using the \code{proDA} package
#' @description
#'   Tests for differential expression using the \code{proDA::proDA()} function.
#' @details
#'   Uses the \code{proDA::proDA()} function. Returned results sorted by p-value.
#'   Returns peptide/precursor/gene-level hypothesis testing results based on 
#'     peptide/precursor/gene-level input. Natively handles missing values.
#'   When \code{config$test_term} resolves to a single design matrix column, that
#'     column is tested by \code{proDA::test_diff()} as a single contrast, and
#'     \code{logFC} is the corresponding coefficient: for a two level factor, the non
#'     reference level minus the reference level declared in
#'     \code{config$reference_levels}. When it resolves to several columns, all of them are 
#'     tested jointly, as a likelihood ratio test of the full model against a reduced model. 
#'     Reported as an F statistic with \code{logFC} \code{NA}, as for the F tests of
#'     \code{h0testr::test_trend()} and \code{h0testr::test_voom()}.
#'   Testing the intercept (\code{config$test_term} \code{"1"}) is supported:
#'     \code{proDA::proDA()} renames that column \code{Intercept}, and the contrast is
#'     looked up under that name.
#'   The one formula this method cannot take that the others can is a
#'     \code{config$test_term} that leaves no parameters in the reduced model.
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
#'   \code{h0testr::init_state()} and \code{h0testr::normalize()} maintain;
#'   passing both is an error unless they agree.
#' @param prior_df Strictly positive count (\code{location_prior_df}) indicating number of dfs for prior.
#' @param maxit Strictly positive count indicating maximum number of iterations for \code{proDA::proDA()} algorithm.
#'   Reachable only by calling this function directly: \code{h0testr::test_h0()} passes no
#'   value for it and no \code{config} key carries one.
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
#'     \code{h0testr::test_h0()} reports as \code{logfc}. The likelihood ratio test reports no
#'     effect size at all, having no single contrast to report, and \code{h0testr::test_h0()}
#'     reports the total swing instead. See \code{h0testr::test_h0()}.
#' @examples
#' pkgs <- c("proDA", "SummarizedExperiment")
#' if(all(vapply(pkgs, requireNamespace, logical(1), quietly=TRUE))) {
#'   ## setup of expression data: ten peptides per gene, a third of them dropped, and no
#'   ##   missing values, so that the example is about the test rather than about missingness:
#'   set.seed(101)
#'   samps <- h0testr::sim_samples(factors=list(grp=c("ctl", "trt"), sex=c("F", "M")),
#'     n_per_cell=3)
#'   sim <- h0testr::sim_design(samps, frm=~grp + sex, test_term="grp", n_genes=100,
#'     peps_per_gene=10, p_drop=0.33, mnar_c0=-Inf, mnar_c1=0, mcar_p=0)
#'   state <- sim$state
#'   config <- sim$config
#'   rm(samps, sim)
#'
#'   ## no grp:sex term here, so testing "grp" is the single coefficient grptrt, and a
#'   ##   fold change is reported:
#'   out <- h0testr::init_state(state, config, minimal=TRUE)
#'
#'   ## actual test:
#'   result <- h0testr::test_proda(out$state, out$config, is_log_transformed=FALSE)
#'   head(result$hits)
#' }
#' @export

test_proda <- function(state, config, is_log_transformed=NULL, prior_df=3, maxit=20) {

  check_config(config)

  f.need_pkgs(c("proDA", "SummarizedExperiment"), "test_proda", config)
  f.check_state(state, config)

  is_log_transformed <- f.is_log_transformed(is_log_transformed, config,
    "test_proda")

  ## proDA::proDA() builds its own design from col_data and the formula:

  state <- f.relevel_state_covariates(state, config, caller="test_proda")
  design <- f.design_test_cols(state, config)
  cols_pick <- colnames(design$X)[design$cols_test]

  fit <- proDA::proDA(state$expression, design=design$parsed$frm, col_data=state$samples,
    data_is_log_transformed=is_log_transformed, location_prior_df=prior_df, max_iter=maxit)

  ## proDA::proDA() renames '(Intercept)' to 'Intercept':

  cols_want <- colnames(design$X)
  cols_want[cols_want %in% "(Intercept)"] <- "Intercept"
  cols_have <- colnames(proDA::design(fit))

  if(!identical(cols_want, cols_have)) {
    f.err("test_proda: the design proDA::proDA() built does not match the design",
      "h0testr derived from config$frm;", "\n", "proDA::design(fit):", cols_have,
      "\n", "h0testr:", cols_want, config=config)
  }

  ## config$contrast goes to the likelihood ratio branch below whatever it weights.
  ##   proDA::test_diff() takes a contrast expression, so one df Wald test is available:

  if(is.null(design$contrast) && length(design$cols_test) %in% 1) {

    ## one design column carries the test, so the Wald test on that coefficient:

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

    ## several design columns carry the test, which no single contrast can express:
    x_red <- design$X_red

    ## no columns left in the reduced model:
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

    ## proDA::test_diff() requires a full rank reduced model;
    ##   reports rank deficient one as colinear:

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

