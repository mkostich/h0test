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
#'   Coefficients carrying that test are columns of the design matrix
#'     assigned to \code{config$test_term} and to every term containing it. Naming a 
#'     variable that also appears in an interaction tests the interaction too: with
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
#'     \code{expression} \tab Numeric matrix with non-negative expression values. \cr
#'     \code{features}   \tab A data.frame with feature meta-data for rows of expression. \cr
#'     \code{samples}    \tab A data.frame with observation meta-data for columns of expression. \cr
#'   } 
#' @param config List with configuration values. Uses the following keys:
#'   \tabular{ll}{
#'     \code{feat_col}       \tab Name of column in \code{feature_file_in} corresponding to \code{rownames(state$expression)}. \cr
#'     \code{obs_col}        \tab Name of column in \code{sample_file_in} corresponding to \code{colnames(state$expression)}. \cr
#'     \code{frm}            \tab Formula (formula) to be fit. \cr
#'     \code{test_term}      \tab Term (scalar character) to be tested for non-zero coefficient. \cr
#'     \code{contrast}      \tab Weighted sum (character scalar) of coefficients of \code{config$frm} to test instead of \code{config$test_term}; "" for none. \cr
#'     \code{reference_levels} \tab Named character vector with the reference level of each factor variable in \code{config$frm} (see examples). \cr
#'   }
#' @return 
#'   A list with components:
#'   \tabular{ll}{
#'     \code{hits}  \tab \code{data.frame} of results from \code{limma::topTable()}. \cr
#'     \code{fit}   \tab Model returned by \code{limma::eBayes()}. \cr
#'   } 
#'   \code{logFC} is an effect size on the scale of \code{state$expression}: for a two level
#'     factor, the difference between its levels, or log fold change when the input is log
#'     transformed; for a \strong{continuous} covariate, the change \strong{per unit} of it,
#'     whose size depends on the units the covariate is recorded in. This is what
#'     \code{h0testr::test_h0()} reports as \code{logfc}. A joint test has no \code{logFC}
#'     column at all, having no single fold change to report, and \code{h0testr::test_h0()}
#'     reports the total swing instead. See \code{h0testr::test_h0()}.
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
#' out <- h0testr::init_state(state, config, minimal=TRUE)
#' 
#' tbl <- h0testr::test_trend(out$state, out$config)
#' print(tbl)
#' @export

test_trend <- function(state, config) {

  check_config(config)
  f.check_state(state, config)

  if(!is.matrix(state$expression)) {
    f.err("test_trend: !is.matrix(state$expression)", config=config)
  }

  ## declared reference level first in each factor covariate, before design is built:

  state <- f.relevel_state_covariates(state, config, caller="test_trend")
  
  ## design and columns of it carrying test, from same helper
  ##   test_lm() and filter_features_by_estimability() use:

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

