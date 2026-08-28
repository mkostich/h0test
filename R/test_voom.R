#' Hypothesis testing using \code{limma::voom}
#' @description
#'   Tests for differential expression using the \code{limma::voom()} function.
#' @details
#'   The \code{limma::voom()} model is fit to \code{config$frm} and a test
#'     is performed for whether the effect of \code{config$test_term} on
#'     \code{state$expression} is zero.
#'   Naming a variable that also appears in an interaction gives a
#'     joint test over the interaction as well. One coefficient gives a moderated
#'     t-test with a \code{logFC} column, several give an F-test with an \code{F}
#'     column.
#'   Note \code{limma::voom()} models a count mean-variance relationship, so it is
#'     appropriate for count-like input rather than for already log-transformed
#'     abundances.
#'   \code{limma::voom()} fits one mean-variance trend across the whole matrix and
#'     has no handling for a missing value, so a feature carrying any is held out of
#'     the fit and does not appear in the result. How many were held
#'     out, and the first few of them, are logged. If fewer than two features are
#'     left, the call is refused: impute first, see \code{h0testr::impute()}, or use a
#'     method that does not need complete features, such as \code{h0testr::test_lm()}.
#'   See documentation for \code{h0testr::new_config()} for detailed description of 
#'     configuration parameters.
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
#'     \code{contrast}      \cr \tab Weighted sum (character scalar) of coefficients of \code{config$frm} to test instead of \code{config$test_term}; "" for none. \cr
#'     \code{reference_levels} \cr \tab Named character vector with the reference level of each factor variable in \code{config$frm} (see examples). \cr
#'   }
#' @param normalize.method Character in
#'   \code{c("TMM", "TMMwsp", "RLE", "upperquartile", "none")}. Reachable only by calling
#'   this function directly: \code{h0testr::test()} passes no value for it and no
#'   \code{config} key carries one, so a workflow run gets the default. The
#'   \code{h0testr::normalize()} step is where a workflow scales its data.
#' @return 
#'   A list with components:
#'   \tabular{ll}{
#'     \code{hits}  \cr \tab \code{data.frame} of results from \code{limma::topTable()}. \cr
#'     \code{fit}   \cr \tab Model returned by \code{limma::eBayes()}. \cr
#'   } 
#'   \code{logFC} is an effect size on the scale of \code{state$expression}: for a two level
#'     factor, the difference between its levels, so a log fold change when the input is log
#'     transformed; for a \strong{continuous} covariate, the change \strong{per unit} of it,
#'     whose size depends on the units the covariate is recorded in. This is what
#'     \code{h0testr::test()} reports as \code{logfc}. A joint test has no \code{logFC}
#'     column at all, having no single fold change to report, and \code{h0testr::test()}
#'     reports the total swing instead. See \code{h0testr::test()}.
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
#' ## set up and check configuration, including covariates:
#' out <- h0testr::initialize(state, config, minimal=TRUE)
#' 
#' result <- h0testr::test_voom(out$state, out$config)
#' head(result$hits)

test_voom <- function(state, config, normalize.method="none") {

  check_config(config)
  f.check_state(state, config)
  
  ## declared reference level first in each factor covariate, before design built:
  state <- f.relevel_state_covariates(state, config, caller="test_voom")

  exprs <- state$expression
  i <- apply(exprs, 1, function(v) any(is.na(v)))

  ## limma::voom() fits one mean-variance trend across the whole matrix;
  ##   has no handling for a missing value:

  if(any(i)) {
    f.msg("WARNING: test_voom: dropping", sum(i), "of", nrow(exprs), "features with",
      "at least one missing value, which limma::voom() cannot fit;", "\n",
      "first few:", paste(utils::head(rownames(exprs)[i], 5), collapse=", "),
      config=config)
  }

  if(sum(!i) < 2) {
    f.err("test_voom: only", sum(!i), "of", nrow(exprs), "features have no missing",
      "value, and limma::voom() needs at least two to fit a mean-variance trend;",
      "impute first, see h0testr::impute(), or use a method that does not need",
      "complete features, such as h0testr::test_lm()", config=config)
  }

  exprs <- exprs[!i, , drop=F]

  design <- f.design_test_cols(state, config)

  obj <- limma::voom(exprs, design$X, plot=F, normalize.method=normalize.method)
  fit <- limma::lmFit(obj, design$X)
  lc <- f.limma_contrast_fit(fit, design, config)
  fit <- limma::eBayes(lc$fit, trend=F)

  ## single coefficient yields t-test and logFC column; several yield F-test:

  tbl <- limma::topTable(fit, coef=lc$coef, number=Inf)

  f.msg("test_voom:", f.test_label(design, config), "; design columns:",
    ncol(design$X), "; test columns:", length(design$cols_test), "; df:",
    design$df_intend, config=config)
  f.msg("tested", nrow(exprs), "of", nrow(state$expression), "features",
    config=config)
  f.msg("found", sum(tbl$adj.P.Val < 0.05, na.rm=T), "hits", config=config)
  
  return(list(hits=tbl, fit=fit))
}

