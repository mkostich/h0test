% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/test.R
\name{f.test}
\alias{f.test}
\title{Hypothesis testing}
\usage{
f.test(state, config)
}
\arguments{
\item{state}{List with elements formatted like the list returned by `f.read_data()`:
\tabular{ll}{
  \code{expression} \cr \tab Numeric matrix with non-negative expression values. \cr
  \code{features}   \cr \tab A data.frame with feature meta-data for rows of expression. \cr
  \code{samples}    \cr \tab A data.frame with observation meta-data for columns of expression. \cr
}}

\item{config}{List with configuration values. Uses the following keys:
\tabular{ll}{
  \code{feat_id_col}    \cr \tab Column (character scalar) in \code{feature_file_in} that corresponds to rows of \code{data_file_in}. \cr
  \code{obs_id_col}     \cr \tab Column (character scalar) in \code{sample_file_in} that corresponds to columns of \code{data_file_in}. \cr
  \code{sample_id_col}  \cr \tab Column (character scalar) in \code{sample_file_in} with unique sample labels. \cr
  \code{frm}            \cr \tab Formula (formula) to be fit \cr
  \code{test_term}      \cr \tab Term (character scalar) to be tested for non-zero coefficient. \cr
  \code{sample_factors} \cr \tab List specifying levels of factor variables in \code{config$frm} (see examples). \cr
  \code{test_method}    \cr \tab Character scalar in \code{c("trend", "voom")}. \cr
}}
}
\value{
A data.frame containing results of test.
}
\description{
Test hypotheses using .
}
\details{
Tests for differential expression using method specified in config.
}
\examples{
set.seed(101)
## no missing values: mnar_c0=-Inf, mnar_c1=0, mcar_p=0
exprs <- h0testr::f.sim2(n_obs1=6, n_obs2=6, n_feats=25, n_sig=5, fold_change=2, 
  mnar_c0=-Inf, mnar_c1=0, mcar_p=0)
exprs <- log2(exprs + 1)
feats <- data.frame(feature_id=rownames(exprs))
samps <- data.frame(observation_id=colnames(exprs), 
  condition=c(rep("placebo", 6), rep("drug", 6)))
state <- list(expression=exprs, features=feats, samples=samps)

config <- h0testr::f.new_config()
config$feat_id_col <- "feature_id"
config$obs_id_col <- "observation_id"
config$sample_id_col <- "observation_id"
config$frm <- ~condition
config$test_term <- "condition"
config$sample_factors <- list(condition=c("placebo", "drug"))
config$test_method <- "trend"

## set up and check covariates:
out <- h0testr::f.preprocess_covariates(state, config)

## save_state=FALSE so test results not written to file:
config <- list(
  log_file="", frm=~condition, test_term="condition", 
  sample_factors=list(condition=c("placebo", "drug")), test_method="trend", 
  feat_id_col="feature_id", save_state=FALSE
)

out <- h0testr::f.test(state, config)
print(out)
}
