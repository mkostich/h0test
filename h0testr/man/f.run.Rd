% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/test.R
\name{f.run}
\alias{f.run}
\title{Run a basic workflow}
\usage{
f.run(config)
}
\arguments{
\item{config}{List with configuration values like those returned by \code{f.new_config()}.}
}
\value{
A list with the following elements:
  \tabular{ll}{
    \code{state}  \cr \tab A list with elements \code{$expression}, \code{$features}, and \code{$samples}. \cr
    \code{config} \cr \tab A list with configuration settings. \cr
    \code{tbl}    \cr \tab A data.frame containing results of test.. \cr
  }
}
\description{
Run a basic workflow according to: 
    \code{config$run_order}.
}
\details{
Run a basic workflow: 
    \code{f.load_data -> config$run_order -> f.test}, where 
      \code{config$run_order} is vector of functions which are run in
      the specified order.
}
\examples{
config <- h0testr::f.new_config()
config$dir_in <- system.file("extdata", package="h0testr")  ## where example data 
config$feature_file_in <- "features.tsv"
config$sample_file_in <- "samples.tsv"
config$data_file_in <- "expression.tsv" 
config$feat_id_col <- "feature_id"
config$obs_id_col <- "observation_id"
config$sample_id_col <- "observation_id"
config$frm <- ~condition
config$test_term <- "condition"
config$test_method <- "trend"
config$sample_factors <- list(condition=c("placebo", "drug"))
config$n_features_min <- 10      ## default 1000 too big for small demo dataset
config$run_order <- c("normalize", "combine_reps", "filter", "impute")
config$save_state <- FALSE

print(config$run_order)

out <- h0testr::f.run(config)   ## run workflow
print(out$tbl)                  ## hit table
}
