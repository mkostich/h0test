% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/config.R
\name{f.report_config}
\alias{f.report_config}
\title{Report configuration}
\usage{
f.report_config(config)
}
\arguments{
\item{config}{List with configuration values.}
}
\description{
Reports configuration settings used for run, and checks if
  \code{config$test_term} is compatible with \code{config$frm}.
}
\details{
Report written to \code{config$log_file}; if \code{config$log_file == ""},
  written to standard out (console or terminal). Only supports non-lists
    and lists of non-lists (not lists of lists) as \code{config} values.
  Throws error if \code{config$test_term} is not compatible with \code{config$frm}.
}
\examples{
config <- f.new_config()

## you must customize frm, test_term, and sample_factors:
config$frm <- ~ age + sex + age:sex
config$test_term <- "age:sex"
config$sample_factors <- list(age=c("young", "old"), sex=c("female", "male"))

## these may or may not need customization, depending on your file formats:
config$feat_id_col <- "gene"
config$sample_id_col <- "sample"
config$obs_id_col <- "observation"

config$log_file <- ""          ## log to console/stdout
f.report_config(config)
}
