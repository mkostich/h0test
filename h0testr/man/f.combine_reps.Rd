% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/normalize.R
\name{f.combine_reps}
\alias{f.combine_reps}
\title{Combine technical replicates}
\usage{
f.combine_reps(state, config)
}
\arguments{
\item{state}{List with elements formatted like the list returned by \code{f.read_data()}:
\tabular{ll}{
  \code{expression} \cr \tab Numeric matrix with non-negative expression values. \cr
  \code{features}   \cr \tab A data.frame with feature meta-data for rows of expression. \cr
  \code{samples}    \cr \tab A data.frame with observation meta-data for columns of expression. \cr
}}

\item{config}{List with configuration values. Uses the following keys:
\tabular{ll}{
  \code{log_file}      \cr \tab Path to log file (character); \code{log_file=""} outputs to console. \cr
  \code{obs_id_col}    \cr \tab Column in observation metadata matching rownames of \code{state$expression} \cr
  \code{sample_id_col} \cr \tab Column in observation metadata with unique sample ids. \cr
}}
}
\value{
A list (the processed state) with the following elements:
  \tabular{ll}{
    \code{expression} \cr \tab Numeric matrix with non-negative expression values. \cr
    \code{features}   \cr \tab A data.frame with feature meta-data for rows of expression. \cr
    \code{samples}    \cr \tab A data.frame with observation meta-data for columns of expression. \cr
  }
}
\description{
Combine expression signals from technical replicates.
}
\details{
Combines signals for each feature across technical replicates by taking median.
}
\examples{
  \dontrun{
    config <- list(obs_id_col="observation_id", sample_id_col="sample_id", log_file="")
    state <- list(expression=exprs, samples=samps, features=feats)
    cat("config$obs_col:", config$obs_col, "\n")
    out <- f.combine_reps(state, config)
    config2 <- out$config
    cat("config2$obs_col:", config2$obs_col, "\n")
    state2 <- out$state
    exprs2 <- state2$expression
    samps2 <- state2$samples
    f.do_stuff(state2, config2)
  }
}
