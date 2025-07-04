% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/initialize.R
\name{f.add_filter_stats}
\alias{f.add_filter_stats}
\title{Add filter statistics}
\usage{
f.add_filter_stats(state, config)
}
\arguments{
\item{state}{List with elements formatted like the list returned by \code{f.read_data()}:
\tabular{ll}{
  \code{expression} \cr \tab Numeric matrix with non-negative expression values. \cr
  \code{features}   \cr \tab A data.frame with feature meta-data for rows of expression. \cr
  \code{samples}    \cr \tab A data.frame with observation meta-data for columns of expression. \cr
}}

\item{config}{List with configuration values. Requires the following keys:
\tabular{ll}{
  \code{feat_id_col}  \cr \tab Name of column (character) in \code{feature_file_in} that corresponds to rows of \code{data_file_in}. \cr
  \code{obs_id_col}   \cr \tab Name of column (character) in \code{sample_file_in} that corresponds to columns of \code{expression}. \cr
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
Adds filtering-related statistics to \code{state$features}, 
    and \code{state$samples}.
}
\details{
Wrapper for \code{f.samples_per_feature()}, \code{f.feature_median_expression()}, 
    \code{f.features_per_sample()}. Also reports quantiles of distributions.
}
\examples{
exprs <- h0testr::f.sim1(n_obs=6, n_feats=8)$mat
feats <- data.frame(feature_id=rownames(exprs))
samps <- data.frame(observation_id=colnames(exprs))
state <- list(expression=exprs, features=feats, samples=samps)
config <- h0testr::f.new_config()
state <- h0testr::f.add_filter_stats(state, config)
print(state)
}
