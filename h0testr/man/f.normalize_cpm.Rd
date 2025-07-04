% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/normalize.R
\name{f.normalize_cpm}
\alias{f.normalize_cpm}
\title{Normalize expression using total expression in each sample.}
\usage{
f.normalize_cpm(state, config, multiplier = 1e+06)
}
\arguments{
\item{state}{List with elements formatted like the list returned by \code{f.read_data()}:
\tabular{ll}{
  \code{expression} \cr \tab Numeric matrix with non-negative expression values. \cr
  \code{features}   \cr \tab A data.frame with feature meta-data for rows of expression. \cr
  \code{samples}    \cr \tab A data.frame with observation meta-data for columns of expression. \cr
}}

\item{config}{List with configuration values. Does not use any keys, so can pass empty list.}

\item{multiplier}{Numeric greater than zero used to scale returned values 
after dividing by total counts in observation. For example, the default 
\code{multiplier=1e6} yields normalized expression as CPM (counts per million).}
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
Normalize expression using variant of CPM normalization which excludes missing values.
}
\details{
Inter-observation normalization, based on total expression
    in each observation. Makes total expression (excluding missing values)
    equal in each observation.
}
\examples{
set.seed(101)
exprs <- h0testr::f.sim1(n_obs=6, n_feats=12)$mat
feats <- data.frame(feature_id=rownames(exprs))
samps <- data.frame(observation_id=colnames(exprs))
state <- list(expression=exprs, features=feats, samples=samps)
config <- list()
state2 <- h0testr::f.normalize_cpm(state, config)
print(state$expression)
print(state2$expression)
apply(state$expression, 2, function(v) sum(v, na.rm=TRUE))
apply(state2$expression, 2, function(v) sum(v, na.rm=TRUE))
}
