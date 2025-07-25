% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/normalize.R
\name{f.normalize_qquantile}
\alias{f.normalize_qquantile}
\title{Old-school 'quantile' normalization.}
\usage{
f.normalize_qquantile(state, config)
}
\arguments{
\item{state}{List with elements formatted like the list returned by \code{f.read_data()}:
\tabular{ll}{
  \code{expression} \cr \tab Numeric matrix with non-negative expression values. \cr
  \code{features}   \cr \tab A data.frame with feature meta-data for rows of expression. \cr
  \code{samples}    \cr \tab A data.frame with observation meta-data for columns of expression. \cr
}}

\item{config}{List with configuration values. Does not use any keys, so can pass empty list.}
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
Normalize expression values using the 'quantile normalization' algorithm.
}
\details{
Inter-observation normalization resulting in nearly identical
    signal distributions across all samples, so all quantiles in \code{0:1} match
    across all samples. Calls \code{limma::normalizeQuantiles()} under the hood.
}
\examples{
set.seed(101)
exprs <- h0testr::f.sim1(n_obs=6, n_feats=12, mnar_c0=-Inf)$mat
feats <- data.frame(feature_id=rownames(exprs))
samps <- data.frame(observation_id=colnames(exprs))
state <- list(expression=exprs, features=feats, samples=samps)
config <- list()
state2 <- h0testr::f.normalize_qquantile(state, config)
print(state$expression)
print(state2$expression)

summary(apply(state$expression, 2, sd, na.rm=TRUE))
summary(apply(state2$expression, 2, sd, na.rm=TRUE))

## afterwards, all quantiles line up:
apply(state$expression, 2, quantile, probs=c(0.1, 0.25, 0.5, 0.75, 0.9), na.rm=TRUE)
apply(state2$expression, 2, quantile, probs=c(0.1, 0.25, 0.5, 0.75, 0.9), na.rm=TRUE)
}
