% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/normalize.R
\name{f.normalize_edger}
\alias{f.normalize_edger}
\title{Normalize expression using \code{edgeR}}
\usage{
f.normalize_edger(state, config, method = NULL, p = NULL)
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
  \code{norm_method}     \cr \tab Normalization method in \code{c("TMM", "TMMwsp", "RLE", "upperquartile")}. \cr
  \code{norm_quantile}   \cr \tab Quantile to use for \code{norm_method \%in\% c("quantile", "upperquartile")}. \cr
}}

\item{method}{Character in set
\code{c("TMM", "TMMwsp", "RLE", "upperquartile", "none")}. Can be set with
\code{config$norm_method}.}

\item{p}{Numeric in closed interval \code{[0, 1]} specifying quantile to use for method 
\code{upperquartile}. Can be set with \code{config$norm_quantile}.}
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
Normalize expression using functionality from \code{edgeR} package.
}
\details{
Inter-observation normalization. \code{edgeR::calcNormFactors()} 
  called under the hood. Returned values on a counts-per-million scale.
}
\examples{
set.seed(101)
exprs <- h0testr::f.sim1(n_obs=6, n_feats=12)$mat
feats <- data.frame(feature_id=rownames(exprs))
samps <- data.frame(observation_id=colnames(exprs))
state <- list(expression=exprs, features=feats, samples=samps)

config <- list(norm_method="upperquartile", norm_quantile=0.75)
state2 <- h0testr::f.normalize_edger(state, config)
print(state$expression)
print(state2$expression)

## norm_quantile ignored for norm_method="TMM":
config <- list(norm_method="TMM", norm_quantile=0)
state2 <- h0testr::f.normalize_edger(state, config)
print(state$expression)
print(state2$expression)
}
