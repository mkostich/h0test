% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/impute.R
\name{f.impute_glm_binom}
\alias{f.impute_glm_binom}
\title{Impute by drawing from `p(missing|intensity)` from binomial glm.}
\usage{
f.impute_glm_binom(state, config, n_pts = NULL, off = 1, f_mid = stats::median)
}
\arguments{
\item{state}{A list with elements like that returned by \code{f.read_data()}:
  \tabular{ll}{
    \code{expression} \cr \tab Numeric matrix with non-negative expression values. \cr
    \code{features}   \cr \tab A data.frame with feature meta-data for rows of expression. \cr
    \code{samples}    \cr \tab A data.frame with observation meta-data for columns of expression. \cr
  } 
It is assumed that state$expression has been previously \code{log2(x+1)} transformed.}

\item{config}{List with configuration values. Uses the following keys:
\tabular{ll}{
  \code{impute_n_pts}  \cr \tab Numeric greater than one. Determines granularity of imputation. Larger values lead to finer grain.
}}

\item{n_pts}{Numeric greater than one. Granularity of prediction 
grid. Larger values lead to less chance of duplicate imputed values.
Larger values require more compute time and memory.}

\item{off}{Numeric offset for calculating 
\code{p.missing = (n.missing + off) / (n.total + off)}.}

\item{f_mid}{Function to use for calculating central tendency of 
feature expression across samples.}
}
\value{
An updated `state` list with the following elements:
  \tabular{ll}{
    \code{expression} \cr \tab Numeric matrix with non-negative expression values. \cr
    \code{features}   \cr \tab A data.frame with feature meta-data for rows of expression. \cr
    \code{samples}    \cr \tab A data.frame with observation meta-data for columns of expression. \cr
  }
}
\description{
Impute missing values by randomly drawing from an estimated density 
    of \code{p(missing|intensity)}.
}
\details{
Impute random draws by drawing from the estimated density of 
  \code{p(missing|intensity)}. Estimates \code{p(missing|intensity)} using 
  a binomial glm with formula \code{cbind(n.missing, n.found) ~ intensity},
  where intensity is assumed to be previously log transformed. If you  
  want \code{0} to be considered missing, and have \code{0} in the data, do 
  something like \code{exprs[exprs \%in\% 0] <- NA} prior to imputing.
}
\examples{
set.seed(101)
exprs <- h0testr::f.sim1(n_obs=6, n_feats=12)$mat
exprs <- log2(exprs + 1)
feats <- data.frame(feature_id=rownames(exprs))
samps <- data.frame(observation_id=colnames(exprs))
state <- list(expression=exprs, features=feats, samples=samps)
config <- list(impute_n_pts=1e7)
state2 <- h0testr::f.impute_glm_binom(state, config)
summary(c(state$expression))   ## note number of NAs
summary(c(state2$expression))  ## note number of NAs
}
