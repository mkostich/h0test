#' Hypothesis testing using limma voom.
#' @description
#' `f.test_voom` tests for differential expression using the `limma::voom()` 
#'   function.
#' @details
#' Tests for differential expression using the `limma::voom()`. Model is fit 
#'   to `frm` and an F-test is performed for whether the effect of `test_term` 
#'   on `exprs` is zero. 
#' @param state List with elements formatted like the list returned by `f.read_data()`:
#'   \itemize{
#'     \item{expression}{Numeric matrix with feature rows and observation columns.}
#'     \item{features}{A data.frame with feature rows corresponding to columns of expression.}
#'     \item{samples}{A data.frame with observation rows corresponding to columns of expression.}
#'   }
#' @param config List with configuration values.
#' @param normalize.method Character in 
#'   `c("TMM", "TMMwsp", "RLE", "upperquartile", "none")`.
#' @return A data.frame containing results of test.
#' @examples
#' state <- list(expression=exprs, features=feats, samples=samps)
#' config <- list(frm=~age+sex+age:sex, test_term="age:sex")
#' tbl <- f.test_voom(state, config)
#' head(tbl)

f.test_voom <- function(state, config, normalize.method="none") {

  if(!is.matrix(state$expression)) {
    f.err("f.test_voom: !is.matrix(state$expression)", config=config)
  }
  
  exprs <- state$expression
  i <- apply(exprs, 1, function(v) any(is.na(v)))
  exprs <- exprs[!i, ]
  
  design <- model.matrix(config$frm, data=state$samples)
  obj <- limma::voom(exprs, design, plot=F, normalize.method=normalize.method)
  fit <- limma::lmFit(obj, design)
  fit <- limma::eBayes(fit, trend=F)
  
  frm0 <- as.formula(paste("~0 + ", config$test_term))
  cols <- colnames(model.matrix(frm0, data=state$samples))
  cols <- cols[cols %in% colnames(design)]
  i <- colnames(design) %in% cols
  tbl <- limma::topTable(fit, coef=which(i), number=Inf)
  
  f.msg("tested", nrow(exprs), "genes", config=config)
  f.msg("found", sum(tbl$adj.P.Val < 0.05, na.rm=T), "hits", config=config)
  
  return(tbl)
}

#' Hypothesis testing using limma trend.
#' @description
#' `f.test_trend` tests for differential expression using the `limma` 
#'   package trend approach.
#' @details
#' Tests for differential expression using the `limma` pakcage by running 
#'   `lmFit()`, then `eBayes()`, then `topTable()`. Model is fit 
#'   to `frm` and an F-test is performed for whether the effect of `test_term` 
#'   on `exprs` is zero. 
#' @param state List with elements formatted like the list returned by `f.read_data()`:
#'   \itemize{
#'     \item{expression}{Numeric matrix with feature rows and observation columns.}
#'     \item{features}{A data.frame with feature rows corresponding to columns of expression.}
#'     \item{samples}{A data.frame with observation rows corresponding to columns of expression.}
#'   }
#' @param config List with configuration values.
#' @param samps A data.frame with rows corresponding to columns of `exprs`, 
#'   and columns containing observation annotations.
#' @param frm Formula composed of terms containing variables corresponding
#'   to columns found in `samps`.
#' @param test_term Character term in `frm` on which hypothesis testing is to 
#'   be performed.
#' @return A data.frame containing results of test.
#' @examples
#' state <- list(expression=exprs, features=feats, samples=samps)
#' config <- list(frm=~age+sex+age:sex, test_term="age:sex")
#' tbl <- f.test_trend(state, config)
#' head(tbl)

f.test_trend <- function(state, config) {

  if(!is.matrix(state$expression)) {
    f.err("f.test_trend: !is.matrix(state$expression)")
  }
  
  design <- model.matrix(config$frm, data=state$samples)
  fit <- limma::lmFit(state$expression, design)
  fit <- limma::eBayes(fit, trend=T)
  
  frm0 <- as.formula(paste("~0 + ", config$test_term))
  cols <- colnames(model.matrix(frm0, data=state$samples))
  cols <- cols[cols %in% colnames(design)]
  i <- colnames(design) %in% cols
  tbl <- limma::topTable(fit, coef=which(i), number=Inf)
  
  f.msg("tested", nrow(state$expression), "genes")
  f.msg("found", sum(tbl$adj.P.Val < 0.05, na.rm=T), "hits")
  
  return(tbl)
}

