#' Hypothesis testing using limma voom.
#' @description
#' `f.test_voom` tests for differential expression using the `limma::voom()` 
#'   function.
#' @details
#' Tests for differential expression using the `limma::voom()`. Model is fit 
#'   to `frm` and an F-test is performed for whether the effect of `test_term` 
#'   on `exprs` is zero. 
#' @param exprs Numeric matrix with feature rows and observation columns.
#' @param samps A data.frame with rows corresponding to columns of `exprs`, 
#'   and columns containing observation annotations.
#' @param frm Formula composed of terms containing variables corresponding
#'   to columns found in `samps`.
#' @param test_term Character term in `frm` on which hypothesis testing is to 
#'   be performed.
#' @param normalize.method Character in 
#'   `c("TMM", "TMMwsp", "RLE", "upperquartile", "none")`.
#' @return A data.frame containing results of test.
#' @examples
#' tbl <- f.test_voom(exprs, samps, frm=~age+sex+age:sex, test_term="age:sex")
#' head(tbl)

f.test_voom <- function(exprs, samps, frm, test_term, normalize.method="none") {
  if(!is.matrix(exprs)) f.err("f.test_voom: !is.matrix(exprs)")
  i <- apply(exprs, 1, function(v) any(is.na(v)))
  exprs <- exprs[!i, ]
  design <- model.matrix(frm, data=samps)
  obj <- limma::voom(exprs, design, plot=F, normalize.method=normalize.method)
  fit <- limma::lmFit(obj, design)
  fit <- limma::eBayes(fit, trend=F)
  cols <- colnames(model.matrix(as.formula(paste("~0 + ", test_term)), data=samps))
  cols <- cols[cols %in% colnames(design)]
  i <- colnames(design) %in% cols
  tbl <- limma::topTable(fit, coef=which(i), number=Inf)
  f.msg("tested", nrow(exprs), "genes")
  f.msg("found", sum(tbl$adj.P.Val < 0.05, na.rm=T), "hits")
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
#' @param exprs Numeric matrix with feature rows and observation columns.
#' @param samps A data.frame with rows corresponding to columns of `exprs`, 
#'   and columns containing observation annotations.
#' @param frm Formula composed of terms containing variables corresponding
#'   to columns found in `samps`.
#' @param test_term Character term in `frm` on which hypothesis testing is to 
#'   be performed.
#' @return A data.frame containing results of test.
#' @examples
#' tbl <- f.test_trend(exprs, samps, frm=~age+sex+age:sex, test_term="age:sex")
#' head(tbl)

f.test_trend <- function(exprs, samps, frm, test_term) {
  if(!is.matrix(exprs)) f.err("f.test_trend: !is.matrix(exprs)")
  design <- model.matrix(frm, data=samps)
  fit <- limma::lmFit(exprs, design)
  fit <- limma::eBayes(fit, trend=T)
  cols <- colnames(model.matrix(as.formula(paste("~0 + ", test_term)), data=samps))
  cols <- cols[cols %in% colnames(design)]
  i <- colnames(design) %in% cols
  tbl <- limma::topTable(fit, coef=which(i), number=Inf)
  f.msg("tested", nrow(exprs), "genes")
  f.msg("found", sum(tbl$adj.P.Val < 0.05, na.rm=T), "hits")
  return(tbl)
}

