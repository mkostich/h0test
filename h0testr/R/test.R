#' Hypothesis testing using \code{limma::voom}
#' @description
#'   Tests for differential expression using the \code{limma::voom()} function.
#' @details
#'   Tests for differential expression using the \code{limma::voom()}. 
#'     Model is fit to \code{config$frm} and an F-test is performed for 
#'     whether the effect of \code{config$test_term} on 
#'     \code{state$expression} is zero. 
#' @param state List with elements formatted like the list returned by \code{f.read_data()}:
#'   \tabular{ll}{
#'     \code{expression} \cr \tab Numeric matrix with non-negative expression values. \cr
#'     \code{features}   \cr \tab A data.frame with feature meta-data for rows of expression. \cr
#'     \code{samples}    \cr \tab A data.frame with observation meta-data for columns of expression. \cr
#'   } 
#' @param config List with configuration values.
#' @param normalize.method Character in 
#'   \code{c("TMM", "TMMwsp", "RLE", "upperquartile", "none")}.
#' @return A data.frame containing results of test.
#' @examples
#'   \dontrun{
#'     state <- list(expression=exprs, features=feats, samples=samps)
#'     config <- list(frm=~age+sex+age:sex, test_term="age:sex")
#'     tbl <- f.test_voom(state, config)
#'     head(tbl)
#'   }

f.test_voom <- function(state, config, normalize.method="none") {
  
  if(!is.matrix(state$expression)) {
    f.err("f.test_voom: !is.matrix(state$expression)", config=config)
  }
  
  exprs <- state$expression
  i <- apply(exprs, 1, function(v) any(is.na(v)))
  exprs <- exprs[!i, ]
  
  design <- stats::model.matrix(config$frm, data=state$samples)
  obj <- limma::voom(exprs, design, plot=F, normalize.method=normalize.method)
  fit <- limma::lmFit(obj, design)
  fit <- limma::eBayes(fit, trend=F)
  
  frm0 <- stats::as.formula(paste("~0 + ", config$test_term))
  cols <- colnames(stats::model.matrix(frm0, data=state$samples))
  cols <- cols[cols %in% colnames(design)]
  i <- colnames(design) %in% cols
  tbl <- limma::topTable(fit, coef=which(i), number=Inf)
  
  f.msg("tested", nrow(exprs), "genes", config=config)
  f.msg("found", sum(tbl$adj.P.Val < 0.05, na.rm=T), "hits", config=config)
  
  return(tbl)
}

#' Hypothesis testing using \code{limma} trend.
#' @description
#'   Test for differential expression using the \code{limma} 
#'   package trend approach.
#' @details
#'   Tests for differential expression using the \code{limma} pakcage by running 
#'     \code{lmFit()}, then \code{eBayes()}, then \code{topTable()}. Model is fit 
#'     to \code{config$frm} and an F-test is performed for whether the effect of 
#'     \code{config$test_term} on \code{state$expression} is zero. 
#' @param state List with elements formatted like the list returned by \code{f.read_data()}:
#'   \tabular{ll}{
#'     \code{expression} \cr \tab Numeric matrix with non-negative expression values. \cr
#'     \code{features}   \cr \tab A data.frame with feature meta-data for rows of expression. \cr
#'     \code{samples}    \cr \tab A data.frame with observation meta-data for columns of expression. \cr
#'   } 
#' @param config List with configuration values.
#' @return A \code{data.frame} containing results of test.
#' @examples
#'   \dontrun{
#'     state <- list(expression=exprs, features=feats, samples=samps)
#'     config <- list(frm=~age+sex+age:sex, test_term="age:sex")
#'     tbl <- f.test_trend(state, config)
#'     head(tbl)
#'   }

f.test_trend <- function(state, config) {

  if(!is.matrix(state$expression)) {
    f.err("f.test_trend: !is.matrix(state$expression)", config=config)
  }
  
  design <- stats::model.matrix(config$frm, data=state$samples)
  fit <- limma::lmFit(state$expression, design)
  fit <- limma::eBayes(fit, trend=T)
  
  frm0 <- stats::as.formula(paste("~0 + ", config$test_term))
  cols <- colnames(stats::model.matrix(frm0, data=state$samples))
  cols <- cols[cols %in% colnames(design)]
  i <- colnames(design) %in% cols
  tbl <- limma::topTable(fit, coef=which(i), number=Inf)
  
  f.msg("tested", nrow(state$expression), "genes", config=config)
  f.msg("found", sum(tbl$adj.P.Val < 0.05, na.rm=T), "hits", config=config)
  
  return(tbl)
}

#' Hypothesis testing
#' @description
#'   Test hypotheses using .
#' @details
#'   Tests for differential expression using method specified in config. 
#' @param state List with elements formatted like the list returned by `f.read_data()`:
#'   \tabular{ll}{
#'     \code{expression} \cr \tab Numeric matrix with non-negative expression values. \cr
#'     \code{features}   \cr \tab A data.frame with feature meta-data for rows of expression. \cr
#'     \code{samples}    \cr \tab A data.frame with observation meta-data for columns of expression. \cr
#'   } 
#' @param config List with configuration values.
#' @return A data.frame containing results of test.
#' @examples
#'   \dontrun{
#'     state <- list(expression=exprs, features=feats, samples=samps)
#'     config <- list(frm=~age+sex+age:sex, test_term="age:sex", test_method="trend")
#'     tbl <- f.test(state, config)
#'     head(tbl)
#'   }

f.test <- function(state, config) {

  if(config$test_method %in% "voom") {
    tbl <- f.test_voom(state, config)
  } else if(config$test_method %in% "trend") {
    tbl <- f.test_trend(state, config)
  } else if(config$test_method %in% "none") {
    f.msg("skipping testing: TEST_METHOD %in% 'none'", config=config)
    return(NULL)
  } else f.err("unexpected TEST_METHOD:", config$test_method, config=config)
  
  rownames(state$features) <- state$features[[config$feat_id_col]]
  if(!all(rownames(tbl) %in% rownames(state$features))) {
    f.err("!all(rownames(tbl) %in% rownames(state$features))", config=config)
  }
  
  tbl <- cbind(state$features[rownames(tbl), ], tbl)
  rownames(tbl) <- NULL
  
  file_out <- paste0(config$dir_out, "/", length(config$run_order) + 3, config$result_mid_out, config$suffix_out)
  f.log("writing results to", file_out, config=config)
  f.save_tsv(tbl, file_out)
  
  return(tbl)
}

#' Run a basic workflow
#' @description
#'   Run a basic workflow according to: 
#'     \code{config$run_order}.
#' @details
#'   Run a basic workflow: 
#'     \code{f.load_data -> config$run_order -> f.test}, where 
#'       \code{config$run_order} is vector of functions which are run in
#'       the specified order.
#' @param config List with configuration values.
#' @return A list with the following elements:
#'   \tabular{ll}{
#'     \code{state}  \cr \tab A list with elements \code{$expression}, \code{$features}, and \code{$samples}. \cr
#'     \code{config} \cr \tab A list with configuration settings. \cr
#'     \code{tbl}    \cr \tab A data.frame containing results of test.. \cr
#'   }
#' @examples
#'   \dontrun{
#'     config <- list(frm=~age+sex+age:sex, test_term="age:sex", test_method="trend")
#'     config$run_order <- c("normalize", "combine_reps", "filter", "impute")
#'     output <- f.run(config)
#'     tbl <- output$tbl
#'     head(tbl)
#'   }

f.run <- function(config) {

  f.log_block("starting f.load_data", config=config)
  out <- f.load_data(config)
  
  for(f in config$run_order) {
    f <- paste0("f.", f)
    f.log_block("starting", f, config=config)
    f <- get(f)
    out <- f(out$state, out$config)
  }
  
  f.log_block("starting f.test", config=out$config)
  tbl <- f.test(out$state, out$config)
  
  return(list(state=out$state, config=out$config, tbl=tbl))
}

