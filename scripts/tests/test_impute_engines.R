## Tests for the imputation engines that wrap an external package, and for the
##   h0testr::impute() dispatch onto them. Each engine is skipped, not failed,
##   when the package it wraps is not installed.

usage <- function(msg=NULL) {
  if(!is.null(msg)) cat("ERROR:", msg, "\n\n", file=stderr())
  cat(
    "Test the imputation engines that wrap an external package:",
    "h0testr::impute_knn(), impute_min_det(), impute_min_prob(), impute_qrilc(),",
    "impute_pca(), impute_missforest(), impute_rf(), impute_glmnet() and",
    "impute_glm_binom(). Checks that each preserves the shape, the dimnames, the",
    "measured values and the metadata tables, resolves every missing value to a",
    "finite value that could be a measurement, and that h0testr::impute()",
    "dispatches onto each with the same result under the same seed. Engines whose",
    "package is not installed are skipped rather than failed.",
    "",
    "Usage: Rscript test_impute_engines.R <r_dir>",
    "",
    "Required positional arguments:",
    "  <r_dir>  Path to the h0testr package R/ source directory; all .R files",
    "             there are sourced (the installed package is not used).",
    "",
    "Output: one PASS/FAIL line per assertion to stdout, then a count of",
    "  passes and failures. Log output, including messages from expected",
    "  errors, is written to a temporary file whose path is reported at the end.",
    "",
    "Exit codes: 0 all assertions passed; 1 one or more failed; 2 usage error.",
    "",
    "Examples:",
    "  Rscript test_impute_engines.R C:/path/to/h0test/h0testr/R",
    "  Rscript test_impute_engines.R ../../h0test/h0testr/R",
    "  Rscript test_impute_engines.R C:/path/to/h0testr/R > t.out 2>&1",
    sep="\n", file=stderr()
  )
  quit(status=2)
}

args <- commandArgs(trailingOnly=TRUE)
if(length(args) != 1) usage("wrong number of arguments")

r_dir <- args[1]
if(!dir.exists(r_dir)) usage(paste("r_dir not a directory:", r_dir))

for(f in list.files(r_dir, pattern="[.]R$", full.names=TRUE)) source(f)

###############################################################################

log_file <- tempfile(fileext=".log")
invisible(file.create(log_file))
n_pass <- 0
n_fail <- 0
n_skip <- 0
t0 <- Sys.time()

report <- function(ok, msg) {
  if(isTRUE(ok)) {
    n_pass <<- n_pass + 1
    cat("PASS:", msg, "\n")
  } else {
    n_fail <<- n_fail + 1
    cat("FAIL:", msg, "\n")
  }
  utils::flush.console()
}

section <- function(msg) {
  cat("\n##", msg, "; elapsed:",
    round(as.numeric(difftime(Sys.time(), t0, units="secs"))), "s\n")
  utils::flush.console()
}

skip <- function(label, pkg) {
  n_skip <<- n_skip + 1
  cat("SKIP:", label, "; package not installed:", pkg, "\n")
  utils::flush.console()
}

threw <- function(expr) inherits(try(expr, silent=TRUE), "try-error")

has <- function(pkg) requireNamespace(pkg, quietly=TRUE)

mk <- function(e) {
  list(expression=e, features=data.frame(feature_id=rownames(e)),
    samples=data.frame(observation_id=colnames(e)))
}

cfg0 <- list(feat_col="feature_id", obs_col="observation_id",
  sample_id_col="observation_id", log_file=log_file, save_state=FALSE,
  verbose=FALSE)
cfg_raw <- c(cfg0, list(is_log_transformed=FALSE))
cfg_log <- c(cfg0, list(is_log_transformed=TRUE))

## the whole-matrix fixture; every feature measured somewhere, so no engine is
##   asked for a median or a fit it cannot have:

set.seed(101)
e_raw <- sim1(n_obs=8, n_feats=60, mcar_p=0.05)$mat
e_raw <- e_raw[apply(e_raw, 1, function(v) sum(!is.na(v))) >= 3, , drop=FALSE]
e_raw <- e_raw[seq_len(min(40, nrow(e_raw))), , drop=FALSE]
e_log <- log2(e_raw + 1)

s_raw <- mk(e_raw)
s_log <- mk(e_log)

## a smaller fixture for the per-feature model fitters, which are slow:

set.seed(202)
e_sm <- sim1(n_obs=10, n_feats=40, mcar_p=0.05)$mat
e_sm <- e_sm[apply(e_sm, 1, function(v) sum(!is.na(v))) >= 4, , drop=FALSE]
e_sm <- e_sm[seq_len(min(20, nrow(e_sm))), , drop=FALSE]
s_sm <- mk(e_sm)

###############################################################################
section("fixtures")

report(any(is.na(e_raw)) && any(!is.na(e_raw)),
  "the fixture carries both measured and missing values")
report(all(e_raw[!is.na(e_raw)] > 0),
  "and every measured value is above zero, so the raw scale is usable")
report(all(apply(e_raw, 2, function(v) any(!is.na(v)))),
  "and no observation is missing everywhere")
report(any(is.na(e_sm)) && all(e_sm[!is.na(e_sm)] > 0),
  "the small fixture carries missing values among positive measurements")

###############################################################################

## the contract every engine owes its caller:

## tol is the relative deviation allowed in a measured value: zero for an engine
##   that touches only the missing cells, and a floating point margin for one
##   that round-trips the whole matrix through log2(x + 1).

chk <- function(label, r, s_in, raw, tol=0) {
  e_in <- s_in$expression
  i <- is.na(e_in)
  m <- r$expression
  report(is.matrix(m), paste(label, "returns a state carrying an expression matrix"))
  if(!is.matrix(m)) return(invisible(NULL))
  report(identical(dim(m), dim(e_in)), paste(label, "keeps the shape of its input"))
  report(identical(dimnames(m), dimnames(e_in)),
    paste(label, "keeps the feature and observation names"))
  report(!any(is.na(m)), paste(label, "fills in every missing value"))
  dev <- abs(m[!i] - e_in[!i]) / pmax(abs(e_in[!i]), 1)
  report(max(dev) <= tol, paste(label, "leaves every measured value alone"))
  report(all(is.finite(m)), paste(label, "imputes no infinite value"))
  if(raw) {
    report(all(m > 0), paste(label, "imputes nothing at or below zero on the raw scale"))
  }
  report(identical(r$features, s_in$features) && identical(r$samples, s_in$samples),
    paste(label, "passes the feature and observation tables through untouched"))
}

## one seeded direct call, then the same call through impute(): identical
##   results say the dispatch forwards everything and that the engine is
##   reproducible under a seed.

engine <- function(label, method, s_in, cfg, raw, seed, direct, tol=0) {
  set.seed(seed)
  r <- try(direct(), silent=TRUE)
  report(!inherits(r, "try-error"), paste(label, "runs"))
  if(inherits(r, "try-error")) return(invisible(NULL))
  chk(label, r, s_in, raw, tol=tol)
  set.seed(seed)
  o <- try(suppressMessages(impute(s_in, cfg, method=method)), silent=TRUE)
  report(!inherits(o, "try-error"), paste0("impute(method='", method, "') runs"))
  if(inherits(o, "try-error")) return(invisible(r))
  report(identical(o$state$expression, r$expression),
    paste0("impute(method='", method, "') reproduces the direct call under one seed"))
  invisible(r)
}

refuses_df <- function(label, fn, cfg) {
  bad <- mk(e_raw)
  bad$expression <- as.data.frame(e_raw)
  report(threw(fn(bad, cfg)),
    paste(label, "refuses a state whose expression is not a matrix"))
}

###############################################################################
section("impute_knn(): impute::impute.knn()")

if(!has("impute")) {
  skip("impute_knn", "impute")
} else {
  r <- engine("impute_knn", "knn", s_raw, cfg_raw, raw=TRUE, seed=11,
    function() impute_knn(s_raw, cfg_raw))
  refuses_df("impute_knn", impute_knn, cfg_raw)

  ## k is clamped to sqrt(n_features), since a neighborhood cannot be larger
  ##   than the pool it is drawn from:
  r_big <- impute_knn(s_raw, cfg_raw, k=nrow(e_raw))
  r_max <- impute_knn(s_raw, cfg_raw, k=round(sqrt(nrow(e_raw))))
  report(identical(r_big$expression, r_max$expression),
    "impute_knn() clamps k to the square root of the feature count")
}

###############################################################################
section("impute_min_det(): imputeLCMD::impute.MinDet()")

if(!has("imputeLCMD")) {
  skip("impute_min_det", "imputeLCMD")
} else {
  r <- engine("impute_min_det", "min_det", s_raw, cfg_raw, raw=TRUE, seed=13,
    function() impute_min_det(s_raw, cfg_raw))
  refuses_df("impute_min_det", impute_min_det, cfg_raw)

  if(!is.null(r)) {
    i_na <- is.na(e_raw)
    fill <- lapply(seq_len(ncol(e_raw)),
      function(j) unique(r$expression[i_na[, j], j]))
    fill <- fill[lengths(fill) > 0]
    report(all(lengths(fill) == 1),
      "impute_min_det() substitutes one deterministic value per observation")
    j_na <- which(colSums(i_na) > 0)
    q01 <- mapply(function(j) stats::quantile(e_raw[, j], 0.01, na.rm=TRUE), j_na)
    report(isTRUE(all.equal(unname(unlist(fill)), unname(q01))),
      "and that value is the 1% quantile of the measurements in its observation")
    report(all(unlist(fill) <
      mapply(function(j) stats::median(e_raw[, j], na.rm=TRUE), j_na)),
      "which sits in the low tail, as a limit of detection should")
  }

  ## a higher quantile of the observed values is a higher floor:
  q_lo <- impute_min_det(s_raw, cfg_raw, impute_quantile=0.01)$expression
  q_hi <- impute_min_det(s_raw, cfg_raw, impute_quantile=0.2)$expression
  i_na <- is.na(e_raw)
  report(all(q_hi[i_na] >= q_lo[i_na]) && any(q_hi[i_na] > q_lo[i_na]),
    "impute_min_det() raises the imputed values with impute_quantile")
}

###############################################################################
section("impute_min_prob(): imputeLCMD::impute.MinProb()")

if(!has("imputeLCMD")) {
  skip("impute_min_prob", "imputeLCMD")
} else {
  engine("impute_min_prob", "min_prob", s_raw, cfg_raw, raw=TRUE, seed=17,
    function() impute_min_prob(s_raw, cfg_raw), tol=1e-10)
  engine("impute_min_prob (log scale)", "min_prob", s_log, cfg_log, raw=FALSE,
    seed=17, function() impute_min_prob(s_log, cfg_log))
  refuses_df("impute_min_prob", impute_min_prob, cfg_raw)

  i_na <- is.na(e_log)
  set.seed(19)
  s1 <- stats::sd(impute_min_prob(s_log, cfg_log, scale.=0.1)$expression[i_na])
  set.seed(19)
  s3 <- stats::sd(impute_min_prob(s_log, cfg_log, scale.=1.5)$expression[i_na])
  report(s3 > s1, "impute_min_prob() widens the draw with scale.")

  report(threw(impute_min_prob(s_log, cfg_log, is_log_transformed=FALSE)),
    "impute_min_prob() refuses an is_log_transformed disagreeing with the config")
}

###############################################################################
section("impute_qrilc(): imputeLCMD::impute.QRILC()")

if(!has("imputeLCMD")) {
  skip("impute_qrilc", "imputeLCMD")
} else {
  r <- engine("impute_qrilc", "qrilc", s_log, cfg_log, raw=FALSE, seed=23,
    function() impute_qrilc(s_log, cfg_log))
  refuses_df("impute_qrilc", impute_qrilc, cfg_log)

  if(!is.null(r)) {
    i_na <- is.na(e_log)
    report(stats::median(r$expression[i_na]) <
      stats::median(e_log[!i_na], na.rm=TRUE),
      "impute_qrilc() imputes below the measured values, as left censoring")
  }

  set.seed(29)
  s1 <- stats::sd(impute_qrilc(s_log, cfg_log, scale.=0.1)$expression[is.na(e_log)])
  set.seed(29)
  s3 <- stats::sd(impute_qrilc(s_log, cfg_log, scale.=1.5)$expression[is.na(e_log)])
  report(s3 > s1, "impute_qrilc() widens the draw with scale.")
}

###############################################################################
section("impute_pca(): pcaMethods::pca()")

if(!has("pcaMethods")) {
  skip("impute_pca", "pcaMethods")
} else {
  engine("impute_pca (bpca)", "bpca", s_raw, cfg_raw, raw=TRUE, seed=31,
    function() impute_pca(s_raw, cfg_raw, method="bpca"), tol=1e-10)
  engine("impute_pca (ppca)", "ppca", s_raw, cfg_raw, raw=TRUE, seed=37,
    function() impute_pca(s_raw, cfg_raw, method="ppca"), tol=1e-10)
  engine("impute_pca (svdImpute)", "svdImpute", s_log, cfg_log, raw=FALSE,
    seed=41, function() impute_pca(s_log, cfg_log, method="svdImpute"))
  refuses_df("impute_pca", impute_pca, cfg_raw)

  report(threw(impute_pca(s_raw, cfg_raw, method="not_a_method")),
    "impute_pca() refuses a method outside bpca, ppca and svdImpute")

  ## n_pcs cannot exceed what the matrix can support, and is clamped:
  n_max <- round(sqrt(nrow(e_log)))
  set.seed(43)
  a <- impute_pca(s_log, cfg_log, method="svdImpute", n_pcs=nrow(e_log))
  set.seed(43)
  b <- impute_pca(s_log, cfg_log, method="svdImpute", n_pcs=n_max)
  report(identical(a$expression, b$expression),
    "impute_pca() clamps n_pcs to the square root of the feature count")
}

###############################################################################
section("impute_missforest(): missForest::missForest()")

if(!has("missForest")) {
  skip("impute_missforest", "missForest")
} else {
  engine("impute_missforest", "missforest", s_sm, cfg_raw, raw=TRUE, seed=47,
    function() impute_missforest(s_sm, cfg_raw))
  refuses_df("impute_missforest", impute_missforest, cfg_raw)
}

###############################################################################
section("impute_rf(): randomForest::randomForest()")

if(!has("randomForest")) {
  skip("impute_rf", "randomForest")
} else {
  engine("impute_rf", "rf", s_sm, cfg_raw, raw=TRUE, seed=53,
    function() impute_rf(s_sm, cfg_raw, verbose=FALSE)$state)
  refuses_df("impute_rf", function(s, cfg) impute_rf(s, cfg, verbose=FALSE), cfg_raw)

  set.seed(59)
  out <- impute_rf(s_sm, cfg_raw, aug_steps=0, verbose=FALSE)
  report(is.list(out) && all(c("state", "log") %in% names(out)),
    "impute_rf() returns the state alongside a log of its per-feature fits")
  report(is.null(out$log) || nrow(out$log) ==
    sum(apply(e_sm, 1, function(v) any(is.na(v)))),
    "with one log row per feature that had a missing value")

  report(threw(impute_rf(s_sm, cfg_raw, aug_steps=-1, verbose=FALSE)),
    "impute_rf() refuses a negative aug_steps")
}

###############################################################################
section("impute_glmnet(): glmnet::cv.glmnet()")

if(!has("glmnet")) {
  skip("impute_glmnet", "glmnet")
} else {
  engine("impute_glmnet", "glmnet", s_sm, cfg_raw, raw=TRUE, seed=61,
    function() impute_glmnet(s_sm, cfg_raw, verbose=FALSE)$state)
  refuses_df("impute_glmnet",
    function(s, cfg) impute_glmnet(s, cfg, verbose=FALSE), cfg_raw)

  report(threw(impute_glmnet(s_sm, cfg_raw, alpha=1.5, verbose=FALSE)),
    "impute_glmnet() refuses an alpha outside zero and one")
  report(threw(impute_glmnet(s_sm, cfg_raw, aug_steps=-1, verbose=FALSE)),
    "impute_glmnet() refuses a negative aug_steps")
}

###############################################################################
section("impute_glm_binom(): missingness against intensity")

cfg_gb <- cfg_raw
cfg_gb$impute_n_pts <- 1e5

engine("impute_glm_binom", "glm_binom", s_raw, cfg_gb, raw=TRUE, seed=67,
  function() impute_glm_binom(s_raw, cfg_gb, n_pts=1e5))
refuses_df("impute_glm_binom",
  function(s, cfg) impute_glm_binom(s, cfg, n_pts=1e3), cfg_raw)

report(threw(impute_glm_binom(s_raw, cfg_raw, n_pts=0)),
  "impute_glm_binom() refuses a prediction grid with no points")
report(threw(impute_glm_binom(s_raw, cfg_raw, n_pts=1e3, min_fit_pts=nrow(e_raw) + 1)),
  "impute_glm_binom() refuses a fit with fewer points than min_fit_pts")

###############################################################################
section("every engine is reachable through impute()")

report(all(c("knn", "min_det", "min_prob", "qrilc", "bpca", "ppca", "svdImpute",
  "missforest", "rf", "glmnet", "glm_binom") %in% impute_methods()),
  "impute_methods() names every engine tested here")

report(threw(suppressMessages(impute(s_raw, cfg_raw, method="not_a_method"))),
  "impute() refuses a method it does not know")

###############################################################################

cat("\n## log file:", log_file, "\n")
cat("\n## passes:", n_pass, "; failures:", n_fail, "; skipped engines:", n_skip,
  "; elapsed:", round(as.numeric(difftime(Sys.time(), t0, units="secs"))), "s\n")

quit(save="no", status=if(n_fail > 0) 1 else 0)
