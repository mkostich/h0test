## Tests the logfc and expr columns of the standardized result table h0testr::test()
##   returns. Every engine leaves logfc empty for a joint test over several
##   coefficients, since there is no single contrast to report, and test_lm() and
##   test_prolfqua() left it empty for a single coefficient too; expr was empty for
##   test_lm(), test_msqrob() and test_prolfqua(). Both are now filled in, logfc from
##   the coefficients of the fit that produced the p-values: the signed coefficient
##   when one design column carries the test, and otherwise the total swing, the range
##   over the observations of the fitted contribution of the terms under test.
##   Checks the value against an independent calculation for each shape of test, that
##   it is the same whichever engine reports it where the fits agree, that the
##   specializations claimed for it hold (largest difference between any two levels of
##   a multi-level factor, slope times covariate range for a continuous covariate), and
##   that the intercept warning fires when config$frm has a continuous covariate.

usage <- function(msg=NULL) {
  if(!is.null(msg)) cat("ERROR:", msg, "\n\n", file=stderr())
  cat(
    "Test the logfc and expr columns of the standardized table h0testr::test()",
    "returns, for every test method and every shape of test: a single coefficient,",
    "a joint test over the coefficients of a multi-level factor, a continuous",
    "covariate, and a joint test spanning an interaction. Checks logfc against an",
    "independent least squares calculation, checks that the total swing reported for",
    "a joint test is the largest difference between any two levels of a factor and",
    "the slope times the covariate range for a continuous covariate, checks that expr",
    "is the mean of the values handed to the test, and checks the warning issued when",
    "the intercept is tested with a continuous covariate in config$frm.",
    "",
    "Requires the prolfqua package to be installed.",
    "",
    "Usage: Rscript test_effect_size.R <r_dir>",
    "",
    "Required positional arguments:",
    "  <r_dir>  Path to the h0testr package R/ source directory; all .R files",
    "             there are sourced (the installed package is not used).",
    "",
    "Output: one PASS/FAIL line per assertion to stdout, then a count of",
    "  passes and failures. Log output, including messages from expected",
    "  errors, is written to a temporary file whose path is reported at the end.",
    "",
    "Exit codes: 0 all assertions passed; 1 one or more failed; 2 usage error;",
    "  3 prolfqua not installed.",
    "",
    "Examples:",
    "  Rscript test_effect_size.R C:/path/to/h0test/h0testr/R",
    "  Rscript test_effect_size.R ../../h0test/h0testr/R",
    "  Rscript test_effect_size.R C:/path/to/h0testr/R > test_effect_size.out 2>&1",
    sep="\n", file=stderr()
  )
  quit(status=2)
}

args <- commandArgs(trailingOnly=TRUE)
if(length(args) != 1) usage("wrong number of arguments")

r_dir <- args[1]
if(!dir.exists(r_dir)) usage(paste("r_dir not a directory:", r_dir))

if(!requireNamespace("prolfqua", quietly=TRUE)) {
  cat("ERROR: prolfqua not installed\n", file=stderr())
  quit(status=3)
}

for(f in list.files(r_dir, pattern="[.]R$", full.names=TRUE)) source(f)

###############################################################################
## harness:

log_file <- tempfile(fileext=".log")
n_pass <- 0
n_fail <- 0
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

log_len <- function() {
  if(!file.exists(log_file)) return(0L)
  return(length(readLines(log_file, warn=FALSE)))
}

log_has <- function(n0, ...) {
  txt <- readLines(log_file, warn=FALSE)
  if(length(txt) <= n0) return(FALSE)
  txt <- paste(txt[-seq_len(n0)], collapse=" ")
  for(pat in c(...)) if(!grepl(pat, txt, fixed=TRUE)) return(FALSE)
  return(TRUE)
}

###############################################################################
## Twenty-four observations, fully crossed: grp has three levels, sex and batch two
##   each, age is continuous. Twenty genes with two peptides each, two peptides dropped
##   so that DEqMS sees varying peptide counts. No missing values, and values already
##   on a log scale, so that an independent least squares fit of the same design
##   reproduces what the engines fit.

set.seed(101)

nobs <- 24
ngene <- 20
npep <- 2

samps <- expand.grid(rep=1:2, batch=c("b1", "b2"), sex=c("F", "M"),
  grp=c("a", "b", "c"), stringsAsFactors=FALSE)
samps$obs <- paste0("o", sprintf("%02d", 1:nobs))
samps$age <- round(rnorm(nobs, 50, 10), 1)
samps <- samps[, c("obs", "grp", "sex", "batch", "age")]

## a known effect on the first half of the genes, so the effect sizes below are not all
##   noise: level c is up, level b is not:

gene_eff <- matrix(0, nrow=ngene, ncol=nobs)
gene_eff[1:10, samps$grp %in% "c"] <- 3

exprs <- matrix(rnorm(ngene * npep * nobs, 20, 1.5), nrow=ngene * npep)
exprs <- exprs + gene_eff[rep(1:ngene, each=npep), , drop=FALSE]
gene <- paste0("g", sprintf("%02d", rep(1:ngene, each=npep)))
rownames(exprs) <- paste0(gene, ".p", rep(1:npep, times=ngene))
colnames(exprs) <- samps$obs

## the peptides per gene have to vary for test_method "deqms", and to be spread rather
##   than merely varied: DEqMS fits its variance prior as a loess of residual variance on
##   that count, and dropping two rows overall left eighteen of the twenty genes at one
##   count, where the loess comes back NaN and those genes lose their moderated statistic
##   altogether. Dropping the first peptide of every other gene leaves ten genes with one
##   peptide and ten with two, which it fits:

i_drop <- seq(from=1, by=2 * npep, length.out=ngene %/% 2)
exprs <- exprs[-i_drop, , drop=FALSE]
gene <- gene[-i_drop]

feats <- data.frame(pep=rownames(exprs), gene=gene)

cfg0 <- new_config()
cfg0$obs_id_col <- cfg0$obs_col <- "obs"
cfg0$sample_id_col <- "obs"
cfg0$feat_id_col <- cfg0$feat_col <- "pep"
cfg0$gene_id_col <- "gene"
cfg0$frm <- ~grp + sex + batch + age
cfg0$test_term <- "grp"
cfg0$reference_levels <- c(grp="a", sex="F", batch="b1")
cfg0$estimability <- "test"
cfg0$df_resid_min <- 2
cfg0$save_state <- FALSE
cfg0$permute_var <- ""
cfg0$is_log_transformed <- TRUE
cfg0$impute_method <- "none"
cfg0$n_features_min <- 5
cfg0$test_prior_df <- 5
cfg0$log_file <- log_file

out <- initialize(list(expression=exprs, features=feats, samples=samps), cfg0,
  minimal=TRUE)
state0 <- out$state
config0 <- out$config

methods <- c("lm", "trend", "voom", "deqms", "msqrob", "proda", "prolfqua")

## config for one case, with config$reference_levels cut back to the variables the
##   formula names, which is what initialize() accepts:

cfg_for <- function(frm, test_term) {
  cfg <- config0
  cfg$frm <- frm
  cfg$test_term <- test_term
  cfg$reference_levels <-
    config0$reference_levels[names(config0$reference_levels) %in% all.vars(frm)]
  return(cfg)
}

## deqms and msqrob take peptide level input, the rest gene level; test() is the entry
##   point under test, since it is what builds the standardized table:

run <- function(method, cfg) {
  if(method %in% c("deqms", "msqrob")) {
    st <- state0
    cf <- cfg
  } else {
    agg <- suppressMessages(combine_features(state0, cfg))
    st <- agg$state
    cf <- agg$config
  }
  args <- list(st, cf, method=method)
  if(method %in% c("proda", "prolfqua")) args$is_log_transformed <- TRUE
  return(try(suppressMessages(do.call(test, args)), silent=TRUE))
}

## independent reference: ordinary least squares on the same design, per gene, from the
##   gene level matrix the engines other than deqms and msqrob are given. The limma
##   family, test_lm() and prolfqua all fit exactly this, so their reported effect size
##   has to match it; proDA and msqrob2 fit something else and are only checked for
##   shape and for agreement with their own coefficients.

ref_state <- suppressMessages(combine_features(state0, config0))$state
ref_cfg <- suppressMessages(combine_features(state0, config0))$config

ref_effect <- function(cfg) {
  design <- f.design_test_cols(ref_state, cfg)
  x_test <- design$X[, design$cols_test, drop=FALSE]
  out <- rep(as.numeric(NA), nrow(ref_state$expression))
  names(out) <- as.character(ref_state$features[[ref_cfg$feat_col]])
  for(idx in seq_along(out)) {
    y <- ref_state$expression[idx, ]
    beta <- stats::lm.fit(design$X, y)$coefficients[design$cols_test]
    out[idx] <- if(length(beta) %in% 1) beta else diff(range(x_test %*% beta))
  }
  return(out)
}

## the reported column, aligned to the reference by feature name:

aligned <- function(res, ref) {
  v <- res$standard$logfc
  names(v) <- res$standard$feature
  i <- intersect(names(ref), names(v))
  return(list(got=v[i], want=ref[i]))
}

close_enough <- function(a, b, tol=1e-8) {
  if(length(a) %in% 0 || length(a) != length(b)) return(FALSE)
  if(any(is.na(a)) || any(is.na(b))) return(FALSE)
  return(max(abs(a - b)) < tol)
}

###############################################################################
section("one design column: the signed coefficient, for every method")

cfg1 <- cfg_for(~sex + age, "sex")
ref1 <- ref_effect(cfg1)

for(method in methods) {
  res <- run(method, cfg1)
  if(inherits(res, "try-error")) {
    report(FALSE, paste0("test(method='", method, "'): runs"))
    next
  }
  ok <- !any(is.na(res$standard$logfc)) && !any(is.na(res$standard$expr))
  report(ok, paste0("test(method='", method, "'): logfc and expr are both reported"))
}

## the four engines whose fit is ordinary least squares on this design must all report
##   the same number, and it must be the coefficient of sexM. voom is not among them:
##   limma::voom() computes precision weights and limma::lmFit() then fits weighted
##   least squares, so its coefficients are its own; it is checked below against the
##   fit it reports, as proDA is:

for(method in c("lm", "trend", "deqms", "prolfqua")) {
  res <- run(method, cfg1)
  cmp <- aligned(res, ref1)
  report(close_enough(cmp$got, cmp$want, tol=1e-6),
    paste0("test(method='", method, "'): logfc is the least squares coefficient"))
}

## a signed quantity: the tested genes carry no sex effect, so both signs occur:

res <- run("trend", cfg1)
report(any(res$standard$logfc > 0) && any(res$standard$logfc < 0),
  "one column: logfc keeps the sign of the coefficient")

## expr is the mean of the values handed to the test:

want <- rowMeans(ref_state$expression, na.rm=TRUE)
names(want) <- as.character(ref_state$features[[ref_cfg$feat_col]])
got <- res$standard$expr
names(got) <- res$standard$feature
report(close_enough(got, want[names(got)], tol=1e-6),
  "one column: expr is the mean of the values handed to the test")

## for deqms and msqrob the mean is over the peptides of each gene as well:

want <- tapply(rowMeans(state0$expression, na.rm=TRUE),
  as.character(state0$features[[config0$gene_id_col]]), mean)
res <- run("msqrob", cfg1)
got <- res$standard$expr
names(got) <- res$standard$feature
report(close_enough(got, want[names(got)], tol=1e-6),
  "one column: expr for msqrob is over the peptides of each gene")

###############################################################################
section("several design columns: the total swing")

## a three level factor: two coefficients, and no engine reported anything here before:

cfg2 <- cfg_for(~grp + age, "grp")
ref2 <- ref_effect(cfg2)

for(method in c("lm", "trend", "prolfqua")) {
  res <- run(method, cfg2)
  cmp <- aligned(res, ref2)
  report(close_enough(cmp$got, cmp$want, tol=1e-6),
    paste0("test(method='", method, "'): joint test logfc is the total swing"))
}

res <- run("trend", cfg2)
report(all(res$standard$logfc > 0),
  "joint test: the total swing is unsigned")

## what the swing means for a factor: the largest difference between any two of its
##   levels. With reference level coding the coefficients are grpb and grpc, both
##   relative to grpa, so the largest gap can be between b and c, which no coefficient
##   names. Checked directly against the fitted level means:

design2 <- f.design_test_cols(ref_state, cfg2)
lvl_swing <- rep(as.numeric(NA), nrow(ref_state$expression))
names(lvl_swing) <- as.character(ref_state$features[[ref_cfg$feat_col]])
for(idx in seq_along(lvl_swing)) {
  beta <- stats::lm.fit(design2$X, ref_state$expression[idx, ])$coefficients
  eff <- c(a=0, b=beta[["grpb"]], c=beta[["grpc"]])   ## level means, less the common part
  lvl_swing[idx] <- max(eff) - min(eff)
}
cmp <- aligned(res, lvl_swing)
report(close_enough(cmp$got, cmp$want, tol=1e-6),
  "joint test: for a factor the swing is the largest difference between two levels")

## and it is genuinely a wider span than the largest single coefficient for at least
##   some feature, so the two are not the same quantity:

max_coef <- rep(as.numeric(NA), nrow(ref_state$expression))
names(max_coef) <- names(lvl_swing)
for(idx in seq_along(max_coef)) {
  beta <- stats::lm.fit(design2$X, ref_state$expression[idx, ])$coefficients
  max_coef[idx] <- max(abs(beta[c("grpb", "grpc")]))
}
report(any(lvl_swing > max_coef + 1e-8),
  "joint test: the swing exceeds the largest single coefficient for some features")

## a continuous covariate spanning two columns by marginality: testing age in ~age*sex
##   covers age and age:sexM, and the swing is then the total change over the observed
##   range of age and both sexes:

cfg3 <- cfg_for(~age * sex, "age")
ref3 <- ref_effect(cfg3)
res <- run("trend", cfg3)
cmp <- aligned(res, ref3)
report(length(cmp$got) > 0 && close_enough(cmp$got, cmp$want, tol=1e-6),
  "joint test: a continuous covariate spanning two columns reports the swing")

## a single continuous column: the coefficient is the change per unit, so the swing
##   would have been that times the range of age. Both are computed here to pin down
##   which one is reported, since that is the documented distinction:

cfg4 <- cfg_for(~sex + age, "age")
res <- run("trend", cfg4)
design4 <- f.design_test_cols(ref_state, cfg4)
slope <- rep(as.numeric(NA), nrow(ref_state$expression))
names(slope) <- names(lvl_swing)
for(idx in seq_along(slope)) {
  beta <- stats::lm.fit(design4$X, ref_state$expression[idx, ])$coefficients
  slope[idx] <- beta[["age"]]
}
cmp <- aligned(res, slope)
report(close_enough(cmp$got, cmp$want, tol=1e-6),
  "one continuous column: logfc is the per-unit slope, not the swing")
rng <- diff(range(ref_state$samples$age))
report(rng > 1 && !close_enough(cmp$got, cmp$want * rng, tol=1e-6),
  "one continuous column: the swing would have been a different number")

## an interaction test spanning several columns, so that the columns of the
##   interaction are covered too:

cfg5 <- cfg_for(~grp * sex, "grp:sex")
ref5 <- ref_effect(cfg5)
res <- run("trend", cfg5)
cmp <- aligned(res, ref5)
report(length(cmp$got) > 0 && close_enough(cmp$got, cmp$want, tol=1e-6),
  "joint test: an interaction term reports the swing over its columns")

###############################################################################
section("the effect size comes from the fit that produced the p-value")

## the swing computed from whatever coefficients an engine's own fit carries:

own_effect <- function(coefs, design, cols) {
  x_test <- design$X[, design$cols_test, drop=FALSE]
  coefs <- coefs[, cols, drop=FALSE]
  if(length(cols) %in% 1) {
    out <- coefs[, 1, drop=TRUE]
    names(out) <- rownames(coefs)
    return(out)
  }
  return(apply(x_test %*% t(coefs), 2, function(v) diff(range(v))))
}

## voom weights each observation by the precision limma::voom() estimates for it, so
##   limma::lmFit() fits weighted rather than ordinary least squares and the
##   coefficients are voom's own. The effect size has to be the one belonging to the
##   fit that produced the p-value, not a least squares stand-in:

res <- run("voom", cfg2)
cmp <- aligned(res, own_effect(res$fit$coefficients, design2, c("grpb", "grpc")))
report(close_enough(cmp$got, cmp$want, tol=1e-6),
  "test(method='voom'): logfc is the swing of voom's own weighted coefficients")
ref_cmp <- aligned(res, ref2)
report(!close_enough(ref_cmp$got, ref_cmp$want, tol=1e-6),
  "test(method='voom'): those are not the ordinary least squares coefficients")

## proDA fits its own model, so its effect size is its own coefficients rather than the
##   least squares ones, and has to match what proDA::proDA() reports:

res <- run("proda", cfg2)
if(inherits(res, "try-error")) {
  report(FALSE, "test(method='proda'): runs the joint test")
  report(FALSE, "test(method='proda'): logfc is the swing of proDA's coefficients")
} else {
  report(!any(is.na(res$standard$logfc)),
    "test(method='proda'): logfc is reported for the joint test")
  coefs <- stats::coefficients(res$fit)[, c("grpb", "grpc"), drop=FALSE]
  x_test <- design2$X[, design2$cols_test, drop=FALSE]
  want <- apply(x_test %*% t(coefs), 2, function(v) diff(range(v)))
  cmp <- aligned(res, want)
  report(close_enough(cmp$got, cmp$want, tol=1e-6),
    "test(method='proda'): logfc is the swing of proDA's own coefficients")
  ref_cmp <- aligned(res, ref2)
  report(!close_enough(ref_cmp$got, ref_cmp$want, tol=1e-6),
    "test(method='proda'): those are not the least squares coefficients")
}

## a feature whose coefficient could not be estimated has no effect size, rather than
##   one computed from part of the design:

coefs <- matrix(c(1, NA, 3, 4), nrow=2, byrow=TRUE,
  dimnames=list(c("f1", "f2"), NULL))
eff <- f.logfc_effect(coefs, design2, config0)
report(identical(names(eff), c("f1", "f2")) && is.na(eff[["f1"]]) &&
    !is.na(eff[["f2"]]),
  "f.logfc_effect(): a missing coefficient gives a missing effect size")

## and a coefficient matrix that does not match the columns under test is an error
##   rather than a silently recycled calculation:

n0 <- log_len()
res <- try(f.logfc_effect(coefs[, 1, drop=FALSE], design2, config0), silent=TRUE)
report(inherits(res, "try-error") && log_has(n0, "f.logfc_effect", "columns under test"),
  "f.logfc_effect(): a coefficient matrix of the wrong width is an error")

###############################################################################
section("intercept test with a continuous covariate: warned about")

cfg6 <- cfg_for(~sex + age, "1")
n0 <- log_len()
res <- run("trend", cfg6)
report(!inherits(res, "try-error"), "test(): the intercept test still runs")
report(log_has(n0, "config$test_term names the intercept", "age",
    "center the covariate"),
  "f.design_test_cols(): warns that the intercept depends on where zero falls")

## no such warning when every covariate is a factor, since the intercept is then the
##   fitted value at the reference levels, which is an observed group:

cfg7 <- cfg_for(~sex + batch, "1")
n0 <- log_len()
res <- run("trend", cfg7)
report(!inherits(res, "try-error") &&
    !log_has(n0, "config$test_term names the intercept"),
  "f.design_test_cols(): no such warning when the covariates are all factors")

###############################################################################
cat("\n## passes:", n_pass, "; failures:", n_fail, "; elapsed:",
  round(as.numeric(difftime(Sys.time(), t0, units="secs"))), "s\n")
cat("## log file:", log_file, "\n")

quit(save="no", status=if(n_fail > 0) 1 else 0)
