## Tests how h0testr handles missing, non-finite and blank covariate values, which is
##   the one covariate-value case the formula/test_term matrix never exercised.
##   stats::model.matrix() drops observations with missing covariate values by default,
##   which would leave the design with fewer rows than state$expression has columns and
##   silently pair observations with the wrong covariates, so every path that builds a
##   design has to refuse such a covariate rather than fit it. Covers init_state(), the
##   two formula-aware filters, f.design_test_cols() and all seven test_*() methods,
##   both for a state that arrives with missing values and for a state whose covariates
##   are edited after init_state() has already checked them.
##   Blank values are covered as well: utils::read.table() reads an empty field in a
##   character column as "" rather than as NA, so an empty cell in the samples file
##   would otherwise become a factor level of its own.

usage <- function(msg=NULL) {
  if(!is.null(msg)) cat("ERROR:", msg, "\n\n", file=stderr())
  cat(
    "Test that missing (NA, NaN), non-finite (Inf) and blank ('') covariate values",
    "are refused rather than silently changing which observations are fit. Checks",
    "that init_state() rejects them up front, naming the covariate and the",
    "offending observations; that a covariate edited to contain NA after",
    "init_state() is still refused by f.design_test_cols(),",
    "filter_features_by_formula(), filter_features_by_estimability() and by each of",
    "test_lm(), test_trend(), test_voom(), test_deqms(), test_msqrob(),",
    "test_proda() and test_prolfqua(); and that blank strings, which",
    "utils::read.table() leaves as '' in a character column rather than reading as",
    "NA, are not accepted as a factor level.",
    "",
    "Requires the prolfqua package to be installed.",
    "",
    "Usage: Rscript test_na_covariates.R <r_dir>",
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
    "  Rscript test_na_covariates.R C:/path/to/h0test/h0testr/R",
    "  Rscript test_na_covariates.R ../../h0test/h0testr/R",
    "  Rscript test_na_covariates.R C:/path/to/h0testr/R > test_na_covariates.out 2>&1",
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

## f.err() writes the informative message to config$log_file and then stops with a
##   bare "Stopping", so what a message says has to be read back out of the log; only
##   the lines an attempt added are considered:

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
## Twelve observations, sex and batch crossed with three replicates each, age
##   continuous. Eight genes with two peptides each, two peptides dropped so that
##   DEqMS sees varying peptide counts. No missing expression values: the subject here
##   is missingness in state$samples, not in state$expression.

set.seed(101)

nobs <- 12
ngene <- 8
npep <- 2

samps <- data.frame(
  obs=paste0("o", sprintf("%02d", 1:nobs)),
  sex=rep(c("F", "M"), each=nobs / 2),
  batch=rep(c("b1", "b2"), times=nobs / 2),
  age=round(rnorm(nobs, 50, 10), 1),
  stringsAsFactors=FALSE
)

gene_eff <- matrix(0, nrow=ngene, ncol=nobs)
gene_eff[1:4, samps$sex %in% "M"] <- 3

exprs <- matrix(rnorm(ngene * npep * nobs, 20, 1.5), nrow=ngene * npep)
exprs <- exprs + gene_eff[rep(1:ngene, each=npep), , drop=FALSE]
gene <- paste0("g", sprintf("%02d", rep(1:ngene, each=npep)))
rownames(exprs) <- paste0(gene, ".p", rep(1:npep, times=ngene))
colnames(exprs) <- samps$obs

## the peptides per gene have to vary for test_method "deqms", and to be spread rather
##   than merely varied: DEqMS fits its variance prior as a loess of residual variance on
##   that count, and dropping two rows overall left all but two genes at one count, where
##   the loess comes back NaN and those genes lose their moderated statistic altogether.
##   Dropping the first peptide of every other gene leaves half the genes with one peptide
##   and half with two, which it fits:

i_drop <- seq(from=1, by=2 * npep, length.out=ngene %/% 2)
exprs <- exprs[-i_drop, , drop=FALSE]
gene <- gene[-i_drop]

feats <- data.frame(pep=rownames(exprs), gene=gene)

cfg0 <- new_config()
cfg0$obs_id_col <- cfg0$obs_col <- "obs"
cfg0$sample_id_col <- "obs"
cfg0$feat_id_col <- cfg0$feat_col <- "pep"
cfg0$gene_id_col <- "gene"
## age is in config$frm so that a continuous covariate is checked too, and so that
##   f.subset_covariates() keeps the column in state$samples:
cfg0$frm <- ~sex + batch + age
cfg0$test_term <- "sex"
cfg0$reference_levels <- c(sex="F", batch="b1")
cfg0$estimability <- "test"
cfg0$df_resid_min <- 2
cfg0$save_state <- FALSE
cfg0$permute_var <- ""
cfg0$is_log_transformed <- TRUE     ## values above are already on a log scale
cfg0$impute_method <- "none"
cfg0$n_features_min <- 5            ## default of 1000 would filter out every observation
cfg0$log_file <- log_file

state_in <- list(expression=exprs, features=feats, samples=samps)

out <- init_state(state_in, cfg0, minimal=TRUE)
state0 <- out$state
config0 <- out$config

## deqms and msqrob take peptide-level input; the rest take gene-level:

run_test <- function(nm, st, cfg) {
  fn <- get(paste0("test_", nm))
  if(!(nm %in% c("deqms", "msqrob"))) {
    agg <- try(suppressMessages(combine_features(st, cfg)), silent=TRUE)
    if(inherits(agg, "try-error")) return(agg)
    st <- agg$state
    cfg <- agg$config
  }
  args <- list(st, cfg)
  if(nm %in% "proda") args <- list(st, cfg, is_log_transformed=TRUE, prior_df=5)
  if(nm %in% "prolfqua") args <- list(st, cfg, is_log_transformed=TRUE)
  return(try(suppressMessages(do.call(fn, args)), silent=TRUE))
}

methods <- c("lm", "trend", "voom", "deqms", "msqrob", "proda", "prolfqua")

###############################################################################
section("init_state(): missing and non-finite covariate values are refused")

## a missing value in a continuous covariate:

st <- state_in
st$samples$age[c(2, 5)] <- NA
n0 <- log_len()
res <- try(suppressMessages(init_state(st, cfg0, minimal=TRUE)), silent=TRUE)
report(inherits(res, "try-error"),
  "init_state(): NA in a continuous covariate is an error")
report(log_has(n0, "f.check_covariate_values", "age", "n missing: 2"),
  "init_state(): the error names the covariate and how many values are missing")
report(log_has(n0, "o02", "o05"),
  "init_state(): the error names the offending observations")

## a missing value in a factor covariate declared in config$reference_levels:

st <- state_in
st$samples$sex[4] <- NA
n0 <- log_len()
res <- try(suppressMessages(init_state(st, cfg0, minimal=TRUE)), silent=TRUE)
report(inherits(res, "try-error") &&
    log_has(n0, "f.check_covariate_values", "sex", "n missing: 1", "o04"),
  "init_state(): NA in a declared factor covariate is an error naming it")

## a missing value in a covariate that is already a factor, so is classified without
##   any declaration in config$reference_levels:

st <- state_in
st$samples$batch <- factor(st$samples$batch)
st$samples$batch[6] <- NA
cfg <- cfg0
cfg$reference_levels <- cfg0$reference_levels["sex"]
n0 <- log_len()
res <- try(suppressMessages(init_state(st, cfg, minimal=TRUE)), silent=TRUE)
report(inherits(res, "try-error") &&
    log_has(n0, "f.check_covariate_values", "batch", "n missing: 1", "o06"),
  "init_state(): NA in an undeclared factor covariate is an error naming it")

## NaN is missing; Inf is not missing but is not fittable either, and is reported
##   separately so that the message matches what is actually in the column:

st <- state_in
st$samples$age[3] <- NaN
n0 <- log_len()
res <- try(suppressMessages(init_state(st, cfg0, minimal=TRUE)), silent=TRUE)
report(inherits(res, "try-error") && log_has(n0, "has missing values", "age", "o03"),
  "init_state(): NaN in a continuous covariate is reported as missing")

st <- state_in
st$samples$age[3] <- Inf
n0 <- log_len()
res <- try(suppressMessages(init_state(st, cfg0, minimal=TRUE)), silent=TRUE)
report(inherits(res, "try-error") &&
    log_has(n0, "has non-finite", "age", "n non-finite: 1", "o03"),
  "init_state(): Inf in a continuous covariate is reported as non-finite")

## an undeclared character covariate is refused for its class before its missing
##   values are reached, since the reference level has to be settled first; once it is
##   declared, the missing values are what stops it:

st <- state_in
st$samples$batch[2] <- NA
cfg <- cfg0
cfg$reference_levels <- cfg0$reference_levels["sex"]
n0 <- log_len()
res <- try(suppressMessages(init_state(st, cfg, minimal=TRUE)), silent=TRUE)
report(inherits(res, "try-error") &&
    log_has(n0, "f.covariate_types", "batch", "config$reference_levels"),
  "init_state(): an undeclared character covariate is refused for its class first")

n0 <- log_len()
res <- try(suppressMessages(init_state(st, cfg0, minimal=TRUE)), silent=TRUE)
report(inherits(res, "try-error") &&
    log_has(n0, "f.check_covariate_values", "batch", "n missing: 1"),
  "init_state(): once declared, its missing values are what stops it")

###############################################################################
section("control: every method runs on the same fixture with no missing covariates")

for(nm in methods) {
  res <- run_test(nm, state0, config0)
  report(!inherits(res, "try-error"),
    paste0("test_", nm, "(): runs with complete covariates"))
}

###############################################################################
section("covariates edited after init_state(): every design-building path refuses")

## init_state() checks the covariates once, so a state edited afterwards (by hand, or
##   by a script that merges in more sample annotation) can still reach the fits with
##   missing values in it. Each of the paths below builds a design, and each has to
##   refuse rather than fit the observations model.matrix() leaves:

state_na <- state0
state_na$samples$age[c(2, 5)] <- NA
config_na <- config0

n0 <- log_len()
res <- try(f.design_test_cols(state_na, config_na), silent=TRUE)
report(inherits(res, "try-error"),
  "f.design_test_cols(): a covariate with NA is an error, not a shortened design")
report(log_has(n0, "design matrix has", "10", "state$expression has", "12"),
  "f.design_test_cols(): the error gives both row counts")
report(log_has(n0, "drops observations with missing covariate values"),
  "f.design_test_cols(): the error says why the counts differ")

n0 <- log_len()
res <- try(suppressMessages(
  filter_features_by_estimability(state_na, config_na)), silent=TRUE)
report(inherits(res, "try-error") && log_has(n0, "design matrix has"),
  "filter_features_by_estimability(): refuses a covariate with NA")

for(nm in methods) {
  res <- run_test(nm, state_na, config_na)
  report(inherits(res, "try-error"),
    paste0("test_", nm, "(): refuses a covariate with NA"))
}

###############################################################################
section("filter_features_by_formula(): counts only over observations that are fit")

## this filter counts observations rather than building a design, so on its own it
##   would complete on a state that no fit can use, reporting a count taken over the
##   observations whose covariates happen to be known (base::table() drops NA
##   silently). It checks the covariate values for that reason:

n0 <- log_len()
res <- try(suppressMessages(
  filter_features_by_formula(state_na, config_na)), silent=TRUE)
report(inherits(res, "try-error"),
  "filter_features_by_formula(): refuses a covariate with NA")
report(log_has(n0, "filter_features_by_formula", "age", "has missing values"),
  "filter_features_by_formula(): the error names itself and the covariate")

## the same check makes a constant covariate an error here, rather than a screen that
##   silently drops every feature for want of a second level to count:

st <- state0
st$samples$batch[] <- factor("b1", levels=levels(st$samples$batch))
n0 <- log_len()
res <- try(suppressMessages(filter_features_by_formula(st, config0)), silent=TRUE)
report(inherits(res, "try-error") &&
    log_has(n0, "filter_features_by_formula", "batch", "is constant"),
  "filter_features_by_formula(): refuses a constant covariate")

## init_state() has already warned about a continuous covariate with few distinct
##   values, so this filter does not repeat that warning:

st <- state0
st$samples$age <- rep(c(40, 60), each=nrow(st$samples) / 2)
n0 <- log_len()
res <- try(suppressMessages(filter_features_by_formula(st, config0)), silent=TRUE)
report(!inherits(res, "try-error") && !log_has(n0, "distinct values"),
  "filter_features_by_formula(): does not repeat the distinct-value warning")

###############################################################################
section("blank strings: '' is not a covariate value")

## utils::read.table(as.is=TRUE), which read_data() uses, reads an empty field in a
##   character column as "" rather than as NA (only na.strings, default "NA", becomes
##   NA), so a samples.tsv with a blank cell arrives as a covariate with a "" value.
##   That is neither NA nor non-finite, so it is checked separately:

tsv <- tempfile(fileext=".tsv")
writeLines(c("obs\tsex\tage", "o01\tF\t50", "o02\t\t60", "o03\tM\t"), tsv)
tmp <- utils::read.table(tsv, header=TRUE, sep="\t", quote="", as.is=TRUE)
unlink(tsv)

report(identical(tmp$sex[2], "") && is.na(tmp$age[3]),
  "read.table(): a blank character field is '', a blank numeric field is NA")

st <- state_in
st$samples$sex[4] <- ""
n0 <- log_len()
res <- try(suppressMessages(init_state(st, cfg0, minimal=TRUE)), silent=TRUE)
report(inherits(res, "try-error") &&
    log_has(n0, "f.check_covariate_values", "sex", "has blank values", "n blank: 1",
      "o04"),
  "init_state(): a blank factor value is refused, not made a level")
report(log_has(n0, "empty cell in the samples file"),
  "init_state(): the error explains where a blank value comes from")

## whitespace only is a blank as well; a value that only looks empty is the harder
##   one to notice by eye in the samples file:

st <- state_in
st$samples$sex[4] <- "  "
n0 <- log_len()
res <- try(suppressMessages(init_state(st, cfg0, minimal=TRUE)), silent=TRUE)
report(inherits(res, "try-error") && log_has(n0, "has blank values", "sex"),
  "init_state(): a whitespace-only factor value is refused too")

## a blank in a continuous covariate cannot arise from read.table(), which reads an
##   empty numeric field as NA; if one is constructed anyway the column is character,
##   so it is refused for its class before its values are reached:

st <- state_in
st$samples$age <- as.character(st$samples$age)
st$samples$age[4] <- ""
n0 <- log_len()
res <- try(suppressMessages(init_state(st, cfg0, minimal=TRUE)), silent=TRUE)
report(inherits(res, "try-error") &&
    log_has(n0, "f.covariate_types", "age", "config$reference_levels"),
  "init_state(): a blank in a character-valued numeric covariate is refused")

###############################################################################
cat("\n## passes:", n_pass, "; failures:", n_fail, "; elapsed:",
  round(as.numeric(difftime(Sys.time(), t0, units="secs"))), "s\n")
cat("## log file:", log_file, "\n")

quit(save="no", status=if(n_fail > 0) 1 else 0)
