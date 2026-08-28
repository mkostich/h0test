## Tests what h0testr::test_voom() does with a feature that carries a missing value.
##   limma::voom() fits one mean-variance trend across the whole matrix and has no
##   handling for a missing value, so test_voom() holds such a feature out of the fit.
##   It used to do that silently: exprs[!i, ] with no count of what went, so the run
##   logged the number of features it tested and nothing anywhere logged the number it
##   did not, and a result table shorter than the state it came from looked like the
##   whole answer. With every feature carrying a missing value the empty matrix reached
##   limma::voom(), which failed from inside on a message about the input to a
##   mean-variance fit rather than about the state.
##   Checked here: that the features held out are exactly the incomplete ones, that the
##   log names how many and the first few of them, that the tested count is reported
##   against the number of features that went in, that the same holds through
##   h0testr::test_h0(), and that a matrix in which nothing is complete is refused with a
##   message saying what to do instead.

usage <- function(msg=NULL) {
  if(!is.null(msg)) cat("ERROR:", msg, "\n\n", file=stderr())
  cat(
    "Test how h0testr::test_voom() handles features carrying missing values: that",
    "the incomplete features are the ones held out of the limma::voom() fit, that",
    "the log says how many were held out and names the first few, that the tested",
    "count is reported against the number of features handed in, that the result of",
    "h0testr::test_h0(method='voom') carries only the tested features, and that a state",
    "in which no feature is complete is refused with a message naming what to do",
    "instead rather than failing inside limma::voom().",
    "",
    "Usage: Rscript test_voom_missing.R <r_dir> [--seed=<integer>]",
    "",
    "Required positional arguments:",
    "  <r_dir>  Path to the h0testr package R/ source directory; all .R files",
    "             there are sourced (the installed package is not used).",
    "",
    "Optional named arguments:",
    "  --seed=<integer>  Seed for the simulated fixture; default 101.",
    "",
    "Output: one 'PASS: <what>' or 'FAIL: <what>' line per assertion on stdout,",
    "  section headers with elapsed seconds, and a final count of passes and",
    "  failures. Exit code 0 if every assertion passed, 1 if any failed, 2 on a",
    "  usage error.",
    "",
    "Examples:",
    "  Rscript test_voom_missing.R C:/path/to/h0testr/R",
    "  Rscript test_voom_missing.R ../../h0test/h0testr/R",
    "  Rscript test_voom_missing.R C:/path/to/h0testr/R --seed=7",
    "",
    sep="\n", file=stderr()
  )
  quit(save="no", status=2)
}

args <- commandArgs(trailingOnly=TRUE)
if(length(args) < 1 || length(args) > 2) usage("wrong number of arguments")

r_dir <- args[1]
if(!dir.exists(r_dir)) usage(paste("r_dir not a directory:", r_dir))

seed <- 101
if(length(args) %in% 2) {
  if(!grepl("^--seed=", args[2])) usage(paste("unrecognized argument:", args[2]))
  seed <- suppressWarnings(as.integer(sub("^--seed=", "", args[2])))
  if(is.na(seed)) usage(paste("--seed not an integer:", args[2]))
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

threw <- function(expr) inherits(try(expr, silent=TRUE), "try-error")
mark <- function() length(readLines(log_file))
logged <- function(pat, since=0) {
  txt <- readLines(log_file)
  if(since >= length(txt)) return(FALSE)
  any(grepl(pat, txt[(since + 1):length(txt)], fixed=TRUE))
}

###############################################################################
## shared data: 20 count-like features over 12 observations in two groups of six, with
##   a crossed second covariate. limma::voom() models a count mean-variance
##   relationship, so the fixture is counts rather than log abundances. Counts are kept
##   above zero because init_state() reads a raw zero as a missing value, and this file
##   is about missing values that are put in deliberately:

set.seed(seed)

nobs <- 12
nfeat <- 20

samps <- data.frame(
  obs=paste0("o", sprintf("%02d", 1:nobs)),
  grp=rep(c("ctl", "trt"), each=nobs / 2),
  sex=rep(c("F", "M"), nobs / 2),
  stringsAsFactors=FALSE
)

exprs <- matrix(stats::rpois(nfeat * nobs, lambda=300) + 20, nrow=nfeat)
exprs[1:6, samps$grp %in% "trt"] <- exprs[1:6, samps$grp %in% "trt"] + 300
rownames(exprs) <- paste0("f", sprintf("%02d", 1:nfeat))
colnames(exprs) <- samps$obs
feats <- data.frame(pep=rownames(exprs), gene=rownames(exprs),
  stringsAsFactors=FALSE)

cfg0 <- new_config()
cfg0$obs_id_col <- cfg0$obs_col <- "obs"
cfg0$sample_id_col <- "obs"
cfg0$feat_id_col <- cfg0$feat_col <- "pep"
cfg0$gene_id_col <- "gene"
cfg0$reference_levels <- c(grp="ctl", sex="F")
cfg0$save_state <- FALSE
cfg0$permute_var <- ""
cfg0$impute_method <- "none"
cfg0$n_features_min <- 5           ## the default of 1000 would filter out everything
cfg0$log_file <- log_file
cfg0$frm <- ~grp + sex
cfg0$test_term <- "grp"

out <- suppressMessages(init_state(list(expression=exprs, features=feats,
  samples=samps), cfg0, minimal=TRUE))
state_full <- out$state
config0 <- out$config

report(!any(is.na(state_full$expression)),
  "fixture: every feature is complete before missing values are put in")

## five features carry a missing value, one of them almost entirely missing, so that
##   "any missing value" and "mostly missing" are not the same set:

gone <- c("f02", "f05", "f11", "f17", "f20")
state_na <- state_full
state_na$expression["f02", 3] <- NA
state_na$expression["f05", c(1, 7)] <- NA
state_na$expression["f11", 12] <- NA
state_na$expression["f17", 2:11] <- NA
state_na$expression["f20", 6] <- NA

kept <- setdiff(rownames(state_na$expression), gone)

report(setequal(gone, rownames(state_na$expression)[apply(state_na$expression, 1,
  function(v) any(is.na(v)))]),
  "fixture: exactly five features carry a missing value")

###############################################################################
section("the features held out are the incomplete ones, and the log says so")

m <- mark()
res <- try(suppressMessages(test_voom(state_na, config0)), silent=TRUE)

report(!inherits(res, "try-error"),
  "test_voom() returns on a state with incomplete features")

if(!inherits(res, "try-error")) {

  report(setequal(rownames(res$hits), kept),
    "the result carries exactly the complete features")
  report(nrow(res$hits) %in% length(kept) &&
    !any(gone %in% rownames(res$hits)),
    "and none of the incomplete ones")

  report(logged("WARNING: test_voom: dropping 5 of 20 features", since=m),
    "the log says how many features were held out, and out of how many")
  report(logged("at least one missing value, which limma::voom() cannot fit",
    since=m),
    "and why they were held out")
  report(logged("f02, f05, f11, f17, f20", since=m),
    "and names the first few of them")
  report(logged("tested 15 of 20 features", since=m),
    "the tested count is reported against the number of features handed in")
}

## the same call on a complete state says nothing about dropping, so the warning is a
##   statement about this state and not a fixed part of the output:

m <- mark()
res_ok <- try(suppressMessages(test_voom(state_full, config0)), silent=TRUE)
report(!inherits(res_ok, "try-error") && nrow(res_ok$hits) %in% nfeat,
  "a complete state gives one row per feature")
report(!logged("WARNING: test_voom: dropping", since=m),
  "and draws no warning about dropped features")
report(logged("tested 20 of 20 features", since=m),
  "and reports the tested count against the same number")

###############################################################################
section("the same holds through test_h0()")

m <- mark()
res_t <- try(suppressMessages(test_h0(state_na, config0, method="voom")), silent=TRUE)

report(!inherits(res_t, "try-error"), "test_h0(method='voom') returns on such a state")

if(!inherits(res_t, "try-error")) {
  report(setequal(as.character(res_t$standard$feature), kept),
    "the standardized table carries exactly the complete features")
  report(logged("WARNING: test_voom: dropping 5 of 20 features", since=m),
    "and the warning reaches the log of a test_h0() run too")
}

###############################################################################
section("a state with too few complete features is refused")

## every feature carrying a missing value used to leave limma::voom() with an empty
##   matrix, which failed from inside it saying "Need at least two genes to fit a
##   mean-variance trend": a statement about the input to a variance fit, with nothing
##   about the missing values that emptied it:

state_none <- state_full
state_none$expression[, 4] <- NA

report(all(apply(state_none$expression, 1, function(v) any(is.na(v)))),
  "fixture: no feature of this state is complete")

m <- mark()
report(threw(suppressMessages(test_voom(state_none, config0))),
  "test_voom() refuses a state in which no feature is complete")
report(logged("test_voom: only 0 of 20 features have no missing value", since=m),
  "and the message says what is wrong with the state")
report(logged("limma::voom() needs at least two to fit a mean-variance trend",
  since=m),
  "and how many it needed")
report(logged("impute first, see h0testr::impute()", since=m),
  "and what to do about it")
report(logged("h0testr::test_lm()", since=m),
  "and names a method that does not need complete features")

## two is limma::voom()'s own bound rather than a round number, so one complete feature
##   is refused for the same reason as none, and two is accepted:

state_one <- state_full
state_one$expression[-1, 4] <- NA
report(sum(apply(state_one$expression, 1, function(v) !any(is.na(v)))) %in% 1,
  "fixture: exactly one feature of this state is complete")

m <- mark()
report(threw(suppressMessages(test_voom(state_one, config0))),
  "one complete feature is refused too, that being fewer than voom can fit")
report(logged("test_voom: only 1 of 20 features have no missing value", since=m),
  "and the message counts the one that was complete")

state_two <- state_full
state_two$expression[-(1:2), 4] <- NA
res_two <- try(suppressMessages(test_voom(state_two, config0)), silent=TRUE)
report(!inherits(res_two, "try-error") && nrow(res_two$hits) %in% 2,
  "two complete features are accepted, so the bound is where voom's is")

## and that method really does handle this state, so the advice is not empty:

res_lm <- try(suppressMessages(test_lm(state_none, config0)), silent=TRUE)
report(!inherits(res_lm, "try-error") && nrow(res_lm$hits) %in% nfeat &&
  any(!is.na(res_lm$hits$pval)),
  "the method it names does test this state")

###############################################################################
cat("\n## passes:", n_pass, "; failures:", n_fail, "; elapsed:",
  round(as.numeric(difftime(Sys.time(), t0, units="secs"))), "s\n")

quit(save="no", status=if(n_fail > 0) 1 else 0)
