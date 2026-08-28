usage <- function(msg=NULL) {
  if(!is.null(msg)) cat("ERROR:", msg, "\n\n", file=stderr())
  cat(
    "Test what h0testr::tune() does when a normalization method refuses the data: that",
    "f.tune_norm() maps the q50, q75 and upperquartile names onto a method and a quantile,",
    "that a refused normalization costs its own cells rather than the sweep, that those",
    "cells carry nhits and ntests NA with step 'normalize', that the surviving",
    "normalizations are tested and report join compatible norm and nquant values, and that",
    "a sweep in which every normalization is refused is an error rather than a silent table",
    "of NA.",
    "",
    "Usage: Rscript test_norm_refused.R <r_dir>",
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
    "  Rscript test_norm_refused.R C:/path/to/h0test/h0testr/R",
    "  Rscript test_norm_refused.R ../../h0test/h0testr/R",
    "  Rscript test_norm_refused.R C:/path/to/h0testr/R > t.out 2>&1",
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
## logged() reads this before anything has written to it:
invisible(file.create(log_file))
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

logged <- function(pattern, since=0) {
  txt <- readLines(log_file)
  if(since > 0) {
    if(length(txt) <= since) return(FALSE)
    txt <- txt[(since + 1):length(txt)]
  }
  any(grepl(pattern, txt, fixed=TRUE))
}

cfg0 <- new_config()
cfg0$save_state <- FALSE
cfg0$dir_in <- file.path(dirname(r_dir), "inst", "extdata")
cfg0$feature_file_in <- "features.tsv"
cfg0$sample_file_in <- "samples.tsv"
cfg0$data_file_in <- "expression.tsv"
cfg0$feat_id_col <- cfg0$gene_id_col <- "feature_id"
cfg0$obs_id_col <- cfg0$sample_id_col <- "observation_id"
cfg0$n_features_min <- 10
cfg0$frm <- ~condition
cfg0$test_term <- "condition"
cfg0$reference_levels <- c(condition="placebo")
cfg0$permute_var <- ""
cfg0$log_file <- log_file

if(!dir.exists(cfg0$dir_in)) usage(paste("extdata not found:", cfg0$dir_in))

## a name tune() accepts and normalize() refuses on this fixture:
bad_norm <- "quantiles.robust"

###############################################################################
section("f.tune_norm() name mapping")

cfg <- cfg0
cfg$normalization_quantile <- 0.9

c50 <- f.tune_norm(cfg, "q50")
report(c50$normalization_method %in% "quantile" && c50$normalization_quantile %in% 0.5,
  "q50 maps onto method quantile at 0.50")

c75 <- f.tune_norm(cfg, "q75")
report(c75$normalization_method %in% "quantile" && c75$normalization_quantile %in% 0.75,
  "q75 maps onto method quantile at 0.75")

cuq <- f.tune_norm(cfg, "upperquartile")
report(cuq$normalization_method %in% "upperquartile" &&
  cuq$normalization_quantile %in% 0.75,
  "upperquartile keeps its name and gets the quantile it needs")

crle <- f.tune_norm(cfg, "RLE")
report(crle$normalization_method %in% "RLE" && crle$normalization_quantile %in% 0.9,
  "a method needing no quantile passes through with the config value untouched")

report(bad_norm %in% union(normalize_methods(), c("q50", "q75")),
  paste(bad_norm, "is a name tune() accepts, so the sweep reaches normalize()"))

###############################################################################
section("a refused normalization costs its own cells")

norms <- c(bad_norm, "none", "q50")
imputes <- c("sample_lod", "none")

m <- mark()
set.seed(101)
res <- try(suppressMessages(tune(cfg0, normalization_methods=norms,
  impute_methods=imputes, test_methods="lm")), silent=TRUE)

report(!inherits(res, "try-error"), "the sweep returns rather than failing")

if(!inherits(res, "try-error")) {

  report(nrow(res) %in% (length(norms) * length(imputes)),
    "every combination has a row, the refused ones included")

  i_bad <- res$norm %in% bad_norm
  report(sum(i_bad) %in% length(imputes),
    "the refused normalization has one row per imputation")

  report(all(is.na(res$nhits[i_bad])) && all(is.na(res$ntests[i_bad])),
    "and they carry no hits and no tests")
  report(all(res$step[i_bad] %in% "normalize"),
    "their step names normalization, not the step below it")
  report(all(!is.na(res$reason[i_bad])), "and a reason is recorded")

  report(all(!is.na(res$ntests[!i_bad])) && all(res$ntests[!i_bad] > 0),
    "every combination under a normalization that worked was tested")
  report(all(is.na(res$step[!i_bad])) && all(is.na(res$reason[!i_bad])),
    "and a tested row leaves step and reason empty")

  report(logged(paste("normalization_method", bad_norm, "failed"), since=m),
    "the log says which normalization was refused")

  i_q50 <- res$norm %in% "quantile"
  report(any(i_q50) && all(res$nquant[i_q50] %in% 0.5),
    "a q50 row reports norm quantile and nquant 0.50, as a tested row does")
  report(!any(res$norm %in% c("q50", "q75")),
    "so tune_check() joins the swept name to the mapped one and not the reverse")
}

###############################################################################
section("a sweep with nothing left to test")

m <- mark()
res2 <- try(suppressMessages(tune(cfg0, normalization_methods=bad_norm,
  impute_methods="none", test_methods="lm")), silent=TRUE)

report(inherits(res2, "try-error"),
  "a sweep in which every normalization is refused is an error")
report(logged("no parameter combination could be tested", since=m),
  "rather than a table of NA the caller might rank")

###############################################################################

cat("\n## log file:", log_file, "\n")
cat("\n## passes:", n_pass, "; failures:", n_fail, "; elapsed:",
  round(as.numeric(difftime(Sys.time(), t0, units="secs"))), "s\n")

quit(save="no", status=if(n_fail > 0) 1 else 0)
