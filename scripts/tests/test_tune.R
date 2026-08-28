## Tests three things about h0testr::tune(): that a method name it cannot use is refused
##   before any data is read, that a sweep which tested nothing returns an empty result
##   table rather than NULL, and that the combinations a deleted guard used to skip are now
##   swept.
##   tune() validated test_methods and neither of the other two lists it assigns from. A
##   name that was not an impute_method reached f.err() three loops deep, after every
##   earlier combination had been normalized, aggregated, imputed and tested, and a name
##   that was not a normalization method was not checked at all: it surfaced from inside
##   normalize(), naming the argument it was handed rather than the sweep that handed it
##   over. Both are now checked alongside test_methods, before load_data().
##   The allowed normalization set is not normalize_methods(): it is that plus "q50" and
##   "q75", which normalize() does not accept and f.tune1() maps onto normalization_method
##   "quantile" at a normalization_quantile of 0.5 and 0.75. Checking against
##   normalize_methods() alone would reject tune()'s own defaults, so that is checked here
##   too. "none" stays legal in both lists, unlike for test_methods, not normalizing and
##   not imputing both being baselines a sweep has a use for.
##   And tune() built its result by rbind() onto NULL, so a sweep whose loops never reached
##   f.tune2() returned NULL, which the write.table() in tune()'s own documented usage
##   errors on. It now returns a zero row data.frame carrying the columns a real result row
##   carries.
##   And a guard in the sweep, labelled OIL_WATER, skipped test_method "msqrob" and "voom"
##   whenever the expression matrix carried a missing value, on the reading that neither
##   engine can take one. Neither can be held to that: msqrob2 counts the non-missing
##   values of a feature and weights by them, and test_voom() holds out the incomplete
##   features and says how many. The guard also read the matrix before imputation, so it
##   fired for every impute_method in its branch rather than only where the missing values
##   survive, and it advanced the loop without recording a row, so the combination was
##   absent from the result rather than present with nhits NA. It is gone, and the six
##   cells it used to take two of are checked here.
##   Nothing else in this suite calls tune(); the sweep's inner steps are covered through
##   f.tune2(), f.tune2_na_row() and tune_check() elsewhere.

usage <- function(msg=NULL) {
  if(!is.null(msg)) cat("ERROR:", msg, "\n\n", file=stderr())
  cat(
    "Test what h0testr::tune() does around its sweep: that a name which is not one of",
    "h0testr::impute_methods() or of h0testr::normalize_methods() is refused before any",
    "data is read, that \"q50\", \"q75\" and \"none\" are accepted despite not being",
    "normalize() methods or being a skip, and that a sweep which tested no combination",
    "returns an empty data.frame with the columns of a real result row rather than NULL,",
    "so that writing it out succeeds, and that a sweep over test_methods \"msqrob\" and",
    "\"voom\" on data carrying missing values now records every combination, a guard",
    "having skipped those cells whether or not the imputation step was about to fill",
    "the missing values in.",
    "",
    "Usage: Rscript test_tune.R <r_dir> [--seed=<integer>]",
    "",
    "Required positional arguments:",
    "  <r_dir>  Path to the h0testr package R/ source directory; all .R files",
    "             there are sourced (the installed package is not used). The demo",
    "             data in ../inst/extdata alongside it is used as the sweep input.",
    "",
    "Optional named arguments:",
    "  --seed=<integer>  Seed for the one sweep that runs a test; default 101.",
    "",
    "Output: one 'PASS: <what>' or 'FAIL: <what>' line per assertion on stdout,",
    "  section headers with elapsed seconds, and a final count of passes and",
    "  failures. Exit code 0 if every assertion passed, 1 if any failed, 2 on a",
    "  usage error.",
    "",
    "Examples:",
    "  Rscript test_tune.R C:/path/to/h0testr/R",
    "  Rscript test_tune.R ../../h0test/h0testr/R",
    "  Rscript test_tune.R C:/path/to/h0testr/R --seed=7",
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
invisible(file.create(log_file))   ## mark() reads it before anything here has logged
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
## shared configuration. The demo data the package ships is 30 features over 12
##   observations in two conditions, which is enough for the one sweep below that runs a
##   test. From the source tree rather than system.file(), the suite sourcing R/ rather
##   than loading an installed package, the same way test_check_config.R reaches it:

extdata <- file.path(dirname(r_dir), "inst", "extdata")
if(!dir.exists(extdata)) extdata <- system.file("extdata", package="h0testr")

set.seed(seed)

base_cfg <- function() {
  cfg <- new_config()
  cfg$log_file <- log_file
  cfg$save_state <- FALSE
  cfg$dir_out <- tempdir()
  cfg$dir_in <- extdata
  cfg$feature_file_in <- "features.tsv"
  cfg$sample_file_in <- "samples.tsv"
  cfg$data_file_in <- "expression.tsv"
  cfg$feat_id_col <- cfg$gene_id_col <- "feature_id"
  cfg$obs_id_col <- cfg$sample_id_col <- "observation_id"
  cfg$frm <- ~condition
  cfg$test_term <- "condition"
  cfg$test_method <- "trend"
  cfg$reference_levels <- c(condition="placebo")
  cfg$n_features_min <- 10          ## the default of 1000 would filter out everything
  cfg$permute_var <- ""
  return(cfg)
}

report(dir.exists(extdata) && file.exists(file.path(extdata, "expression.tsv")),
  "the package ships the demo data this file needs")

###############################################################################
section("a bad method name is refused before any data is read")

## the point of validating up front is that nothing has been computed yet when the run
##   stops, so the check is that the refusal happens with dir_in pointing nowhere: if the
##   name were checked after load_data(), the failure would be about the missing file
##   instead. The control below confirms that dir_in really is unusable:

cfg_nodir <- base_cfg()
cfg_nodir$dir_in <- file.path(tempdir(), "no_such_dir_for_test_tune")
report(!dir.exists(cfg_nodir$dir_in), "fixture: dir_in names a directory that is absent")

m <- mark()
report(threw(suppressMessages(tune(cfg_nodir, normalization_methods="log2",
  impute_methods="not_an_impute_method", test_methods="lm"))),
  "tune() refuses a name that is not one of impute_methods()")
report(logged("tune: unexpected impute_methods:", since=m),
  "and the message says which list the name was in")
report(logged("not_an_impute_method", since=m), "and names it")
report(logged("allowed:", since=m), "and lists what is allowed instead")

## load_data() opens by reporting the whole configuration, one "key : value" line per
##   entry, so the absence of those lines is what says the sweep stopped before it. Read
##   failures themselves warn out of utils::read.table() rather than reaching the log, so
##   looking for those would pass whether or not load_data() had run:

report(!logged("data_file_in :", since=m),
  "and stops before load_data(), which would have reported the configuration first")

m <- mark()
report(threw(suppressMessages(tune(cfg_nodir,
  normalization_methods="not_a_norm_method", impute_methods="none",
  test_methods="lm"))),
  "tune() refuses a name that is not a normalization method")
report(logged("tune: unexpected normalization_methods:", since=m),
  "and the message says which list that name was in")
report(logged("not_a_norm_method", since=m), "and names it")
report(!logged("data_file_in :", since=m), "and stops before load_data() here too")

## the neighbouring check the two above were mirrored from, so that all three are known
##   to be reached from the same place:

m <- mark()
report(threw(suppressMessages(tune(cfg_nodir, normalization_methods="log2",
  impute_methods="none", test_methods="not_a_test_method"))),
  "the test_methods check still refuses a name that is not an engine")
report(logged("tune: unexpected test_methods:", since=m),
  "and says which list it was in")

## every name in each enumeration passes, checked as a set rather than one sweep per
##   name: with dir_in unusable a valid list can only fail at load_data(), so what is
##   asserted is that no name drew the "unexpected" refusal:

## read_data() fails by warning and then erroring out of utils::read.table() rather than
##   through f.err(), so the control is that the call threw at all, the log having nothing
##   in it about an unexpected name:

m <- mark()
res_all <- try(suppressWarnings(suppressMessages(tune(cfg_nodir,
  normalization_methods=union(normalize_methods(), c("q50", "q75")),
  impute_methods=impute_methods(), test_methods=test_methods()))), silent=TRUE)
report(!logged("unexpected", since=m),
  "no name any of the three enumerations offers is refused")
report(inherits(res_all, "try-error"),
  "and the run then fails reading the data instead, so the control holds")

## the two that are not normalize() methods, which validating against normalize_methods()
##   alone would have rejected, and which are tune()'s own defaults:

m <- mark()
invisible(try(suppressWarnings(suppressMessages(tune(cfg_nodir,
  normalization_methods=c("q50", "q75"), impute_methods="none",
  test_methods="lm"))), silent=TRUE))
report(!logged("unexpected", since=m),
  "\"q50\" and \"q75\" are accepted, f.tune1() mapping each onto \"quantile\"")

m <- mark()
invisible(try(suppressWarnings(suppressMessages(tune(cfg_nodir,
  normalization_methods="none", impute_methods="none",
  test_methods="lm"))), silent=TRUE))
report(!logged("unexpected", since=m),
  "\"none\" is accepted in both lists, unlike for test_methods, being a baseline")

## and "none" is still refused as a test_method, that check carrying a reason of its own:

m <- mark()
report(threw(suppressMessages(tune(cfg_nodir, normalization_methods="none",
  impute_methods="none", test_methods="none"))),
  "test_methods \"none\" is still refused, a sweep over not testing measuring nothing")

## an empty list is not a bad list: setdiff() of nothing is nothing, and the sweep below
##   depends on an empty impute_methods getting past this check:

m <- mark()
invisible(try(suppressWarnings(suppressMessages(tune(cfg_nodir,
  normalization_methods="log2", impute_methods=character(0),
  test_methods="lm"))), silent=TRUE))
report(!logged("unexpected", since=m), "an empty method list is not refused as a bad name")

###############################################################################
section("a sweep that tested nothing returns an empty table, not NULL")

## an empty impute_methods is the cheapest way to reach the end of the sweep with no rows:
##   the normalization and test_method loops run, so load_data() and f.tune1() both do
##   their work, and the impute loop that calls f.tune2() never iterates. Before this,
##   rslt was still the NULL it started as:

cfg_e <- base_cfg()

m <- mark()
res_e <- try(suppressMessages(tune(cfg_e, normalization_methods="log2",
  impute_methods=character(0), test_methods="lm")), silent=TRUE)

report(!inherits(res_e, "try-error"), "tune() returns on a sweep with no combinations")

if(!inherits(res_e, "try-error")) {

  report(!is.null(res_e), "and what it returns is not NULL")
  report(is.data.frame(res_e), "and is a data.frame")
  report(nrow(res_e) %in% 0, "with no rows, no combination having been tested")
  report(identical(names(res_e), names(f.tune2_na_row(cfg_e))),
    "and the columns of a real result row, in the same order")
  report(logged("no parameter combination was tested", since=m),
    "and the log says no combination was tested")

  ## the failure this fixes: tune()'s own documented usage writes the result out, and
  ##   write.table() errors on NULL:

  f_out <- tempfile(fileext=".tsv")
  ok <- !threw(utils::write.table(res_e, f_out, quote=FALSE, sep="\t", row.names=FALSE))
  report(ok, "and writing it out succeeds, which is what the documented usage does")
  report(ok && file.exists(f_out) && length(readLines(f_out)) %in% 1,
    "leaving a file of just the header row")
  if(file.exists(f_out)) invisible(file.remove(f_out))
}

## and the shape is the shape of a real result, checked against a one cell sweep rather
##   than asserted, so that the empty table cannot drift away from what tune() returns
##   when it does test something:

m <- mark()
res_1 <- try(suppressMessages(tune(cfg_e, normalization_methods="log2",
  impute_methods="none", test_methods="lm")), silent=TRUE)

report(!inherits(res_1, "try-error"), "a one combination sweep returns")

if(!inherits(res_1, "try-error") && !inherits(res_e, "try-error")) {
  report(nrow(res_1) %in% 1, "with one row, that being the one combination")
  report(identical(names(res_1), names(res_e)),
    "and the same columns as the empty table, so the two cannot drift apart")
  report(res_1$norm %in% "log2" && res_1$impute %in% "none" && res_1$test %in% "lm",
    "and the row records the combination that was swept")
  report(!logged("no parameter combination was tested", since=m),
    "and nothing is said about an untested sweep")
}

###############################################################################
section("the combinations the deleted OIL_WATER guard skipped are swept")

## a guard in tune(), labelled OIL_WATER, skipped test_method "msqrob" and "voom" whenever
##   the expression matrix carried a missing value, so those cells were absent from the
##   result rather than present with nhits NA, which is what every other skip in the sweep
##   produces and what tune_check() joins on. Three things were wrong with it. It read the
##   matrix before imputation, so it fired for every impute_method in its branch and not
##   only where the missing values survive: "sample_lod" fills them and the cell was
##   dropped anyway. Neither engine needs it even where they do survive, which is why it
##   was deleted rather than narrowed to impute_method "none". And it advanced the loop
##   without recording anything.
##   The demo data is the whole fixture: it carries missing values where the sweep sees
##   them, and some but not all of its features are complete, so both halves are exercised
##   with nothing simulated.

cfg_g <- base_cfg()

out_g <- suppressMessages(load_data(cfg_g))
e_g <- out_g$state$expression
n_complete <- sum(apply(e_g, 1, function(v) !any(is.na(v))))

report(sum(is.na(e_g)) > 0,
  "fixture: the demo data carries missing values where the sweep sees them")
report(n_complete > 2 && n_complete < nrow(e_g),
  "fixture: some but not all of its features are complete, so voom holds some out")

m <- mark()
res_g <- try(suppressWarnings(suppressMessages(tune(cfg_g,
  normalization_methods="log2", impute_methods=c("none", "sample_lod"),
  test_methods=c("msqrob", "voom", "lm")))), silent=TRUE)

report(!inherits(res_g, "try-error"), "the sweep over both engines returns")

if(!inherits(res_g, "try-error")) {

  cell <- function(im, tm) res_g[res_g$impute %in% im & res_g$test %in% tm, , drop=FALSE]

  report(nrow(res_g) %in% 6,
    "every one of the six combinations is in the result, none skipped")
  report(!any(is.na(res_g$ntests)) && !any(is.na(res_g$nhits)),
    "and every one of them ran, none recorded as untested")
  report(all(sapply(c("none", "sample_lod"), function(im)
      all(sapply(c("msqrob", "voom", "lm"), function(tm) nrow(cell(im, tm)) %in% 1)))),
    "one row per combination, keyed by impute and test rather than by position")

  ## msqrob2 counts the non-missing values of a feature and weights by them, so msqrob
  ##   tests every feature whether or not it is complete. This is why the guard could not
  ##   just be narrowed to impute_method "none", asserted here rather than argued:

  report(cell("none", "msqrob")$ntests %in% nrow(e_g),
    "msqrob tests every feature with the missing values left in, needing no imputation")
  report(cell("sample_lod", "msqrob")$ntests %in% nrow(e_g),
    "and every feature with them filled, so the guard cost it both of its cells")

  ## voom holds out the incomplete features and says how many, which is a different thing
  ##   from the cell being skipped: the row is there, and reports what it tested:

  report(cell("none", "voom")$ntests %in% n_complete,
    "voom tests just the complete features with the missing values left in")
  report(cell("sample_lod", "voom")$ntests %in% nrow(e_g),
    "and every feature once imputation has filled them")

  ## the control: lm was never guarded, so its two cells are what the other four now look
  ##   like rather than being the exception:

  report(cell("none", "lm")$ntests %in% nrow(e_g) &&
      cell("sample_lod", "lm")$ntests %in% nrow(e_g),
    "the unguarded control method tested every feature in both of its cells")

  report(!logged("OIL_WATER", since=m), "and nothing in the log mentions the guard")
}

###############################################################################
cat("\n## passes:", n_pass, "; failures:", n_fail, "; elapsed:",
  round(as.numeric(difftime(Sys.time(), t0, units="secs"))), "s\n")

quit(save="no", status=if(n_fail > 0) 1 else 0)
