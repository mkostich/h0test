## Tests for the capability pre-checks h0testr::check_config() performs, which refuse at
##   configuration time the combinations that no other setting can rescue: a
##   config$test_method that is not one of h0testr::test_methods(), and
##   config$test_ridge=TRUE on a mean model whose non-intercept column count the config
##   already settles as fewer than two. Also covers f.frm_min_noint_cols(), the helper that
##   bounds that column count, and f.note_ignored_settings(), the note h0testr::test()
##   makes about a setting the method it resolved does not consult, which lives there
##   rather than in check_config() so that it is said once per run and not once per step.
##   The refusals that depend on the data rather than on the configuration stay with the
##   engines and are tested there.

usage <- function(msg=NULL) {
  if(!is.null(msg)) cat("ERROR:", msg, "\n\n", file=stderr())
  cat(
    "Test the capability pre-checks in h0testr::check_config(): that a test_method",
    "outside h0testr::test_methods() is refused, that config$test_ridge=TRUE is refused",
    "on a mean model the config already settles as having fewer than two non-intercept",
    "columns and accepted when the count is two or unknown, that a setting the method",
    "h0testr::test() resolved does not consult is reported by it as a NOTE rather than",
    "refused and only when it differs from the h0testr::new_config() default, that a",
    "config$run_order naming a step h0testr::run() has no function for is refused while a",
    "step named twice is warned about by run() instead, and that the checks that were",
    "there before still fire.",
    "",
    "Usage: Rscript test_check_config.R <r_dir>",
    "",
    "Required positional arguments:",
    "  <r_dir>  Path to the h0testr package R/ source directory; all .R files",
    "             there are sourced (the installed package is not used).",
    "",
    "Optional named arguments: none.",
    "",
    "Needs no packages beyond those h0testr itself loads: the only engine run is the",
    "  limma fit behind test_method 'trend', to confirm the note reaches the log and to",
    "  walk one whole workflow through h0testr::run() on the demo data the package ships.",
    "",
    "Output: one PASS/FAIL line per assertion to stdout, then a count of passes and",
    "  failures and the elapsed time. Messages from expected errors are written to a",
    "  temporary log file, whose path is reported at the end.",
    "",
    "Exit codes: 0 all assertions passed; 1 one or more failed; 2 usage error.",
    "",
    "Examples:",
    "  Rscript test_check_config.R C:/path/to/h0test/h0testr/R",
    "  Rscript test_check_config.R ../../h0test/h0testr/R",
    "  Rscript test_check_config.R C:/path/to/h0testr/R > cfg.out 2>&1",
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

mark <- function() if(file.exists(log_file)) length(readLines(log_file)) else 0

log_since <- function(since=0) {
  if(!file.exists(log_file)) return(character(0))
  txt <- readLines(log_file)
  if(since >= length(txt)) return(character(0))
  return(txt[(since + 1):length(txt)])
}

logged <- function(pat, since=0) any(grepl(pat, log_since(since), fixed=TRUE))
log_count <- function(pat, since=0) sum(grepl(pat, log_since(since), fixed=TRUE))

## every config below logs to the same file, so that an expected error's message can be
##   read back rather than only its having stopped:

## save_state FALSE and dir_out a temporary directory for every fixture here, not just the
##   ones that run a workflow: new_config() ships save_state TRUE and dir_out ".", so a
##   fixture that reaches test() writes its results table into the working directory, named
##   after length(config$run_order) + 3 and so easy to mistake for a real run's output:

base_cfg <- function() {
  cfg <- new_config()
  cfg$log_file <- log_file
  cfg$save_state <- FALSE
  cfg$dir_out <- tempdir()
  return(cfg)
}

###############################################################################
section("config$test_method against test_methods()")

report(isTRUE(check_config(base_cfg())), "the default config passes")

ok <- TRUE
for(nm in test_methods()) {
  cfg <- base_cfg()
  cfg$test_method <- nm
  if(!isTRUE(try(check_config(cfg), silent=TRUE))) ok <- FALSE
}
report(ok, "every name h0testr::test_methods() offers passes")

cfg <- base_cfg()
cfg$test_method <- ""
report(isTRUE(check_config(cfg)),
  "an empty test_method passes, meaning unset for test(method=) to supply")

cfg <- base_cfg()
cfg$test_method <- "trrend"
m0 <- mark()
report(threw(check_config(cfg)), "a test_method that is not one of them is refused")
report(logged("unexpected test_method: trrend", m0),
  "and the message names the offending value")
report(logged("prolfqua_lmer", m0), "and lists the ones that are allowed")

cfg <- base_cfg()
cfg$test_method <- NA_character_
report(threw(check_config(cfg)),
  "an NA test_method is refused, absence being how a setting is left out")

## the refusal has to be reachable from the config alone, so that a sweep setting one
##   config and varying the method still fails at the first step rather than at the test:

cfg <- base_cfg()
cfg$test_method <- "no_such_method"
m0 <- mark()
report(threw(check_config(cfg)) && logged("no_such_method", m0),
  "the same refusal reaches any function that calls check_config() first")

###############################################################################
section("f.frm_min_noint_cols(), the column count the config settles")

cfg <- base_cfg()
cfg$frm <- ~1
report(f.frm_min_noint_cols(cfg) %in% 0, "no term at all is no non-intercept column")

cfg$frm <- ~grp + sex
report(f.frm_min_noint_cols(cfg) %in% 2,
  "two terms are at least two columns, whatever their levels")

cfg$frm <- ~grp
report(is.na(f.frm_min_noint_cols(cfg)),
  "one unresolved term is unknown, its column count being data")

cfg$covariate_types <- c(grp="numeric")
report(f.frm_min_noint_cols(cfg) %in% 1, "a numeric term is one column")

cfg$covariate_types <- c(grp="factor")
cfg$factor_levels <- list(grp=c("ctl", "trt"))
report(f.frm_min_noint_cols(cfg) %in% 1,
  "a two level factor with an intercept is one column")

cfg$frm <- ~0 + grp
report(f.frm_min_noint_cols(cfg) %in% 2,
  "and two columns when the intercept is suppressed")

cfg$frm <- ~grp
cfg$factor_levels <- list(grp=c("ctl", "trt", "hi"))
report(f.frm_min_noint_cols(cfg) %in% 2,
  "a three level factor with an intercept is two columns")

cfg$frm <- ~grp:sex
cfg$covariate_types <- c(grp="factor", sex="factor")
report(is.na(f.frm_min_noint_cols(cfg)),
  "an interaction alone is left unknown rather than guessed at")

###############################################################################
section("config$test_ridge against the mean model")

ridge_cfg <- function(frm, levs=c("ctl", "trt")) {
  cfg <- base_cfg()
  cfg$test_method <- "msqrob_agg"
  cfg$test_ridge <- TRUE
  cfg$frm <- frm
  cfg$test_term <- "grp"
  cfg$covariate_types <- c(grp="factor")
  cfg$factor_levels <- list(grp=levs)
  return(cfg)
}

m0 <- mark()
report(threw(check_config(ridge_cfg(~grp))),
  "the penalty is refused on a two level factor that keeps its intercept")
report(logged("at least two non-intercept columns", m0),
  "and the message says what the penalty needs")
report(logged("~0 + ...", m0), "and how to write the same model with two columns")

report(isTRUE(check_config(ridge_cfg(~0 + grp))),
  "and accepted once the intercept is suppressed")

report(isTRUE(check_config(ridge_cfg(~grp, levs=c("ctl", "trt", "hi")))),
  "and accepted on a three level factor with an intercept")

cfg <- ridge_cfg(~grp + sex)
cfg$covariate_types <- c(grp="factor", sex="factor")
report(isTRUE(check_config(cfg)), "and accepted on a second term, whatever its levels")

## before initialize() has classified the covariates the count is not knowable from the
##   config, so the combination is left to msqrob2 rather than refused on a guess:

cfg <- ridge_cfg(~grp)
cfg$covariate_types <- NULL
cfg$factor_levels <- NULL
report(isTRUE(check_config(cfg)),
  "an unresolved single term is left to msqrob2 rather than refused")

cfg <- ridge_cfg(~1)
cfg$test_term <- ""
m0 <- mark()
report(threw(check_config(cfg)), "a mean model with no covariate at all is refused")
report(logged("a covariate for the penalty to apply to", m0),
  "and that message asks for a covariate rather than for a suppressed intercept")

cfg <- ridge_cfg(~grp)
cfg$covariate_types <- c(grp="numeric")
report(threw(check_config(cfg)), "a single numeric term is refused, being one column")

cfg <- ridge_cfg(~grp)
cfg$test_ridge <- FALSE
report(isTRUE(check_config(cfg)), "the same model passes with the penalty off")

## only the one method reads the setting, so on any other the model is beside the point
##   and the setting being ignored is what to say instead:

cfg <- ridge_cfg(~grp)
cfg$test_method <- "msqrob"
m0 <- mark()
report(isTRUE(check_config(cfg)),
  "test_ridge on a method that does not read it is not a refusal")
report(!logged("does not consult", m0),
  "and check_config() says nothing about it either, that being test()'s to say")
report(length(f.note_ignored_settings("msqrob", cfg)) %in% 1,
  "which it does, as a note that the setting has no effect")

###############################################################################
section("settings the chosen method does not consult")

## a setting only some engines read is not refused, and is not check_config()'s business:
##   every step of the workflow calls check_config(), which would repeat the note once per
##   step, so test() makes it once per run instead. f.note_ignored_settings() returns what
##   it said, so most of this needs no fit; the run of test() at the end of the section
##   confirms that the note reaches the log from there:

cfg <- base_cfg()
report(length(f.note_ignored_settings("trend", cfg)) %in% 0,
  "a default config has nothing to say")

## config$test_trend is deliberately not in this table: test() resolves it for every
##   method, its trend argument overriding it, so what to say about a method that cannot
##   honor it depends on the resolved value rather than on the config alone. That is
##   f.note_trend()'s, and 1/test_trend_flag.R covers it:

cfg <- base_cfg()
cfg$test_trend <- TRUE
report(length(f.note_ignored_settings("lm", cfg)) %in% 0,
  "test_trend is not this note's business, being resolved by test() instead")

cfg <- base_cfg()
report(length(f.note_ignored_settings("trend", cfg)) %in% 0,
  "the default test_random_obs is not noted, having not been asked for")

cfg <- base_cfg()
cfg$test_random_obs <- FALSE
got <- f.note_ignored_settings("trend", cfg)
report(length(got) %in% 1 && grepl("config$test_random_obs=FALSE", got, fixed=TRUE),
  "but a non-default test_random_obs on a least squares method is")

for(nm in c("prolfqua_lmer", "msqrob_agg")) {
  cfg <- base_cfg()
  cfg$test_random_obs <- FALSE
  report(length(f.note_ignored_settings(nm, cfg)) %in% 0,
    paste("and not on method", nm, "which reads it"))
}

cfg <- base_cfg()
cfg$test_ridge <- TRUE
report(length(f.note_ignored_settings("msqrob_agg", cfg)) %in% 0,
  "test_ridge on the one method that reads it is not noted")
report(length(f.note_ignored_settings("msqrob", cfg)) %in% 1,
  "and is noted on the aggregate-then-fit path, which does not")

## config$test_moderate is read by the least squares prolfqua path and by nothing else, so
##   a FALSE set on any other method used to do nothing and say nothing, while the same
##   mistake with config$test_ridge drew this note. The mixed prolfqua path makes a note of
##   its own about the same key, but only when it is TRUE, so the two are complementary
##   rather than duplicates: 1/test_prolfqua_lmer.R covers that half:

cfg <- base_cfg()
cfg$test_moderate <- FALSE
report(length(f.note_ignored_settings("prolfqua", cfg)) %in% 0,
  "test_moderate on the one method that reads it is not noted")
got <- f.note_ignored_settings("trend", cfg)
report(length(got) %in% 1 && grepl("config$test_moderate=FALSE", got, fixed=TRUE),
  "but a FALSE is noted on a limma path, which moderates whatever the key says")
report(length(f.note_ignored_settings("lm", cfg)) %in% 1,
  "and on lm, which does not moderate at all")
report(length(f.note_ignored_settings("prolfqua_lmer", cfg)) %in% 1,
  "and on the mixed path, whose denominator is not one variance to shrink")
report(grepl("read by test_method prolfqua", got, fixed=TRUE),
  "the note naming the method that does read it")

cfg <- base_cfg()
cfg$test_moderate <- TRUE
report(length(f.note_ignored_settings("trend", cfg)) %in% 0,
  "while an explicit TRUE is not noted, that being what new_config() ships")

cfg <- base_cfg()
cfg$test_random_obs <- FALSE
cfg$test_ridge <- TRUE
m0 <- mark()
got <- f.note_ignored_settings("voom", cfg)
report(length(got) %in% 2, "two ignored settings are reported together")
report(logged("those settings have", m0), "in one note that reads as a plural")

cfg <- base_cfg()
cfg$test_ridge <- NULL
cfg$test_random_obs <- NULL
cfg$test_trend <- NULL
report(length(f.note_ignored_settings("trend", cfg)) %in% 0,
  "and a config that never mentioned the settings says nothing")

## the note is keyed on the method test() resolved rather than on config$test_method,
##   test(method=) overriding it, and this is where it is said:

set.seed(101)
exprs <- sim1(n_obs=12, n_feats=30)$mat
exprs <- log2(exprs + 1)

state_t <- list(
  expression=exprs,
  features=data.frame(feature_id=rownames(exprs), stringsAsFactors=FALSE),
  samples=data.frame(observation_id=colnames(exprs),
    grp=c(rep("ctl", 6), rep("trt", 6)), stringsAsFactors=FALSE)
)

## config$test_ridge is left out of this fixture on purpose: config$test_method is
##   prolfqua_lmer here, and the penalty is refused above on a two level factor that keeps
##   its intercept, which would stop the run in check_config() before the note. Any of the
##   three settings shows the same thing:

cfg <- base_cfg()
cfg$feat_id_col <- cfg$gene_id_col <- "feature_id"
cfg$obs_id_col <- cfg$sample_id_col <- "observation_id"
cfg$frm <- ~grp
cfg$test_term <- "grp"
cfg$reference_levels <- c(grp="ctl")
cfg$test_method <- "prolfqua_lmer"       ## a method that does read test_random_obs
cfg$test_random_obs <- FALSE
cfg$is_log_transformed <- TRUE

out_t <- try(suppressMessages(initialize(state_t, cfg, minimal=TRUE)), silent=TRUE)
report(!inherits(out_t, "try-error"), "the fixture for the run below initializes")

m0 <- mark()
res_t <- try(suppressMessages(test(out_t$state, out_t$config, method="trend")),
  silent=TRUE)

report(!inherits(res_t, "try-error") && nrow(res_t$standard) %in% nrow(exprs),
  "test(method=) overriding config$test_method still runs")
report(logged("config$test_random_obs=FALSE", m0),
  "and notes the setting against the method it resolved, not the one configured")
report(logged("NOTE: test:", m0), "the note saying which function it comes from")
report(log_count("does not consult", m0) %in% 1,
  "said once for the run rather than once per step of it")

###############################################################################
section("config$run_order against the steps run() can walk")

## run() fetches each step with get(), so a name that is not a step fails as an "object
##   not found" at the point that step would have run, with every step before it already
##   done and its output already written. The set is fixed, so the configuration settles it:

cfg <- base_cfg()
report(isTRUE(check_config(cfg)), "the default run_order passes")

report(setequal(f.run_order_steps(),
  c("normalize", "combine_replicates", "combine_features", "filter", "impute")),
  "five steps are on offer")
report(all(new_config()$run_order %in% f.run_order_steps()),
  "and the default run_order is drawn from them")

ok <- TRUE
for(step in f.run_order_steps()) {
  cfg <- base_cfg()
  cfg$run_order <- step
  if(!isTRUE(try(check_config(cfg), silent=TRUE))) ok <- FALSE
  if(!is.function(try(get(step), silent=TRUE))) ok <- FALSE
}
report(ok, "each one passes on its own, and each names a function that exists")

cfg <- base_cfg()
cfg$run_order <- c("normalize", "combine_featurez", "filter")
m0 <- mark()
report(threw(check_config(cfg)), "a step that is not one of them is refused")
report(logged("combine_featurez", m0), "and the message names the offending step")
report(logged("normalize", m0), "and lists the ones that are allowed")

cfg <- base_cfg()
cfg$run_order <- c("normalize", NA_character_)
report(threw(check_config(cfg)),
  "an NA step is refused too, absence being how a step is left out")

cfg <- base_cfg()
cfg$run_order <- c("normalize", "")
report(threw(check_config(cfg)), "as is an empty one")

cfg <- base_cfg()
cfg$run_order <- "test"
report(threw(check_config(cfg)),
  "and a function of the package that is not a pipeline step, get() finding it anyway")

cfg <- base_cfg()
cfg$run_order <- character(0)
report(isTRUE(check_config(cfg)),
  "no steps at all is a workflow, being load_data() and then test()")

cfg <- base_cfg()
cfg$run_order <- c("normalize", "combine_features", "normalize")
m0 <- mark()
report(isTRUE(check_config(cfg)),
  "a step named twice is not refused, running one again being a thing to mean")
report(!logged("more than once", m0),
  "and check_config() says nothing about it, run() saying it once instead")

###############################################################################
section("run() and a step named twice")

## the warning belongs to run(), which walks config$run_order once, rather than to
##   check_config(), which every step of it calls. The demo data shipped with the package
##   is enough to run the whole workflow, this being the only place the suite does:

## from the source tree rather than system.file(), the suite sourcing R/ rather than
##   loading an installed package:

extdata <- file.path(dirname(r_dir), "inst", "extdata")
if(!dir.exists(extdata)) extdata <- system.file("extdata", package="h0testr")

cfg_r <- base_cfg()                      ## save_state FALSE and dir_out temporary there
cfg_r$dir_in <- extdata
cfg_r$feature_file_in <- "features.tsv"
cfg_r$sample_file_in <- "samples.tsv"
cfg_r$data_file_in <- "expression.tsv"
cfg_r$feat_id_col <- cfg_r$gene_id_col <- "feature_id"
cfg_r$obs_id_col <- cfg_r$sample_id_col <- "observation_id"
cfg_r$frm <- ~condition
cfg_r$test_term <- "condition"
cfg_r$test_method <- "trend"
cfg_r$reference_levels <- c(condition="placebo")
cfg_r$n_features_min <- 10
cfg_r$run_order <- c("normalize", "combine_replicates", "filter", "impute")

report(dir.exists(cfg_r$dir_in) && file.exists(file.path(cfg_r$dir_in, "expression.tsv")),
  "the package ships the demo data this section needs")

m0 <- mark()
res_r <- try(suppressMessages(run(cfg_r)), silent=TRUE)
report(!inherits(res_r, "try-error") && nrow(res_r$standard) > 0,
  "run() completes on a run_order that names each step once")
report(!logged("more than once", m0), "with nothing said about repeats")

cfg_r2 <- cfg_r
cfg_r2$run_order <- c("normalize", "combine_replicates", "filter", "impute", "filter")
m0 <- mark()
res_r2 <- try(suppressMessages(run(cfg_r2)), silent=TRUE)
report(!inherits(res_r2, "try-error"), "a repeated step runs rather than being refused")
report(log_count("more than once", m0) %in% 1,
  "and is warned about once, run() walking the order once")
report(logged("overwrites the files", m0),
  "the warning saying what the repeat costs, which is the saved state of the earlier one")

###############################################################################
section("test_method 'none', which skips the test step")

## "none" was reachable by an explicit test(method="none") and by no configuration:
##   check_config() allowed test_method only in c(test_methods(), ""), and test_methods()
##   does not name it. It is now a legal key value, so that run() can be used as a
##   preprocessing workflow, and it is still not one of test_methods(), which names the
##   engines and is looped over. check_config() is the one place that says what the key
##   may hold:

report(!("none" %in% test_methods()),
  "test_methods() names the engines and not this")

cfg <- base_cfg()
cfg$test_method <- "none"
ok_n <- try(check_config(cfg), silent=TRUE)     ## try(), so a refusal fails one assertion
report(isTRUE(ok_n), "but check_config() accepts it as a key value")

cfg$test_method <- "nonesuch"
m0 <- mark()
report(threw(check_config(cfg)), "while a name that is neither an engine nor it is refused")
report(logged("meaning skip the test step", m0),
  "and the refusal says what 'none' is for, that being where a reader would look")

## the value returns before anything is read from state. It used to be the last branch of
##   test()'s dispatch chain, by which point a design had been built and rank checked, so
##   a state whose test term is not estimable failed at the step it was told to skip --
##   which is exactly the state someone would want to stop before testing:

state_n <- list(
  expression=matrix(1:24 + 0, nrow=4,
    dimnames=list(paste0("f", 1:4), paste0("o", 1:6))),
  features=data.frame(feature_id=paste0("f", 1:4), stringsAsFactors=FALSE),
  samples=data.frame(observation_id=paste0("o", 1:6),
    grp=rep("ctl", 6), stringsAsFactors=FALSE)      ## one level: nothing to test
)

cfg_n <- base_cfg()
cfg_n$feat_id_col <- cfg_n$gene_id_col <- "feature_id"
cfg_n$obs_id_col <- cfg_n$sample_id_col <- "observation_id"
cfg_n$frm <- ~grp
cfg_n$test_term <- "grp"
cfg_n$test_method <- "none"

m0 <- mark()
res_n <- try(suppressMessages(test(state_n, cfg_n)), silent=TRUE)
report(!inherits(res_n, "try-error") && is.null(res_n),
  "test() returns NULL for it, on a state whose test term cannot be estimated")
report(logged("skipping testing", m0), "saying so in the log")

cfg_n2 <- cfg_n
cfg_n2$test_method <- "trend"
report(threw(suppressMessages(test(state_n, cfg_n2))),
  "the same state under a real method fails, which is what 'none' steps around")

res_a <- try(suppressMessages(test(state_n, cfg_n2, method="none")), silent=TRUE)
report(!inherits(res_a, "try-error") && is.null(res_a),
  "the argument reaches it too, overriding a config that names an engine")

## and run() with it: the whole workflow, then no test:

cfg_rn <- cfg_r
cfg_rn$test_method <- "none"
m0 <- mark()
res_rn <- try(suppressMessages(run(cfg_rn)), silent=TRUE)
report(!inherits(res_rn, "try-error"), "run() completes with test_method 'none'")

if(!inherits(res_rn, "try-error")) {
  report(is.null(res_rn$original) && is.null(res_rn$standard) && is.null(res_rn$fit),
    "with no result in any of the three result slots")
  report(!is.null(res_rn$state) && nrow(res_rn$state$expression) > 0,
    "and the processed state in place, which is what it was run for")
  report(identical(dim(res_rn$state$expression), dim(res_r$state$expression)),
    "the same state the tested run produced, the steps before the test being the same")
  report(logged("skipping testing", m0), "and the skip recorded in the log")
}

## tune() refuses it, and refuses a typo for the same reason: the sweep assigns each name
##   to config$test_method in turn, and a sweep over not testing has nothing to compare:

m0 <- mark()
report(threw(suppressMessages(tune(cfg_r, test_methods=c("trend", "none")))),
  "tune() refuses 'none' among its methods")
report(logged("unexpected test_methods:", m0), "naming what it would not accept")
report(logged("nothing to compare", m0), "and why a sweep cannot use it")
report(threw(suppressMessages(tune(cfg_r, test_methods="trrend"))),
  "and refuses a typo, which used to fail only once that combination had been computed")

###############################################################################
section("the checks that were there before")

cfg <- base_cfg()
cfg$contrast <- "grptrt"
m0 <- mark()
report(threw(check_config(cfg)) && logged("one run tests one hypothesis", m0),
  "config$contrast and config$test_term both set is still refused")

cfg <- base_cfg()
cfg$estimability <- "nope"
report(threw(check_config(cfg)), "an unexpected estimability is still refused")

cfg <- base_cfg()
cfg$impute_quantile <- "0.01"
report(threw(check_config(cfg)), "and a value of the wrong type still is too")

cfg <- base_cfg()
cfg$no_such_param <- TRUE
report(threw(check_config(cfg)), "as is an unrecognized parameter name")

report(isTRUE(check_config(list())), "and an empty config is still accepted")

###############################################################################
section("NA and NaN in the numeric parameters")

## is.numeric() is TRUE for NA_real_ and for NaN, so both used to pass the class and length
##   test of the scalar_counts and scalar_props loops and reach the comparison after it,
##   which reported them as R's own "missing value where TRUE/FALSE needed" with nothing to
##   say which parameter it was about. Measured before the check was added, and the reason
##   it was added: h0testr::normalize_loess() and h0testr::normalize_mscoreutils() started
##   calling check_config(), which put config$normalization_span and
##   config$normalization_quantile on that path:

for(nom in c("n_features_min", "impute_k", "impute_n_pts")) {      ## scalar_counts
  for(v in list(NA_real_, NaN)) {
    cfg <- base_cfg()
    cfg[[nom]] <- v
    i <- mark()
    report(threw(check_config(cfg)) && logged("param is NA or NaN", i) && logged(nom, i),
      paste0("a count parameter of ", v, " is refused, naming ", nom))
    report(!logged("missing value where", i), "  and R's own message does not leak out")
  }
}

for(nom in c("normalization_quantile", "impute_quantile", "normalization_span")) {
  for(v in list(NA_real_, NaN)) {                                  ## scalar_props
    cfg <- base_cfg()
    cfg[[nom]] <- v
    i <- mark()
    report(threw(check_config(cfg)) && logged("param is NA or NaN", i) && logged(nom, i),
      paste0("a proportion parameter of ", v, " is refused, naming ", nom))
    report(!logged("missing value where", i), "  and R's own message does not leak out")
  }
}

## a logical NA was already caught by the class test, and says so:

cfg <- base_cfg()
cfg$normalization_quantile <- NA
i <- mark()
report(threw(check_config(cfg)) && logged("param not scalar proportion", i),
  "a logical NA is still refused by the class test, as before")

## and the values on either side of the checks are untouched:

cfg <- base_cfg()
cfg$normalization_quantile <- 0.5
cfg$n_features_min <- 3
report(isTRUE(check_config(cfg)), "a valid count and proportion are still accepted")

for(v in list(-0.5, 1.5)) {
  cfg <- base_cfg()
  cfg$normalization_quantile <- v
  i <- mark()
  report(threw(check_config(cfg)) && logged("param not between 0 and 1", i),
    paste0("a proportion of ", v, " is still refused by the range check"))
}

###############################################################################
section("Inf, and NA in the parameters the NA check missed")

## an infinite count passed every check it was given: Inf == round(Inf), and Inf is not < 0,
##   so config$n_features_min=Inf was accepted outright. Measured before the check went in,
##   which is the only reason it is known to have been reachable:

for(nom in c("n_features_min", "impute_k", "impute_n_pts")) {
  for(v in list(Inf, -Inf)) {
    cfg <- base_cfg()
    cfg[[nom]] <- v
    i <- mark()
    report(threw(check_config(cfg)) && logged("param not finite", i) && logged(nom, i),
      paste0("a count parameter of ", v, " is refused as not finite, naming ", nom))
    report(!logged("missing value where", i), "  and R's own message does not leak out")
  }
}

## config$impute_scale had neither check: NA_real_ and NaN reached the comparison, and Inf
##   was accepted:

for(v in list(NA_real_, NaN)) {
  cfg <- base_cfg()
  cfg$impute_scale <- v
  i <- mark()
  report(threw(check_config(cfg)) && logged("param is NA or NaN", i) &&
      logged("impute_scale", i),
    paste0("an impute_scale of ", v, " is refused, naming the parameter"))
  report(!logged("missing value where", i), "  and R's own message does not leak out")
}

for(v in list(Inf, -Inf)) {
  cfg <- base_cfg()
  cfg$impute_scale <- v
  i <- mark()
  report(threw(check_config(cfg)) && logged("param not finite", i) &&
      logged("impute_scale", i),
    paste0("an impute_scale of ", v, " is refused as not finite"))
}

## config$probs is a vector, so one NA anywhere in it made any(x < 0) itself NA:

for(v in list(c(0.5, NA_real_), c(NA_real_, 0.5), c(0.5, NaN))) {
  cfg <- base_cfg()
  cfg$probs <- v
  i <- mark()
  report(threw(check_config(cfg)) && logged("param has an NA or NaN", i) &&
      logged("probs", i),
    paste0("probs of c(", paste(v, collapse=", "), ") is refused, naming the parameter"))
  report(!logged("missing value where", i), "  and R's own message does not leak out")
}

## an infinite element of probs is left to the range check, which names it as what it is:

for(v in list(c(0.5, Inf), c(0.5, -Inf))) {
  cfg <- base_cfg()
  cfg$probs <- v
  i <- mark()
  report(threw(check_config(cfg)) && logged("proportions out of range", i),
    paste0("probs of c(", paste(v, collapse=", "),
      ") is refused as out of range, not as not finite"))
}

## and the proportions loop is deliberately not given a finiteness check either, for the
##   same reason: "not between 0 and 1" says more about a proportion than "not finite":

for(v in list(Inf, -Inf)) {
  cfg <- base_cfg()
  cfg$normalization_quantile <- v
  i <- mark()
  report(threw(check_config(cfg)) && logged("param not between 0 and 1", i) &&
      !logged("param not finite", i),
    paste0("a normalization_quantile of ", v, " is still refused as out of range"))
}

## the valid values on either side of all of these are untouched, and the checks that come
##   after the new ones still fire:

cfg <- base_cfg()
cfg$n_features_min <- 3
cfg$impute_k <- 5
cfg$impute_scale <- 2
cfg$probs <- c(0.25, 0.75)
report(isTRUE(check_config(cfg)),
  "valid counts, a valid impute_scale and a valid probs vector are still accepted")

cfg <- base_cfg()
cfg$n_features_min <- 3.7
i <- mark()
report(threw(check_config(cfg)) && logged("param not an integer", i),
  "a non-integer count is still refused by the check after the new ones")

cfg <- base_cfg()
cfg$impute_scale <- -2
i <- mark()
report(threw(check_config(cfg)) && logged("param not non-negative", i),
  "and a negative impute_scale by the check after those")

###############################################################################

cat("\n## passes:", n_pass, "; failures:", n_fail, "; elapsed:",
  round(as.numeric(difftime(Sys.time(), t0, units="secs"))), "s\n")
cat("## log file:", log_file, "\n")

quit(save="no", status=if(n_fail > 0) 1 else 0)
