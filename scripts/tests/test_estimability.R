## Tests for h0testr::filter_features_by_estimability(): the degrees of freedom
##   it computes per feature, the three nested config$estimability levels, the
##   df_resid_min screen, and its placement within h0testr::filter().

usage <- function(msg=NULL) {
  if(!is.null(msg)) cat("ERROR:", msg, "\n\n", file=stderr())
  cat(
    "Test h0testr::filter_features_by_estimability(): df_test and df_resid per",
    "feature, the nested config$estimability levels ('test', 'term', 'full'),",
    "the df_resid_min screen, the drop-reason accounting, and the call from",
    "h0testr::filter().",
    "",
    "Usage: Rscript test_estimability.R <r_dir>",
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
    "  Rscript test_estimability.R C:/path/to/h0test/h0testr/R",
    "  Rscript test_estimability.R ../../h0test/h0testr/R",
    "  Rscript test_estimability.R C:/path/to/h0testr/R > test_estimability.out 2>&1",
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

## returns TRUE iff expr threw an error:
threw <- function(expr) inherits(try(expr, silent=TRUE), "try-error")

###############################################################################
## Twelve observations: grp has three levels, four observations each; sex has
##   two levels, alternating. So with frm ~grp+sex the design has four columns
##   (intercept, grpb, grpc, sexM) and df_intend for grp is 2.

nobs <- 12
samps <- data.frame(
  obs=paste0("o", 1:nobs),
  grp=rep(c("a", "b", "c"), each=4),
  sex=rep(c("F", "M"), 6)
)

## one row per case, each chosen for the (df_test, df_resid, df_deficit) it
##   produces; the observations each feature was measured in:

obs_seen <- list(
  f1_complete      = 1:12,             ## df_test 2, df_resid 8, deficit 0
  f2_one_group     = 1:4,              ## df_test 0: only grp 'a' measured
  f3_two_groups    = 1:8,              ## df_test 1 < 2: grp 'c' missing
  f4_one_sex       = c(1, 3, 5, 7, 9, 11),  ## df_test 2, but sex collapsed
  f5_no_resid_df   = c(1, 2, 5, 9),    ## df_test 2, df_resid 0
  f6_resid_df_2    = c(1, 2, 5, 6, 9, 10),  ## df_test 2, df_resid 2
  f7_all_missing   = integer(0)        ## df_test 0, df_resid 0
)

set.seed(101)
exprs <- matrix(NA_real_, nrow=length(obs_seen), ncol=nobs,
  dimnames=list(names(obs_seen), samps$obs))
for(nom in names(obs_seen)) {
  i <- obs_seen[[nom]]
  if(length(i) > 0) exprs[nom, i] <- runif(length(i), min=10, max=20)
}

feats <- data.frame(pep=rownames(exprs), gene=rownames(exprs))

cfg0 <- list(
  obs_id_col="obs", sample_id_col="obs", feat_id_col="pep", gene_id_col="gene",
  frm=~grp + sex, test_term="grp",
  reference_levels=c(grp="a", sex="F"),
  estimability="test", df_resid_min=2,
  df_test_col="df_test", df_resid_col="df_resid",
  log_file=log_file
)

out <- initialize(list(expression=exprs, features=feats, samples=samps), cfg0,
  minimal=TRUE)
state0 <- out$state
config0 <- out$config

report(identical(levels(state0$samples$grp), c("a", "b", "c")),
  "initialize() ordered grp levels as declared")

###############################################################################
section("degrees of freedom per feature")

state <- filter_features_by_estimability(state0, config0)

## kept, in original order: f1, f3, f4, f6

report(identical(state$features$pep,
  c("f1_complete", "f3_two_groups", "f4_one_sex", "f6_resid_df_2")),
  "estimability 'test' keeps exactly the four testable features")
report(identical(as.integer(state$features$df_test), c(2L, 1L, 2L, 2L)),
  "df_test correct for survivors")
report(identical(as.integer(state$features$df_resid), c(8L, 5L, 3L, 2L)),
  "df_resid correct for survivors")
report(nrow(state$expression) == nrow(state$features),
  "expression and features stay aligned")
report(identical(rownames(state$expression), state$features$pep),
  "expression rownames match kept features")

## df columns describe survivors only, following add_filter_stats():
report(!any(c("df_test", "df_resid") %in% names(state0$features)),
  "input features left unmodified")

###############################################################################
section("nested estimability levels")

cfg <- config0
cfg$estimability <- "term"
state <- filter_features_by_estimability(state0, cfg)
report(identical(state$features$pep,
  c("f1_complete", "f4_one_sex", "f6_resid_df_2")),
  "estimability 'term' additionally drops the partially estimable feature")

cfg$estimability <- "full"
state <- filter_features_by_estimability(state0, cfg)
report(identical(state$features$pep, c("f1_complete", "f6_resid_df_2")),
  "estimability 'full' additionally drops the collapsed covariate feature")

## the collapsed-covariate feature has the full test df, so only the full-rank
##   requirement can catch it:
report(TRUE == {
  cfg2 <- config0
  cfg2$estimability <- "test"
  s <- filter_features_by_estimability(state0, cfg2)
  "f4_one_sex" %in% s$features$pep &&
    s$features$df_test[s$features$pep %in% "f4_one_sex"] == 2
}, "collapsed covariate keeps full df_test, so df_test cannot reveal it")

###############################################################################
section("df_resid_min")

cfg <- config0
cfg$df_resid_min <- 3
state <- filter_features_by_estimability(state0, cfg)
report(!("f6_resid_df_2" %in% state$features$pep),
  "df_resid_min=3 drops the feature with df_resid 2")
report("f4_one_sex" %in% state$features$pep,
  "df_resid_min=3 keeps the feature with df_resid 3")

cfg$df_resid_min <- 0
state <- filter_features_by_estimability(state0, cfg)
report("f5_no_resid_df" %in% state$features$pep,
  "df_resid_min=0 keeps the feature with df_resid 0")

###############################################################################
section("arguments override config")

state <- filter_features_by_estimability(state0, config0, estimability="full")
report(identical(state$features$pep, c("f1_complete", "f6_resid_df_2")),
  "estimability argument overrides config$estimability")

state <- filter_features_by_estimability(state0, config0, df_resid_min=3)
report(!("f6_resid_df_2" %in% state$features$pep),
  "df_resid_min argument overrides config$df_resid_min")

###############################################################################
section("drop-reason accounting in the log")

n0 <- length(readLines(log_file))
state <- filter_features_by_estimability(state0, config0)
txt <- readLines(log_file)
txt <- txt[(n0 + 1):length(txt)]

report(any(grepl("filtering 3 features by estimability; keeping 4", txt,
  fixed=TRUE)), "log reports the number dropped and kept")
report(any(grepl("dropped, test_term_not_estimable : 2", txt, fixed=TRUE)),
  "log attributes 2 drops to an inestimable test_term")
report(any(grepl("dropped, too_few_residual_df : 1", txt, fixed=TRUE)),
  "log attributes 1 drop to too few residual df")
report(any(grepl("test_term only partly estimable: 3", txt, fixed=TRUE)),
  "log counts partly estimable features whether or not they were dropped")
## f2 and f7 lose grp columns, f3 loses grpc, f4 loses sexM; f5 and f6 are full
##   rank despite being sparse:
report(any(grepl("requested model is not full rank: 4", txt, fixed=TRUE)),
  "log counts rank deficient features whether or not they were dropped")

###############################################################################
section("bad settings")

cfg <- config0
cfg$estimability <- "bogus"
report(threw(filter_features_by_estimability(state0, cfg)),
  "unrecognized estimability rejected by check_config-style guard")

cfg <- config0
cfg$df_resid_min <- -1
report(threw(filter_features_by_estimability(state0, cfg)),
  "negative df_resid_min rejected")

cfg <- config0
cfg$frm <- NULL
report(threw(filter_features_by_estimability(state0, cfg)),
  "missing config$frm rejected")

cfg <- config0
cfg$test_term <- "not_a_term"
report(threw(filter_features_by_estimability(state0, cfg)),
  "test_term absent from frm rejected")

###############################################################################
section("interaction test_term and single-variable formula")

## test_term naming a variable that also appears in an interaction drops both
##   terms, so df_intend is larger:

cfg <- config0
cfg$frm <- ~grp + sex + grp:sex
cfg$test_term <- "grp:sex"
state <- try(filter_features_by_estimability(state0, cfg), silent=TRUE)
report(!inherits(state, "try-error"),
  "interaction test_term handled")
if(!inherits(state, "try-error")) {
  report(!("f4_one_sex" %in% state$features$pep),
    "feature with one sex cannot support a grp:sex interaction")
}

cfg <- config0
cfg$frm <- ~grp
cfg$reference_levels <- c(grp="a")
state <- try(filter_features_by_estimability(state0, cfg), silent=TRUE)
report(!inherits(state, "try-error"), "single-variable formula handled")
if(!inherits(state, "try-error")) {
  report("f4_one_sex" %in% state$features$pep,
    "dropping sex from frm makes the collapsed covariate irrelevant")
}

###############################################################################
section("reduced model is the full design minus columns")

## Both of these once failed, because test_lm() rebuilt its reduced model from
##   formula text: model.matrix() re-codes the remaining factors to full rank,
##   restoring the span that was meant to be removed. Column subsetting keeps the
##   reduced model a strict sub-model, and makes test_lm() agree with
##   filter_features_by_estimability() on how many df the test has.

cfg <- config0
cfg$save_state <- FALSE
cfg$feat_col <- "pep"
cfg$obs_col <- "obs"

## testing the intercept against a factor: reduced model is [grpb, grpc], rank
##   2, so this is a 1 df test of whether the reference level's mean is zero;
##   rebuilding '~0 + grp' instead gives rank 3 and nothing to test:

## test_lm() is run on the filtered state throughout this section, which is what
##   the estimability filter is there for. It used to be the only way to run it at
##   all here: stats::lm() cannot fit a feature with no observations, and that
##   error propagated out of the apply() in test_lm() and cost every other feature
##   its result. test_lm() now guards on the same two ranks the filter screens on
##   and leaves such a feature NA; see test_lm_ftest.R, which checks that the
##   features it leaves NA are exactly the ones this filter drops:

cfg1 <- cfg
cfg1$test_term <- "1"
d <- try(f.design_test_cols(state0, cfg1), silent=TRUE)
report(!inherits(d, "try-error") && d$df_intend %in% 1,
  "test_term '1' with a factor in frm is a 1 df test, not a vacuous one")

res <- try(test_lm(filter_features_by_estimability(state0, cfg1), cfg1),
  silent=TRUE)
report(!inherits(res, "try-error"), "test_lm() handles test_term '1'")
if(!inherits(res, "try-error")) {
  report(all(res$hits$pval < 0.05),
    "test_term '1' rejects a zero mean for log-scale abundances")
  report(!any(is.na(res$hits$pval)), "test_term '1' yields no NA p-values")
}

## no intercept, and test_term names the only term: reduced model has no
##   parameters at all, so the test is against zero rather than against a common
##   mean. Runs, agrees with the filter, and warns:

cfg2 <- cfg
cfg2$frm <- ~0 + grp
cfg2$reference_levels <- c(grp="a")
cfg2$test_term <- "grp"

n0 <- length(readLines(log_file))
d <- try(f.design_test_cols(state0, cfg2), silent=TRUE)
report(!inherits(d, "try-error") && d$df_intend %in% 3,
  "frm ~0 + grp testing grp is a 3 df test against zero")

txt <- readLines(log_file)
txt <- txt[(n0 + 1):length(txt)]
report(any(grepl("reduced model with no parameters", txt, fixed=TRUE)),
  "empty reduced model warns instead of failing")

res <- try(test_lm(filter_features_by_estimability(state0, cfg2), cfg2),
  silent=TRUE)
report(!inherits(res, "try-error"), "test_lm() handles frm without an intercept")
if(!inherits(res, "try-error")) {
  report(all(c("grpa", "grpb", "grpc") %in% names(res$hits)),
    "no-intercept design reports the coefficients actually fitted")
}

## test_lm() and the estimability filter must report the same df for the same
##   configuration, since one screens what the other tests:

cfg3 <- cfg
cfg3$frm <- ~0 + grp + sex
cfg3$test_term <- "grp"
d <- f.design_test_cols(state0, cfg3)
st <- filter_features_by_estimability(state0, cfg3, df_resid_min=0)
report(all(st$features$df_test[st$features$pep %in% "f1_complete"] %in%
  d$df_intend),
  "filter and test agree on df for a complete feature, no intercept")

## testing a variable inside an interaction tests the interaction too, so its
##   coefficients belong in the output:

cfg4 <- cfg
cfg4$frm <- ~grp * sex
cfg4$test_term <- "grp"
res <- try(test_lm(filter_features_by_estimability(state0, cfg4), cfg4),
  silent=TRUE)
report(!inherits(res, "try-error") &&
  any(grepl("sexM", setdiff(names(res$hits), c("pep", "gene")))),
  "coefficients of higher-order terms under test are reported")

###############################################################################
section("called from filter()")

cfg <- config0
cfg$n_samples_min <- 2
cfg$n_features_min <- 1
cfg$n_samples_expr_col <- "n_samps_expr"
cfg$median_raw_col <- "median_raw"
cfg$n_features_expr_col <- "n_feats_expr"
cfg$save_state <- FALSE
cfg$feat_col <- "pep"
cfg$obs_col <- "obs"

out <- filter(state0, cfg)
report(all(c("df_test", "df_resid") %in% names(out$state$features)),
  "filter() writes the df columns into state$features")
report(all(out$state$features$df_resid >= 2),
  "filter() applied the df_resid_min screen")

out2 <- filter(state0, cfg, filter_by_estimability=FALSE)
report(!any(c("df_test", "df_resid") %in% names(out2$state$features)),
  "filter_by_estimability=FALSE skips the filter")
report(nrow(out2$state$features) >= nrow(out$state$features),
  "skipping the estimability filter keeps at least as many features")

###############################################################################

cat("\n## passes:", n_pass, "; failures:", n_fail, "; elapsed:",
  round(as.numeric(difftime(Sys.time(), t0, units="secs"))), "s\n")
cat("## log file:", log_file, "\n")

quit(save="no", status=if(n_fail > 0) 1 else 0)
