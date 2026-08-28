## Tests for scale tracking: config$is_log_transformed as the single record of
##   what scale state$expression is on, config$log_from_raw as the narrower
##   claim that an exact 0 in the matrix cannot be a measurement, and the
##   f.check_state() tripwire that enforces it.

usage <- function(msg=NULL) {
  if(!is.null(msg)) cat("ERROR:", msg, "\n\n", file=stderr())
  cat(
    "Test scale tracking in h0testr: the config$is_log_transformed default and",
    "validation, init_state() filling it in, normalize() setting it and",
    "config$log_from_raw, the f.check_state() tripwire against zeros used as a",
    "missingness sentinel, and resolution of the is_log_transformed argument",
    "against the config value in impute() and the test functions.",
    "",
    "Usage: Rscript test_scale.R <r_dir>",
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
    "  Rscript test_scale.R C:/path/to/h0test/h0testr/R",
    "  Rscript test_scale.R ../../h0test/h0testr/R",
    "  Rscript test_scale.R C:/path/to/h0testr/R > test_scale.out 2>&1",
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

## TRUE iff the log holds the text, which is where f.err() puts its detail;
##   since= skips lines written before a mark(), so that a message left by an
##   earlier assertion is not read as evidence for a later one:

mark <- function() length(readLines(log_file))

logged <- function(pattern, since=0) {
  txt <- readLines(log_file)
  if(since > 0) {
    if(length(txt) <= since) return(FALSE)
    txt <- txt[(since + 1):length(txt)]
  }
  any(grepl(pattern, txt, fixed=TRUE))
}

###############################################################################
## shared data; sim1() draws strictly positive values with NA for dropouts:

set.seed(101)
exprs <- sim1(n_obs=12, n_feats=60, mcar_p=0.1)$mat
feats <- data.frame(pep=rownames(exprs), gene=rownames(exprs))
samps <- data.frame(obs=colnames(exprs), grp=rep(c("ctl", "trt"), 6))

mk_state <- function(e) list(expression=e, features=feats, samples=samps)

cfg0 <- list(
  obs_id_col="obs", sample_id_col="obs", feat_id_col="pep", gene_id_col="gene",
  frm=~grp, test_term="grp", reference_levels=c(grp="trt"), log_file=log_file,
  feat_col="pep", obs_col="obs", save_state=FALSE,
  n_samples_expr_col="n_samps_expr", median_raw_col="median_raw",
  n_features_expr_col="n_feats_expr"
)

###############################################################################
section("config defaults and validation")

cfg <- new_config()
report(identical(cfg$is_log_transformed, FALSE),
  "new_config() defaults is_log_transformed to FALSE")
report(identical(cfg$impute_floor_offset, -1),
  "new_config() defaults impute_floor_offset to -1")
report(is.null(cfg$zeros_to_na), "new_config() no longer has zeros_to_na")
report(isTRUE(check_config(cfg)), "check_config() accepts the new defaults")

## is_log_transformed answers a yes/no question about the data, so anything
##   that is not TRUE or FALSE would have to be guessed at:

for(bad in list("TRUE", 1, NA, c(TRUE, TRUE))) {
  cfg <- cfg0
  cfg$is_log_transformed <- bad
  report(threw(check_config(cfg)),
    paste("check_config() rejects is_log_transformed:",
      paste(as.character(bad), collapse=",")))
}

## impute_floor_offset is an offset below an observed value, so a positive one
##   would put the bound it derives above the value it is measured from:

cfg <- cfg0
cfg$impute_floor_offset <- 0.5
report(threw(check_config(cfg)),
  "check_config() rejects a positive impute_floor_offset")

cfg$impute_floor_offset <- -Inf
report(threw(check_config(cfg)),
  "check_config() rejects a non-finite impute_floor_offset")

cfg$impute_floor_offset <- "-1"
report(threw(check_config(cfg)),
  "check_config() rejects a non-numeric impute_floor_offset")

cfg$impute_floor_offset <- 0
report(isTRUE(check_config(cfg)), "check_config() accepts an offset of 0")

cfg$impute_floor_offset <- -2.5
report(isTRUE(check_config(cfg)), "check_config() accepts a negative offset")

###############################################################################
section("init_state() records the scale")

out <- init_state(mk_state(exprs), cfg0, minimal=TRUE)
report(identical(out$config$is_log_transformed, FALSE),
  "init_state() fills in is_log_transformed when unset")

cfg <- cfg0
cfg$is_log_transformed <- TRUE
out <- init_state(mk_state(exprs), cfg, minimal=TRUE)
report(identical(out$config$is_log_transformed, TRUE),
  "init_state() leaves a declared TRUE alone")

###############################################################################
section("normalize() sets the scale and the zero guarantee")

out <- init_state(mk_state(exprs), cfg0, minimal=TRUE)
st0 <- out$state
cf0 <- out$config

out <- normalize(st0, cf0, method="RLE")
report(identical(out$config$is_log_transformed, TRUE),
  "normalize(RLE) sets is_log_transformed")
report(identical(out$config$log_from_raw, TRUE),
  "normalize(RLE) claims the zero guarantee")
report(!any(out$state$expression %in% 0),
  "normalize(RLE) leaves no exact zeros to claim it about")

out <- normalize(st0, cf0, method="log2")
report(identical(out$config$is_log_transformed, TRUE),
  "normalize(log2) sets is_log_transformed")
report(identical(out$config$log_from_raw, TRUE),
  "normalize(log2) claims the zero guarantee")

## 'none' asks for no normalization, and says nothing about the scale, so it
##   leaves is_log_transformed as the caller declared it, either way. It leaves
##   log_from_raw alone for the same reason: it did not touch the matrix, so
##   whatever the guarantee was before it is still true after it. Setting FALSE
##   there made a step that changes nothing switch off f.check_state()'s zero
##   check for the rest of the run:

out <- normalize(st0, cf0, method="none")
report(identical(out$config$is_log_transformed, FALSE),
  "normalize(none) leaves raw input raw")
report(is.null(out$config$log_from_raw),
  "normalize(none) does not invent a zero guarantee")

cf1 <- cf0
cf1$is_log_transformed <- TRUE
out <- normalize(st0, cf1, method="none")
report(identical(out$config$is_log_transformed, TRUE),
  "normalize(none) leaves declared log input log")

cf1$log_from_raw <- TRUE
out <- normalize(st0, cf1, method="none")
report(identical(out$config$log_from_raw, TRUE),
  "normalize(none) leaves a claimed zero guarantee claimed")

cf1$log_from_raw <- FALSE
out <- normalize(st0, cf1, method="none")
report(identical(out$config$log_from_raw, FALSE),
  "normalize(none) leaves a withheld zero guarantee withheld")

## the consequence: after a real normalization, 'none' no longer disarms the
##   check that catches a missing value written in as a measurement:

out <- normalize(st0, cf0, method="RLE")
out <- normalize(out$state, out$config, method="none")
st2 <- out$state
st2$expression[1, 1] <- 0
report(threw(f.check_state(st2, out$config)),
  "a zero written in after normalize(none) is still caught")

report(threw(normalize(st0, cf1, method="RLE")),
  "normalize() refuses to transform declared log input again")

## vsn output is arsinh-scaled rather than log2: on a log-like scale, so not to
##   be transformed again, but with no guarantee about exact zeros:

if(requireNamespace("vsn", quietly=TRUE)) {
  out <- try(normalize(st0, cf0, method="vsn"), silent=TRUE)
  if(!inherits(out, "try-error")) {
    report(identical(out$config$is_log_transformed, TRUE),
      "normalize(vsn) sets is_log_transformed")
    report(identical(out$config$log_from_raw, FALSE),
      "normalize(vsn) claims no zero guarantee")
  } else {
    cat("SKIP: vsn present but normalize(vsn) failed on this toy data\n")
  }
} else {
  cat("SKIP: vsn not installed\n")
}

###############################################################################
section("f.check_state() tripwire on exact zeros")

out <- normalize(st0, cf0, method="RLE")
st_log <- out$state
cf_log <- out$config

report(!threw(f.check_state(st_log, cf_log)),
  "f.check_state() passes a matrix with no exact zeros")

## a 0 planted where log2(x + 1) of a measured value cannot land: exactly what
##   combining an all-missing replicate group with sum(na.rm=TRUE) used to write

st_bad <- st_log
st_bad$expression[5, 3] <- 0
report(threw(f.check_state(st_bad, cf_log)),
  "f.check_state() rejects a planted exact zero")
report(logged("exact zeros in state$expression"),
  "tripwire message names the problem")
report(logged(rownames(exprs)[5]), "tripwire message names the feature")
report(logged(colnames(exprs)[3]), "tripwire message names the observation")

## unset means normalize() made no such guarantee, so there is nothing to check:

cf_none <- cf_log
cf_none$log_from_raw <- NULL
report(!threw(f.check_state(st_bad, cf_none)),
  "tripwire is off when log_from_raw is unset")

cf_off <- cf_log
cf_off$log_from_raw <- FALSE
report(!threw(f.check_state(st_bad, cf_off)),
  "tripwire is off when log_from_raw is FALSE")

## NA is a missing value, not a zero, and is what the tripwire is protecting:

st_na <- st_log
st_na$expression[5, 3] <- NA
report(!threw(f.check_state(st_na, cf_log)),
  "tripwire does not mistake NA for a zero")

## and it fires through an ordinary step, not just when called directly. Every
##   workflow step checks the state it returns, so a sentinel written by one is
##   caught there; impute() also checks the state it is handed, which is what a
##   sentinel written upstream runs into first. Checked from both sides, since
##   the assertion would otherwise hold for the wrong reason if the step never
##   ran at all:

cf_step <- cf_log
cf_step$n_samples_min <- 2
cf_step$n_features_min <- 5
cf_step$impute_method <- "unif_sample_lod"
cf_step$impute_quantile <- 0

report(!threw(impute(st_log, cf_step)),
  "the step itself runs on a state with no planted zero")
m <- mark()
report(threw(impute(st_bad, cf_step)),
  "a workflow step handed a planted zero stops at the tripwire")
report(logged("exact zeros in state$expression", since=m),
  "the step stopped at the tripwire, not somewhere else")

## and on the way out: combine_replicates() with sum(na.rm=TRUE) is the writer
##   that used to record an all-missing replicate group as a measurement of 0,
##   so its own exit check is the one that would have caught it:

cf_rep <- cf_log
cf_rep$sample_id_col <- "obs"
report(!threw(combine_replicates(st_log, cf_rep, fn=sum)),
  "combine_replicates(sum) leaves no zeros for its exit check to find")

###############################################################################
section("f.is_log_transformed() resolution")

cf_bare <- cfg0                    ## no is_log_transformed
cf_raw <- cfg0
cf_raw$is_log_transformed <- FALSE
cf_lg <- cfg0
cf_lg$is_log_transformed <- TRUE

report(identical(f.is_log_transformed(TRUE, cf_bare, "t"), TRUE),
  "argument answers when the config does not")
report(identical(f.is_log_transformed(NULL, cf_lg, "t"), TRUE),
  "config answers when the argument is not given")
report(identical(f.is_log_transformed("", cf_raw, "t"), FALSE),
  "empty character argument counts as not given")
report(identical(f.is_log_transformed(TRUE, cf_lg, "t"), TRUE),
  "argument and config agreeing is fine")

report(threw(f.is_log_transformed(TRUE, cf_raw, "t")),
  "argument contradicting the config is an error")
report(logged("cannot both be right"),
  "contradiction message explains the conflict")

report(threw(f.is_log_transformed(NULL, cf_bare, "t")),
  "neither argument nor config set is an error")
report(logged("both"), "message says both are unset")

for(bad in list("TRUE", 1, NA, c(TRUE, TRUE))) {
  report(threw(f.is_log_transformed(bad, cf_bare, "t")),
    paste("rejects an argument of:", paste(as.character(bad), collapse=",")))
}

###############################################################################
section("impute() and test_h0() read the scale from config")

cf <- cf_log
cf$n_samples_min <- 2
cf$n_features_min <- 5
cf$impute_method <- "unif_sample_lod"
cf$impute_quantile <- 0

out <- filter_state(st_log, cf)
st_f <- out$state
cf_f <- out$config

## previously impute() required the argument from every caller; now the scale
##   normalize() recorded is enough:

out <- try(impute(st_f, cf_f), silent=TRUE)
report(!inherits(out, "try-error"), "impute() runs with no is_log_transformed")
if(!inherits(out, "try-error")) {
  report(!any(is.na(out$state$expression)), "impute() filled the matrix")
}

report(!threw(impute(st_f, cf_f, is_log_transformed=TRUE)),
  "impute() accepts an argument agreeing with the config")
report(threw(impute(st_f, cf_f, is_log_transformed=FALSE)),
  "impute() rejects an argument contradicting the config")

## the scale used to be guessed from config$normalization_method, which is not
##   a statement about scale at all; with no declaration anywhere it is now an
##   error rather than a guess:

cf_guess <- cf_f
cf_guess$is_log_transformed <- NULL
cf_guess$normalization_method <- "RLE"
report(threw(f.is_log_transformed(NULL, cf_guess, "test_prolfqua")),
  "normalization_method is no longer read as a scale declaration")

###############################################################################
section("combine_features() requires log scale")

## both aggregators fit an additive model, which holds for mass spectrometry
##   signal only after a log transform, so raw input is refused:

st_agg <- mk_state(exprs)
cf_agg <- cfg0
cf_agg$is_log_transformed <- FALSE

m1 <- mark()
report(threw(combine_features(st_agg, cf_agg, method="medianPolish")),
  "combine_features(medianPolish) refuses raw input")
report(threw(combine_features(st_agg, cf_agg, method="robustSummary")),
  "combine_features(robustSummary) refuses raw input")

## the message has to say how to get out of the situation, not just that the
##   situation is bad:

report(logged("run normalize() before combine_features()", since=m1),
  "the refusal points at normalize()")
report(logged("config$is_log_transformed <- TRUE", since=m1),
  "the refusal offers declaring an already transformed input")
report(logged("'none' to skip aggregation", since=m1),
  "the refusal offers method='none'")

## method "none" fits nothing, so it is exempt:

report(!threw(combine_features(st_agg, cf_agg, method="none")),
  "combine_features(none) is exempt from the requirement")

## and a declared log scale is accepted:

cf_agg$is_log_transformed <- TRUE
st_log_agg <- mk_state(log2(exprs + 1))
report(!threw(combine_features(st_log_agg, cf_agg, method="medianPolish")),
  "combine_features(medianPolish) accepts declared log scale input")

## an unset flag is not TRUE, so it is refused too, rather than assumed:

cf_agg$is_log_transformed <- NULL
report(threw(combine_features(st_log_agg, cf_agg, method="medianPolish")),
  "an unset is_log_transformed is refused rather than assumed")

###############################################################################
section("combine_features() refuses to rescale before aggregating")

## rescaling divided each feature by its own mean, a multiplicative operation,
##   and aggregation only runs on log data, where that has no reading: a mean
##   near zero explodes the feature and a negative one flips the sign of its
##   contrasts. Refused for both aggregators:

cf_rs <- cfg0
cf_rs$is_log_transformed <- TRUE
st_rs <- mk_state(log2(exprs + 1))

m5 <- mark()
report(threw(combine_features(st_rs, cf_rs, method="medianPolish", rescale=TRUE)),
  "combine_features(medianPolish) refuses rescale=TRUE")
report(threw(combine_features(st_rs, cf_rs, method="robustSummary", rescale=TRUE)),
  "combine_features(robustSummary) refuses rescale=TRUE")

report(logged("divided each feature by its mean", since=m5),
  "the refusal says what the rescaling did")
report(logged("config$feature_aggregation_scaled, has been removed", since=m5),
  "the refusal says the config key that used to set it is gone")
report(!("feature_aggregation_scaled" %in% names(new_config())),
  "new_config() no longer carries the key")

## method "none" aggregates nothing, so there is nothing to rescale before; the
##   setting is reported as ignored rather than made into an error:

m6 <- mark()
report(!threw(combine_features(st_rs, cf_rs, method="none", rescale=TRUE)),
  "rescale=TRUE with method='none' is not an error")
report(logged("no rescaling is done", since=m6),
  "method='none' reports the setting as ignored")

## the default path is untouched:

report(!threw(combine_features(st_rs, cf_rs, method="medianPolish")),
  "the default rescale=FALSE still aggregates")

## and the reason removal loses nothing: the log scale form of the rescaling,
##   subtracting the per-feature mean, is absorbed exactly by medianPolish()'s
##   own per-feature effect, so it moves each group's overall level and no
##   per-sample effect. Genes of three peptides each, since a group of one is
##   returned without fitting:

e_grp <- log2(exprs + 1)
feats_grp <- data.frame(pep=rownames(e_grp),
  gene=rep(paste0("g", 1:20), each=3))
st_a <- list(expression=e_grp, features=feats_grp, samples=samps)
st_b <- st_a
st_b$expression <- t(apply(e_grp, 1, function(v) v - mean(v, na.rm=TRUE)))

a <- combine_features(st_a, cf_rs, method="medianPolish")$state$expression
b <- combine_features(st_b, cf_rs, method="medianPolish")$state$expression

ctr <- function(m) m - apply(m, 1, mean, na.rm=TRUE)
report(isTRUE(all.equal(ctr(a), ctr(b))),
  "centering each feature first leaves every per-sample effect unchanged")
report(!isTRUE(all.equal(a, b)),
  "it does move each group's overall level, which the trended prior reads")

## whereas the division that was there changes the per-sample effects outright:

st_c <- st_a
st_c$expression <- t(apply(e_grp, 1, function(v) v / mean(v, na.rm=TRUE)))
d <- combine_features(st_c, cf_rs, method="medianPolish")$state$expression
report(max(abs(ctr(a) - ctr(d)), na.rm=TRUE) > 0.01,
  "dividing by the per-feature mean does not leave them unchanged")

###############################################################################
section("combine_features() no longer shifts a group to be positive")

## gene g1 sits far below the others in o1. medianPolish returns that sample at
##   about -3, and earlier code added 2 * abs(min) + 1 to the whole vector to
##   lift it positive, which moved every sample in the group, not just the one:

e_neg <- matrix(c(-3, 5, 6, 7,
                  -3, 5, 6, 7), nrow=2, byrow=TRUE,
  dimnames=list(c("p1", "p2"), paste0("o", 1:4)))
st_neg <- list(
  expression=e_neg,
  features=data.frame(pep=rownames(e_neg), gene=c("g1", "g1")),
  samples=data.frame(obs=colnames(e_neg), grp=c("ctl", "ctl", "trt", "trt"))
)
cf_neg <- list(
  feat_col="pep", feat_id_col="pep", gene_id_col="gene", obs_col="obs",
  obs_id_col="obs", sample_id_col="obs", frm=~grp, test_term="grp",
  is_log_transformed=TRUE, save_state=FALSE, log_file=log_file
)

out <- combine_features(st_neg, cf_neg, method="medianPolish")
report(min(out$state$expression, na.rm=TRUE) < 0,
  "medianPolish keeps a negative summary instead of shifting it positive")
report(isTRUE(all.equal(as.numeric(out$state$expression), c(-3, 5, 6, 7))),
  "the summary is the column effect itself, unshifted")

out <- combine_features(st_neg, cf_neg, method="robustSummary")
report(min(out$state$expression, na.rm=TRUE) < 0,
  "robustSummary keeps a negative summary instead of shifting it positive")

###############################################################################
section("f.pos_mat() is conditional on scale and reports rather than nudges")

## on the log scale a value at or below zero is an ordinary small intensity:

m_log <- matrix(c(-1, 0, 1, 2), nrow=2)
report(identical(f.pos_mat(m_log, cf_neg, TRUE, "probe"), m_log),
  "f.pos_mat() leaves non-positive log scale values alone")

## on the raw scale it cannot be a measurement, and used to be nudged upward
##   into a fabricated one; now it is reported, naming the imputer. cf_raw
##   declares the raw scale, so that the argument agrees with the config and the
##   positivity check is what fires, rather than the contradiction check:

cf_raw <- cf_neg
cf_raw$is_log_transformed <- FALSE

m2 <- mark()
report(threw(f.pos_mat(m_log, cf_raw, FALSE, "impute_probe")),
  "f.pos_mat() refuses non-positive raw scale values")
report(logged("impute_probe", since=m2),
  "the refusal names the imputer that produced the values")
report(logged("2 of 4 imputed values", since=m2),
  "the refusal counts the offending values")

## finiteness is not a property of the scale either, and an infinite value is
##   now reported rather than replaced. It used to be assigned a draw between the
##   largest finite value and twice it, which put an unmeasured feature at the
##   top of the scale when a missing value means it belongs at the bottom, sent
##   -Inf to the top as well, and moved the value downward whenever the largest
##   finite value was negative, the nudge being multiplicative:

m_inf <- matrix(c(1, 2, 3, Inf), nrow=2)
m3 <- mark()
report(threw(f.pos_mat(m_inf, cf_neg, TRUE, "impute_probe")),
  "f.pos_mat() refuses an infinite value on the log scale")
report(logged("impute_probe", since=m3),
  "the refusal names the imputer that produced the infinite value")
report(logged("1 of 4 imputed values", since=m3),
  "the refusal counts the infinite values")

m3 <- mark()
report(threw(f.pos_mat(m_inf, cf_raw, FALSE, "impute_probe")),
  "f.pos_mat() refuses an infinite value on the raw scale too")

## a -Inf on the raw scale is at or below zero, so the positivity check sees it
##   first and the finiteness check never does:

m_ninf <- matrix(c(1, 2, 3, -Inf), nrow=2)
m3 <- mark()
report(threw(f.pos_mat(m_ninf, cf_raw, FALSE, "impute_probe")),
  "f.pos_mat() refuses -Inf on the raw scale")
report(logged("at or below zero", since=m3),
  "the raw scale refusal of -Inf comes from the positivity check")

m3 <- mark()
report(threw(f.pos_mat(m_ninf, cf_neg, TRUE, "impute_probe")),
  "f.pos_mat() refuses -Inf on the log scale")
report(logged("infinite", since=m3),
  "the log scale refusal of -Inf comes from the finiteness check")

## impute_sample_lod() was the one imputer that never called f.pos_mat(), so an
##   observation with nothing measured came back as a column of Inf:

st_inf <- st_neg
st_inf$expression[, 1] <- NA
m3 <- mark()
report(threw(suppressWarnings(impute_sample_lod(st_inf, cf_neg))),
  "impute_sample_lod() refuses an observation with no measured value")
report(logged("impute_sample_lod", since=m3),
  "that refusal names impute_sample_lod()")

###############################################################################
section("impute() hands the resolved scale to imputers that lack the argument")

## impute_unif_sample_lod() has no is_log_transformed argument of its own, so it
##   reads config; impute() resolves the scale and records it there, and without
##   that an argument-only declaration would leave the imputer with nothing:

st_i <- mk_state(exprs)
cf_i <- cfg0
cf_i$is_log_transformed <- NULL
out_i <- try(impute(st_i, cf_i, method="unif_sample_lod",
  is_log_transformed=FALSE), silent=TRUE)
report(!inherits(out_i, "try-error"),
  "impute() propagates an argument-only scale to impute_unif_sample_lod()")

###############################################################################
section("the unif_ imputers draw from an explicit, scale aware floor")

i_na <- is.na(exprs)

## raw input: the floor is zero abundance, so imputed values land in (0, LOD],
##   where LOD is quantile 0 of the per-feature minima:

cf_f <- cfg0
cf_f$is_log_transformed <- FALSE
cf_f$impute_quantile <- 0

lod_raw <- min(apply(exprs, 1, min, na.rm=TRUE))
set.seed(11)
out_f <- impute_unif_global_lod(mk_state(exprs), cf_f)
v_f <- out_f$expression[i_na]

report(all(v_f > 0) && all(v_f <= lod_raw),
  "raw imputed values fall between zero and the global LOD")

## the offset is a log scale notion, since raw data have a real zero; same seed,
##   so an offset that was being consumed would show up as different draws:

cf_f2 <- cf_f
cf_f2$impute_floor_offset <- -5
set.seed(11)
out_f2 <- impute_unif_global_lod(mk_state(exprs), cf_f2)
report(identical(out_f2$expression, out_f$expression),
  "config$impute_floor_offset is ignored for raw input")

## log input: zero is not the bottom of the scale, so the floor goes below the
##   dimmest measured value instead. These data are wholly negative, where a
##   literal zero floor would have inverted the interval:

e_log <- log2(exprs) - 20
min_log <- min(e_log, na.rm=TRUE)

cf_g <- cfg0
cf_g$is_log_transformed <- TRUE
cf_g$impute_quantile <- 0
cf_g$impute_floor_offset <- -1

set.seed(11)
out_g <- impute_unif_global_lod(mk_state(e_log), cf_g)
v_g <- out_g$expression[i_na]

report(!any(is.na(out_g$expression)), "log input imputes without NaN")
report(all(v_g >= min_log - 1) && all(v_g <= min_log),
  "log imputed values fall between the offset floor and the dimmest measurement")

## the offset is what the floor is made of, so a more negative one reaches lower:

cf_h <- cf_g
cf_h$impute_floor_offset <- -5
set.seed(11)
out_h <- impute_unif_global_lod(mk_state(e_log), cf_h)
v_h <- out_h$expression[i_na]

report(min(v_h) < min(v_g) && all(v_h >= min_log - 5),
  "config$impute_floor_offset sets how far below the data the floor sits")

## an offset of zero leaves no width at all at quantile 0, which is an error
##   rather than a silent NaN:

cf_z <- cf_g
cf_z$impute_floor_offset <- 0
m3 <- mark()
report(threw(impute_unif_global_lod(mk_state(e_log), cf_z)),
  "an imputation interval with no width is refused")
report(logged("impute_floor_offset", since=m3),
  "the refusal names the setting that gives the interval its width")

## features whose minimum was at or below zero used to be dropped before the LOD
##   quantile was taken, which on the log scale discarded exactly the dimmest
##   features. Here the dimmest feature is the one below zero, so it sets the LOD:

e_mix <- log2(exprs)
lod_mix <- min(min(e_mix, na.rm=TRUE), 0) - 5
e_mix[1, ] <- e_mix[1, ] - min(e_mix[1, ], na.rm=TRUE) + lod_mix
set.seed(11)
out_m <- impute_unif_global_lod(mk_state(e_mix), cf_g)

report(lod_mix < 0 && all(out_m$expression[i_na] <= lod_mix),
  "a feature dimmer than zero is no longer dropped from the LOD estimate")

## sample LOD: one upper bound per observation, but a single floor, since the
##   bottom of the scale is a property of the data rather than of a column:

set.seed(11)
out_s <- impute_unif_sample_lod(mk_state(e_log), cf_g)
lod_col <- apply(e_log, 2, min, na.rm=TRUE)

ok_col <- sapply(seq_len(ncol(e_log)), function(k) {
  v <- out_s$expression[is.na(e_log[, k]), k]
  length(v) %in% 0 || (all(v <= lod_col[k]) && all(v >= min_log - 1))
})

report(all(ok_col),
  "each observation is imputed below its own LOD and above the shared floor")

## an observation measured nowhere has no LOD to impute against:

e_bad <- e_log
e_bad[, 3] <- NA
m4 <- mark()
report(threw(impute_unif_sample_lod(mk_state(e_bad), cf_g)),
  "an observation with no measured value is refused")
report(logged(colnames(e_log)[3], since=m4),
  "the refusal names the offending observation")

###############################################################################
section("step by step run needs no scale bookkeeping")

## the whole pipeline, with the scale carried in config from normalize() on and
##   never passed by hand:

cf <- cfg0
cf$normalization_method <- "RLE"
cf$impute_method <- "unif_sample_lod"
cf$impute_quantile <- 0
cf$test_method <- "trend"
cf$n_samples_min <- 2
cf$n_features_min <- 5
cf$feature_aggregation <- "none"

out <- init_state(mk_state(exprs), cf, minimal=TRUE)
out$state <- add_filter_stats(out$state, out$config)
out <- normalize(out$state, out$config)
out <- combine_replicates(out$state, out$config, fn=sum)
out <- filter_state(out$state, out$config)
out <- impute(out$state, out$config)
result <- try(test_h0(out$state, out$config), silent=TRUE)

report(!inherits(result, "try-error"), "pipeline completes with no scale argument")
if(!inherits(result, "try-error")) {
  report(nrow(result$standard) > 0, "pipeline produced a hit table")
  report(!any(is.na(result$standard$pval)), "no NA p-values")
}
report(identical(out$config$is_log_transformed, TRUE),
  "config carries the scale to the end of the pipeline")
report(identical(out$config$log_from_raw, TRUE),
  "config carries the zero guarantee to the end of the pipeline")

###############################################################################

cat("\n## log file:", log_file, "\n")
cat("\n## passes:", n_pass, "; failures:", n_fail, "; elapsed:",
  round(as.numeric(difftime(Sys.time(), t0, units="secs"))), "s\n")

quit(save="no", status=if(n_fail > 0) 1 else 0)
