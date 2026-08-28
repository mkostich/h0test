## Tests for the missingness contract established by h0testr::init_state():
##   NA is the only indicator of a missing value. Raw zeros are converted to NA
##   and negative raw values are an error, unless the input is declared to be
##   already transformed.

usage <- function(msg=NULL) {
  if(!is.null(msg)) cat("ERROR:", msg, "\n\n", file=stderr())
  cat(
    "Test the missingness contract established by h0testr::init_state():",
    "conversion of raw zeros to NA, rejection of negative raw values, the",
    "config$is_log_transformed escape, and the cross-check in",
    "h0testr::normalize() against a contradictory method argument.",
    "",
    "Usage: Rscript test_missingness.R <r_dir>",
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
    "  Rscript test_missingness.R C:/path/to/h0test/h0testr/R",
    "  Rscript test_missingness.R ../../h0test/h0testr/R",
    "  Rscript test_missingness.R C:/path/to/h0testr/R > test_missingness.out 2>&1",
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
## shared data; sim1() draws strictly positive values with NA for dropouts, so
##   zeros and negatives here are only the ones planted below:

set.seed(101)
exprs <- sim1(n_obs=6, n_feats=10, mcar_p=0.1)$mat
feats <- data.frame(pep=rownames(exprs), gene=rownames(exprs))
samps <- data.frame(obs=colnames(exprs), grp=rep(c("ctl", "trt"), 3))

mk_state <- function(e) list(expression=e, features=feats, samples=samps)

cfg0 <- list(
  obs_id_col="obs", sample_id_col="obs", feat_id_col="pep", gene_id_col="gene",
  frm=~grp, test_term="grp", reference_levels=c(grp="trt"), log_file=log_file
)

###############################################################################
section("raw input")

e <- exprs
e[1, 1] <- 0
e[2, 2] <- 0
n_na0 <- sum(is.na(e))

out <- init_state(mk_state(e), cfg0, minimal=TRUE)
report(sum(is.na(out$state$expression)) == n_na0 + 2,
  "raw zeros converted to NA")
report(!any(out$state$expression == 0, na.rm=TRUE), "no zeros left")
report(identical(out$config$is_log_transformed, FALSE),
  "config$is_log_transformed recorded as FALSE")
report(identical(dim(out$state$expression), dim(e)),
  "conversion does not change the shape of the matrix")

e <- exprs
report(!threw(init_state(mk_state(e), cfg0, minimal=TRUE)),
  "input with no zeros passes unchanged")

###############################################################################
section("negative raw values are an error")

e <- exprs
e[3, 3] <- -1.5
report(threw(init_state(mk_state(e), cfg0, minimal=TRUE)),
  "negative raw value is an error")

## f.err() stops with "Stopping" and writes the detail to the log:
txt <- readLines(log_file)
report(any(grepl("negative values in state$expression", txt, fixed=TRUE)),
  "error message names the problem")
report(any(grepl("-1.5", txt, fixed=TRUE)),
  "error message reports the offending value")
report(any(grepl(rownames(exprs)[3], txt, fixed=TRUE)),
  "error message names the offending feature")

###############################################################################
section("already transformed input")

e <- exprs
e[3, 3] <- -1.5
e[1, 1] <- 0

cfg <- cfg0
cfg$is_log_transformed <- TRUE
out <- try(init_state(mk_state(e), cfg, minimal=TRUE), silent=TRUE)
report(!inherits(out, "try-error"),
  "is_log_transformed=TRUE accepts negative values")
if(!inherits(out, "try-error")) {
  report(out$state$expression[1, 1] %in% 0,
    "is_log_transformed=TRUE leaves zeros alone")
  report(identical(out$config$is_log_transformed, TRUE),
    "config$is_log_transformed stays TRUE")
}

## config$normalization_method no longer says anything about the scale of the
##   input: 'none' asks for no normalization, which is a separate question from
##   what scale the values are already on. Pinned because it used to imply the
##   input was already transformed, while the rest of the workflow went on
##   treating it as raw:

cfg <- cfg0
cfg$normalization_method <- "none"
report(threw(init_state(mk_state(e), cfg, minimal=TRUE)),
  "normalization_method 'none' alone does not exempt negative values")
txt <- readLines(log_file)
report(any(grepl("set config$is_log_transformed", txt, fixed=TRUE)),
  "error message points at config$is_log_transformed")

## the parameter this one replaced is gone from the vocabulary, so a config
##   written against an older version is refused by name, as is the other key
##   retired at the same time:

cfg <- cfg0
cfg$zeros_to_na <- FALSE
report(threw(check_config(cfg)), "retired config$zeros_to_na is an error")
txt <- readLines(log_file)
report(any(grepl("unexpected parameter name: zeros_to_na", txt, fixed=TRUE)),
  "error message names the offending parameter")

cfg <- cfg0
cfg$feature_aggregation_scaled <- FALSE
report(threw(check_config(cfg)), "retired config$feature_aggregation_scaled is an error")
txt <- readLines(log_file)
report(any(grepl("unexpected parameter name: feature_aggregation_scaled", txt, fixed=TRUE)),
  "error message names that one too")

###############################################################################
section("initialized=TRUE skips the check")

## an already initialized state may hold negative transformed values, so
##   re-checking it would reject valid data:

e <- exprs
e[3, 3] <- -1.5
cfg <- cfg0
cfg$n_samples_expr_col <- "n_samps_expr"
cfg$median_raw_col <- "median_raw"
cfg$n_features_expr_col <- "n_feats_expr"
report(!threw(init_state(mk_state(e), cfg, initialized=TRUE)),
  "initialized=TRUE skips the negative value check")

###############################################################################
section("normalize() cross-check")

cfg <- cfg0
cfg$is_log_transformed <- TRUE
cfg$normalization_method <- "none"
cfg$feat_col <- "pep"
cfg$obs_col <- "obs"
cfg$save_state <- FALSE
out <- init_state(mk_state(exprs), cfg, minimal=TRUE)

report(threw(normalize(out$state, out$config, method="RLE")),
  "normalize() rejects a method that would transform log data again")
txt <- readLines(log_file)
report(any(grepl("would transform it again", txt, fixed=TRUE)),
  "normalize() error explains the contradiction")
report(!threw(normalize(out$state, out$config, method="none")),
  "normalize() allows method 'none' when is_log_transformed is TRUE")

## the ordinary case: raw input, so any method is fine:
out <- init_state(mk_state(exprs), cfg0, minimal=TRUE)
cfg <- out$config
cfg$feat_col <- "pep"
cfg$obs_col <- "obs"
cfg$save_state <- FALSE
report(!threw(normalize(out$state, cfg, method="RLE")),
  "normalize() allows any method when is_log_transformed is FALSE")

###############################################################################
section("downstream screens count NA, not sign")

## after log2(x + 1) a value can be 0 only if the raw value was 0, and those are
##   already NA; a negative value (e.g. vsn output) is real data and has to be
##   counted as measured:

e <- matrix(c(-1, 0, 5, 2, NA, NA), nrow=2, byrow=TRUE,
  dimnames=list(c("p1", "p2"), c("o1", "o2", "o3")))
st <- list(
  expression=e,
  features=data.frame(pep=rownames(e), gene=rownames(e)),
  samples=data.frame(obs=colnames(e), grp=c("ctl", "trt", "trt"))
)
cfg <- list(feat_col="pep", obs_col="obs", log_file=log_file)

report(identical(unname(samples_per_feature(st, cfg)), c(3L, 1L)),
  "samples_per_feature() counts non-NA values regardless of sign")
report(identical(unname(features_per_sample(st, cfg)), c(2L, 1L, 1L)),
  "features_per_sample() counts non-NA values regardless of sign")

## p1 is (-1, 0, 5) and p2 is (2, NA, NA), so the medians are 0 and 2: the
##   median is over the measured values, and a zero among them is one of them:

m <- feature_median_expression(st, cfg)
report(identical(unname(m), c(0, 2)),
  "feature_median_expression() takes the median over the measured values only")
report(identical(names(m), rownames(st$expression)),
  "and names the result by feature")
report(length(m) == nrow(st$expression) && is.numeric(m),
  "returning one number per feature")

e_na <- rbind(st$expression, p3=c(NA, NA, NA))
st_na <- st
st_na$expression <- e_na
st_na$features <- rbind(st$features, data.frame(pep="p3", gene="p3"))
m_na <- feature_median_expression(st_na, cfg)
report(is.na(m_na[["p3"]]) && identical(unname(m_na[1:2]), c(0, 2)),
  "a feature measured in no sample gets NA rather than zero")

st_df <- st
st_df$expression <- as.data.frame(st$expression)
report(threw(feature_median_expression(st_df, cfg)),
  "feature_median_expression() refuses an expression that is not a matrix")

st2 <- filter_features(st, cfg, n_samples_min=2, remove_constant=FALSE,
  filter_by_formula=FALSE)
report(identical(st2$features$pep, "p1"),
  "filter_features() keeps the feature with a negative and a zero value")

st3 <- filter_observations(st, cfg, n_features_min=2, remove_constant=FALSE)
report(identical(st3$samples$obs, "o1"),
  "filter_observations() counts non-NA values regardless of sign")

###############################################################################
section("combine_replicates() propagates NA for an all-missing replicate group")

## p1 is measured in both replicates of s1, in neither replicate of s2, and in
##   one replicate of s3. fn is applied with na.rm=TRUE, so the s2 group is
##   combined from no values at all, and only stats::median returns NA for that
##   on its own; sum returns 0 and mean returns NaN:

e <- matrix(c(4, 6, NA, NA, 8, NA,
              1, 3,  5,  7, 9, 11), nrow=2, byrow=TRUE,
  dimnames=list(c("p1", "p2"), paste0("o", 1:6)))
st <- list(
  expression=e,
  features=data.frame(pep=rownames(e), gene=rownames(e)),
  samples=data.frame(obs=colnames(e), sid=rep(c("s1", "s2", "s3"), each=2),
    grp=rep(c("ctl", "trt", "trt"), each=2))
)
cfg <- list(
  feat_col="pep", feat_id_col="pep", gene_id_col="gene", obs_id_col="obs",
  sample_id_col="sid", frm=~grp, test_term="grp", save_state=FALSE,
  log_file=log_file
)

for(nom in c("median", "sum", "mean")) {
  out <- combine_replicates(st, cfg, fn=get(nom))
  ex <- out$state$expression
  report(is.na(ex["p1", "s2"]),
    paste0("fn=", nom, ": all-missing replicate group combines to NA"))
  report(!isTRUE(ex["p1", "s2"] %in% 0),
    paste0("fn=", nom, ": all-missing replicate group is not recorded as 0"))
  report(!any(is.nan(ex)),
    paste0("fn=", nom, ": no NaN left in the combined matrix"))
  report(identical(unname(ex["p2", ]), c(get(nom)(c(1, 3)), get(nom)(c(5, 7)),
    get(nom)(c(9, 11)))),
    paste0("fn=", nom, ": fully observed groups combine as fn dictates"))
}

## the partially observed group uses only the value that was measured:

out <- combine_replicates(st, cfg, fn=sum)
report(isTRUE(out$state$expression["p1", "s3"] %in% 8),
  "fn=sum: partially observed group combines from the observed value only")
report(isTRUE(out$state$expression["p1", "s1"] %in% 10),
  "fn=sum: fully observed group is unaffected by the NA mask")

###############################################################################
section("combine_features() keeps a feature whose values are all zero")

## the median polish and robust summary paths drop features with no measured
##   values; a feature measured as 0 in every sample has been measured, and on a
##   log scale 0 is an ordinary value:

## gene g1 has one feature measured as 0 everywhere and one with real values.
##   When zeros counted as missing, p1 was dropped and the summary came from p2
##   alone; p1 was measured, so both features contribute. The data here are log
##   scale, which is what makes a 0 an ordinary value, and which
##   combine_features() now requires, so cfg declares it below:

e <- matrix(c(0, 0, 0, 0,
              5, 6, 7, 8), nrow=2, byrow=TRUE,
  dimnames=list(c("p1", "p2"), paste0("o", 1:4)))
st <- list(
  expression=e,
  features=data.frame(pep=rownames(e), gene=c("g1", "g1")),
  samples=data.frame(obs=colnames(e), grp=c("ctl", "ctl", "trt", "trt"))
)
cfg <- list(
  feat_col="pep", feat_id_col="pep", gene_id_col="gene", obs_col="obs",
  obs_id_col="obs", sample_id_col="obs", frm=~grp, test_term="grp",
  is_log_transformed=TRUE, save_state=FALSE, log_file=log_file
)

out <- try(combine_features(st, cfg, method="robustSummary"), silent=TRUE)
report(!inherits(out, "try-error") &&
  isTRUE(all.equal(unname(out$state$expression["g1", ]), c(2.5, 3, 3.5, 4))),
  "combine_features(robustSummary) summarizes from the all-zero feature too")
report(!inherits(out, "try-error") &&
  !isTRUE(all.equal(unname(out$state$expression["g1", ]), c(5, 6, 7, 8))),
  "combine_features(robustSummary) no longer drops the all-zero feature")

## an all-NA feature is still dropped, so a gene keeps being summarized from the
##   features that were measured:

st$expression <- e
st$expression["p1", ] <- NA
out <- try(combine_features(st, cfg, method="robustSummary"), silent=TRUE)
report(!inherits(out, "try-error") &&
  isTRUE(all.equal(unname(out$state$expression["g1", ]), c(5, 6, 7, 8))),
  "combine_features(robustSummary) still drops an all-NA feature")

###############################################################################
section("combine_features(robustSummary) falls back for an all-zero sample column")

## MsCoreUtils::robustSummary() drops a sample whose every measured value in the
##   group is exactly 0 and reports it as NA, which is indistinguishable from a
##   sample where the gene was never measured. Neither feature here is constant
##   -- four distinct observed values each -- so no upstream filter catches it:

e <- matrix(c(0, 5, 6, 7,
              0, 8, 9, 10), nrow=2, byrow=TRUE,
  dimnames=list(c("p1", "p2"), paste0("o", 1:4)))
st <- list(
  expression=e,
  features=data.frame(pep=rownames(e), gene=c("g1", "g1")),
  samples=data.frame(obs=colnames(e), grp=c("ctl", "ctl", "trt", "trt"))
)

report(identical(unname(apply(e, 1, function(v) length(unique(v[!is.na(v)])))),
  c(4L, 4L)),
  "the offending features are not constant, so no upstream filter removes them")
report(nrow(f.prefilter_features(st)$expression) == 2,
  "the offending features survive f.prefilter_features()")
report(nrow(filter_features(st, c(cfg, list(n_samples_min=2)), n_samples_min=2,
  remove_constant=TRUE, filter_by_formula=FALSE)$expression) == 2,
  "the offending features survive filter_features(remove_constant=TRUE)")

## such a group is summarized with MsCoreUtils::medianPolish() instead, which
##   fits the same additive model by medians and so has no design to shed:

m0 <- length(readLines(log_file))
out <- try(combine_features(st, cfg, method="robustSummary"), silent=TRUE)
report(!inherits(out, "try-error"),
  "combine_features(robustSummary) no longer errors on an all-zero sample column")
report(!inherits(out, "try-error") && !anyNA(out$state$expression["g1", ]),
  "the affected sample gets a summary rather than a silent NA")
report(!inherits(out, "try-error") && isTRUE(all.equal(out$state$expression["g1", ],
  MsCoreUtils::medianPolish(e, na.rm=TRUE, maxiter=30))),
  "the fallback summary is exactly what medianPolish() returns for the group")
report(!inherits(out, "try-error") &&
  identical(out$state$features$combine_method, "medianPolish"),
  "state$features$combine_method records the route taken for the gene")
txt <- readLines(log_file)[-seq_len(m0)]
report(any(grepl("exactly 0", txt, fixed=TRUE)) &&
  any(grepl("medianPolish", txt, fixed=TRUE)),
  "the log explains why the group was routed to medianPolish()")
report(any(grepl(" 1 of 1 gene group", txt, fixed=TRUE)) &&
  any(grepl("1 such (gene, sample) cell", txt, fixed=TRUE)),
  "the log counts the affected gene groups and (gene, sample) cells")
report(any(grepl("g1", txt, fixed=TRUE)), "the log names the affected gene")

## a group in which every value is zero would empty robustSummary()'s design and
##   stop inside MASS::rlm(); the same check routes it to medianPolish() first:

st$expression <- matrix(0, nrow=2, ncol=4,
  dimnames=list(c("p1", "p2"), paste0("o", 1:4)))
out <- try(combine_features(st, cfg, method="robustSummary"), silent=TRUE)
report(!inherits(out, "try-error") && all(out$state$expression["g1", ] %in% 0),
  "an all-zero gene group is summarized as zero rather than refused")

## a sample that was genuinely not measured is not affected; NA is the right
##   answer there, and robustSummary() reaches it honestly:

st$expression <- e
st$expression[, "o1"] <- NA
out <- try(combine_features(st, cfg, method="robustSummary"), silent=TRUE)
report(!inherits(out, "try-error") && is.na(out$state$expression["g1", "o1"]),
  "an unmeasured sample is summarized as NA")
report(!inherits(out, "try-error") &&
  identical(out$state$features$combine_method, "robustSummary"),
  "an unaffected group keeps the robustSummary route")

## a single retained feature bypasses robustSummary()'s design entirely, so an
##   all-zero column there is not affected and does not need the fallback:

st$expression <- matrix(c(0, 5, 6, 7, NA, NA, NA, NA), nrow=2, byrow=TRUE,
  dimnames=list(c("p1", "p2"), paste0("o", 1:4)))
out <- try(combine_features(st, cfg, method="robustSummary"), silent=TRUE)
report(!inherits(out, "try-error") &&
  identical(out$state$features$combine_method, "robustSummary"),
  "a single-feature group with a zero keeps the robustSummary route")

## only the affected groups are rerouted; the rest are summarized as before:

e2 <- rbind(e, p3=c(3, 4, 5, 6), p4=c(4, 5, 6, 8))
st2 <- list(
  expression=e2,
  features=data.frame(pep=rownames(e2), gene=c("g1", "g1", "g2", "g2")),
  samples=data.frame(obs=colnames(e2), grp=c("ctl", "ctl", "trt", "trt"))
)
out <- try(combine_features(st2, cfg, method="robustSummary"), silent=TRUE)
report(!inherits(out, "try-error") &&
  identical(out$state$features$combine_method, c("medianPolish", "robustSummary")),
  "one route is recorded per gene, not one for the whole call")
report(!inherits(out, "try-error") && isTRUE(all.equal(out$state$expression["g2", ],
  MsCoreUtils::robustSummary(e2[c("p3", "p4"), ]))),
  "an unaffected gene in the same call is still summarized by robustSummary()")

## the route is recorded for method="medianPolish" too, and under whatever name
##   config$combine_method_col gives it; an existing column is left alone:

out <- try(combine_features(st2, cfg, method="medianPolish"), silent=TRUE)
report(!inherits(out, "try-error") &&
  identical(out$state$features$combine_method, c("medianPolish", "medianPolish")),
  "method='medianPolish' records its own route for every gene")

cfg2 <- c(cfg, list(combine_method_col="agg_route"))
out <- try(combine_features(st2, cfg2, method="robustSummary"), silent=TRUE)
report(!inherits(out, "try-error") &&
  ("agg_route" %in% names(out$state$features)) &&
  !("combine_method" %in% names(out$state$features)),
  "config$combine_method_col names the column")

st3 <- st2
st3$features$combine_method <- "mine"
m0 <- length(readLines(log_file))
out <- try(combine_features(st3, cfg, method="robustSummary"), silent=TRUE)
report(!inherits(out, "try-error") &&
  identical(out$state$features$combine_method, c("mine", "mine")),
  "an existing column of that name is preserved rather than overwritten")
report(any(grepl("already has a column named combine_method",
  readLines(log_file)[-seq_len(m0)], fixed=TRUE)),
  "preserving an existing column is logged")

## method='none' aggregates nothing, so it has no route to record:

out <- try(combine_features(st2, cfg, method="none"), silent=TRUE)
report(!inherits(out, "try-error") &&
  !("combine_method" %in% names(out$state$features)),
  "method='none' adds no route column")

## medianPolish is unaffected and keeps working on the same input:

st$expression <- matrix(0, nrow=2, ncol=4,
  dimnames=list(c("p1", "p2"), paste0("o", 1:4)))
out <- try(combine_features(st, cfg, method="medianPolish"), silent=TRUE)
report(!inherits(out, "try-error") &&
  all(is.finite(out$state$expression["g1", ])),
  "medianPolish summarizes an all-zero gene group without complaint")

###############################################################################
section("no stray warnings escape either aggregator")

## helper: runs a call, returning its value and every warning raised:
with_warnings <- function(expr) {
  ws <- NULL
  v <- withCallingHandlers(try(expr, silent=TRUE),
    warning=function(e) { ws <<- c(ws, conditionMessage(e)); invokeRestart("muffleWarning") })
  list(v=v, warnings=ws)
}

## MsCoreUtils::robustSummary()'s formals are (x, ...), so an na.rm passed in goes
##   to MASS::rlm(), which has none either and warns on every group summarized.
##   robustSummary() masks NA itself, so the values do not depend on it:

st4 <- list(
  expression=rbind(p1=c(15.1, 16.2, 17.3, 16.8), p2=c(14.2, NA, 16.1, 15.5)),
  features=data.frame(pep=c("p1", "p2"), gene=c("g1", "g1")),
  samples=data.frame(obs=paste0("o", 1:4), grp=c("ctl", "ctl", "trt", "trt"))
)
colnames(st4$expression) <- paste0("o", 1:4)
got <- with_warnings(combine_features(st4, cfg, method="robustSummary"))
report(!any(grepl("do not match", got$warnings, fixed=TRUE)),
  "combine_features(robustSummary) raises no 'some of ... do not match' warning")
report(is.null(got$warnings),
  "combine_features(robustSummary) raises no warning at all on an ordinary group")
report(!inherits(got$v, "try-error") &&
  isTRUE(all.equal(got$v$state$expression["g1", ],
    MsCoreUtils::robustSummary(st4$expression))),
  "dropping the stray na.rm leaves robustSummary()'s values unchanged")

## stats::medpolish() warns when it runs out of sweeps. Its criterion compares the
##   sum of absolute residuals between successive sweeps, so a group whose residual
##   sum decays by a constant fraction never satisfies it however long it runs: here
##   the sum halves exactly each sweep, so the warning fires at any maxiter while
##   the summary itself settles by sweep 30. Reported through the log instead:

e5 <- rbind(p1=c(0, 19.3901, 19.0056, 19.0114, 17.9371, 17.5267),
  p2=c(0, NA, NA, NA, NA, NA))
colnames(e5) <- paste0("o", 1:6)
st5 <- list(
  expression=e5,
  features=data.frame(pep=c("p1", "p2"), gene=c("g1", "g1")),
  samples=data.frame(obs=colnames(e5), grp=c("ctl", "ctl", "ctl", "trt", "trt", "trt"))
)
report(!is.null(with_warnings(MsCoreUtils::medianPolish(e5, na.rm=TRUE,
  maxiter=30))$warnings),
  "the fixture does make medpolish() run out of sweeps, so the case is real")

m0 <- length(readLines(log_file))
got <- with_warnings(combine_features(st5, cfg, method="robustSummary"))
txt <- readLines(log_file)[-seq_len(m0)]
report(!inherits(got$v, "try-error"),
  "combine_features(robustSummary) summarizes such a group without error")
report(!any(grepl("did not converge", got$warnings, fixed=TRUE)),
  "medpolish()'s non-convergence warning does not escape to the caller")
report(any(grepl("did not meet stats::medpolish()'s convergence criterion", txt,
  fixed=TRUE)),
  "it is reported through the log instead")
report(any(grepl("describes the sweeps rather than the", txt, fixed=TRUE)),
  "the log says the criterion, not the fit, is what went unmet")
report(any(grepl("1 of 1 gene group(s) summarized by", txt, fixed=TRUE)),
  "the log counts the affected groups against those actually median polished")
report(any(grepl("30 sweeps", txt, fixed=TRUE)),
  "the log names the sweep limit that was reached")
report(any(grepl("g1", txt, fixed=TRUE)), "the log names the affected gene")

## the same reporting covers method="medianPolish", which has always been able to
##   reach this and warned once per gene before:

m0 <- length(readLines(log_file))
got <- with_warnings(combine_features(st5, cfg, method="medianPolish"))
txt <- readLines(log_file)[-seq_len(m0)]
report(!any(grepl("did not converge", got$warnings, fixed=TRUE)),
  "method='medianPolish' does not let the warning escape either")
report(any(grepl("f.combine_features_median_polish : 1 of 1 gene group", txt,
  fixed=TRUE)),
  "method='medianPolish' reports it under its own name")

## a group that converges says nothing at all:

m0 <- length(readLines(log_file))
got <- with_warnings(combine_features(st4, cfg, method="medianPolish"))
txt <- readLines(log_file)[-seq_len(m0)]
report(is.null(got$warnings) &&
  !any(grepl("convergence criterion", txt, fixed=TRUE)),
  "a converging group is not reported")

## the summaries are returned as they stand, unaltered by the interception:

report(!inherits(got$v, "try-error") &&
  isTRUE(all.equal(got$v$state$expression["g1", ],
    MsCoreUtils::medianPolish(st4$expression, na.rm=TRUE, maxiter=30))),
  "intercepting the warning does not change the summary returned")

## any other warning from the fit is still the caller's to see:

p <- with_warnings(f.median_polish(e5, 30))
report(is.null(p$warnings) && isFALSE(p$v$converged),
  "f.median_polish() muffles the convergence warning and reports it as a flag")
report(isTRUE(f.median_polish(st4$expression, 30)$converged),
  "f.median_polish() reports convergence when medpolish() does converge")

###############################################################################
section("imputation counts only NA as missing")

## a matrix holding zeros but no NAs has nothing to impute; when zeros counted as
##   missing, these two took their zeros as missing values and modelled them:

e <- matrix(c(0, 2, 4, 6, 0, 3, 5, 7), nrow=2, byrow=TRUE,
  dimnames=list(c("p1", "p2"), paste0("o", 1:4)))
st <- list(
  expression=e,
  features=data.frame(pep=rownames(e), gene=rownames(e)),
  samples=data.frame(obs=colnames(e), grp=c("ctl", "ctl", "trt", "trt"))
)
cfg <- list(
  feat_col="pep", feat_id_col="pep", gene_id_col="gene", obs_id_col="obs",
  sample_id_col="obs", frm=~grp, test_term="grp", log_file=log_file
)

out <- try(impute_loess_logit(st, cfg), silent=TRUE)
report(!inherits(out, "try-error") && identical(out$expression, e),
  "impute_loess_logit() returns a zero-containing NA-free matrix unchanged")

out <- try(impute_lls(st, cfg, is_log_transformed=TRUE), silent=TRUE)
report(!inherits(out, "try-error") && identical(out$expression, e),
  "impute_lls() returns a zero-containing NA-free matrix unchanged")

section("filtering everything away is reported as such")

## f.check_state() used to reach its dimnames check on an emptied matrix and
##   complain that rownames were missing, which names the wrong problem: R drops
##   rownames when a matrix is subset to zero rows.

set.seed(101)
e_sparse <- sim1(n_obs=6, n_feats=12, mcar_p=0.75)$mat
st_sparse <- list(
  expression=e_sparse,
  features=data.frame(feature_id=rownames(e_sparse)),
  samples=data.frame(observation_id=colnames(e_sparse))
)
cfg_sparse <- list(feat_col="feature_id", obs_col="observation_id",
  log_file=log_file)

n0 <- length(readLines(log_file))
report(threw(prefilter(st_sparse, cfg_sparse)),
  "prefilter() errors when its filters remove every feature and observation")
txt <- readLines(log_file)
txt <- txt[(n0 + 1):length(txt)]
report(any(grepl("nothing is left to work with", txt, fixed=TRUE)),
  "the error says nothing is left rather than blaming missing rownames")
report(!any(grepl("has no rownames", txt, fixed=TRUE)),
  "the misleading rownames message is not what fires")
report(any(grepl("0 features and 0 observations", txt, fixed=TRUE)),
  "the error reports both remaining dimensions")

## regression guard: the sentinel idiom must not come back. NA is the only
##   indicator of a missing value, so an expression value is never tested
##   against 0 to decide whether it is missing:

src <- unlist(lapply(list.files(r_dir, pattern="[.]R$", full.names=TRUE),
  readLines))
src <- src[!grepl("^\\s*#", src)]
report(!any(grepl("is.na(v) | v %in% 0", src, fixed=TRUE)),
  "no 'is.na(v) | v %in% 0' idiom remains in the package sources")
report(!any(grepl("is.na(y) | y %in% 0", src, fixed=TRUE)),
  "no 'is.na(y) | y %in% 0' idiom remains in the package sources")
report(!any(grepl("sum(v > 0, na.rm", src, fixed=TRUE)),
  "no 'sum(v > 0, na.rm=T)' measured-value count remains in the sources")

###############################################################################

cat("\n## passes:", n_pass, "; failures:", n_fail, "; elapsed:",
  round(as.numeric(difftime(Sys.time(), t0, units="secs"))), "s\n")
cat("## log file:", log_file, "\n")

quit(save="no", status=if(n_fail > 0) 1 else 0)
