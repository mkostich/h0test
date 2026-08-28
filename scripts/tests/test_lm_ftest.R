## Tests what h0testr::test_lm() reports for the comparison of its two nested least
##   squares fits, and what it reports for a feature that cannot support that comparison.
##   Three things changed here and each is checked against an independent recomputation:
##   the models are compared with the exact F test of stats::anova() rather than with the
##   likelihood ratio chi-square of lmtest::lrtest(), which treats the residual variance
##   as known and is anti-conservative at the replication these designs have; the F
##   statistic is carried into the standardized table, whose stat column was left empty
##   although the engine had it; and a feature the design cannot support keeps an NA
##   p-value rather than being reported as 1.0, which was a claim about data that was
##   never tested and which inflated the count stats::p.adjust() divides by. Such a
##   feature used to reach stats::lm() and, when it had no observations at all, error out
##   of the apply() taking every other feature's result with it. The set of features left
##   NA is checked to be exactly the set h0testr::filter_features_by_estimability()
##   drops, since one screens what the other tests. Also checks that fdr.method now
##   defaults to "BH", as every other test method here uses, and that "BY" is still
##   available on request.

usage <- function(msg=NULL) {
  if(!is.null(msg)) cat("ERROR:", msg, "\n\n", file=stderr())
  cat(
    "Test the model comparison h0testr::test_lm() reports: that the p-value and the",
    "statistic are those of the exact F test of the two nested fits rather than of the",
    "likelihood ratio chi-square; that the statistic reaches the standardized table",
    "h0testr::test() returns; that a feature the design cannot support is left NA",
    "rather than reported as a p-value of 1.0, for exactly the features the",
    "estimability filter drops, and that such a feature does not count against the",
    "features that were tested; and that fdr.method defaults to \"BH\".",
    "",
    "Usage: Rscript test_lm_ftest.R <r_dir> [--seed=<integer>]",
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
    "  Rscript test_lm_ftest.R C:/path/to/h0testr/R",
    "  Rscript test_lm_ftest.R ../../h0test/h0testr/R",
    "  Rscript test_lm_ftest.R C:/path/to/h0testr/R --seed=7",
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

mark <- function() length(readLines(log_file))
logged <- function(pat, since=0) {
  txt <- readLines(log_file)
  if(since >= length(txt)) return(FALSE)
  any(grepl(pat, txt[(since + 1):length(txt)], fixed=TRUE))
}

close_enough <- function(got, want, tol=1e-9) {
  isTRUE(all.equal(as.numeric(got), as.numeric(want), tolerance=tol))
}

###############################################################################
## the reference implementation the engine is checked against: the same two fits, but
##   built with plain stats::lm() and compared with plain stats::anova(), with no helper
##   of the package involved, so an error in the guard the engine applies cannot hide
##   itself in the reference too. Anything the pair cannot support comes back NA, which
##   is what the engine is asserted to report. Also returns the p-value the likelihood
##   ratio chi-square would have given: for Gaussian least squares fits twice the
##   difference in log likelihood is n * log(RSS_reduced / RSS_full), on the degrees of
##   freedom of the comparison, which is what lmtest::lrtest() reported:

ref_lm <- function(y, X, X_red) {

  out <- c(pval=NA_real_, stat=NA_real_, chisq=NA_real_)

  i <- !is.na(y)
  yy <- y[i]
  xf <- X[i, , drop=FALSE]
  xr <- X_red[i, , drop=FALSE]

  ff <- try(stats::lm(yy ~ xf + 0), silent=TRUE)
  fr <- try(if(ncol(xr) %in% 0) stats::lm(yy ~ 0) else stats::lm(yy ~ xr + 0),
    silent=TRUE)
  if(inherits(ff, "try-error") || inherits(fr, "try-error")) return(out)

  a <- try(stats::anova(fr, ff), silent=TRUE)
  if(inherits(a, "try-error") || !("F" %in% names(a))) return(out)

  fstat <- a[["F"]][2]
  if(length(fstat) < 1 || is.na(fstat) || !is.finite(fstat)) return(out)

  out["stat"] <- fstat
  out["pval"] <- a[["Pr(>F)"]][2]

  rss <- a[["RSS"]]
  out["chisq"] <- stats::pchisq(length(yy) * log(rss[1] / rss[2]),
    df=a[["Df"]][2], lower.tail=FALSE)

  return(out)
}

ref_all <- function(state, config) {
  d <- f.design_test_cols(state, config)
  out <- t(apply(state$expression, 1, ref_lm, d$X, d$X_red))
  rownames(out) <- rownames(state$expression)
  return(out)
}

###############################################################################
## shared data: 24 features over 12 observations in two groups of six, with a crossed
##   second covariate, so the design for ~grp + sex has three columns and rank three.
##   The first eight features carry a group effect, so real p-values as well as null
##   ones are compared. Feature ids say what each is for; the missingness that makes the
##   last four untestable is put in below, after initialize(), so that nothing it does
##   can quietly repair the fixture:

set.seed(seed)

nobs <- 12
nfeat <- 24

samps <- data.frame(
  obs=paste0("o", sprintf("%02d", 1:nobs)),
  grp=rep(c("ctl", "trt"), each=nobs / 2),
  sex=rep(c("F", "M"), nobs / 2),
  stringsAsFactors=FALSE
)

feat_eff <- matrix(0, nrow=nfeat, ncol=nobs)
feat_eff[1:8, samps$grp %in% "trt"] <- 2

exprs <- matrix(stats::rnorm(nfeat * nobs, 20, 1.5), nrow=nfeat) + feat_eff
rownames(exprs) <- c(paste0("f", sprintf("%02d", 1:20)),
  "x_all_na", "x_one_group", "x_saturated", "x_one_resid_df")
colnames(exprs) <- samps$obs
feats <- data.frame(pep=rownames(exprs), gene=rownames(exprs),
  stringsAsFactors=FALSE)

cfg0 <- new_config()
cfg0$obs_id_col <- cfg0$obs_col <- "obs"
cfg0$sample_id_col <- "obs"
cfg0$feat_id_col <- cfg0$feat_col <- "pep"
cfg0$gene_id_col <- "gene"
cfg0$reference_levels <- c(grp="ctl", sex="F")
cfg0$estimability <- "test"
cfg0$df_resid_min <- 1
cfg0$save_state <- FALSE
cfg0$permute_var <- ""
cfg0$is_log_transformed <- TRUE      ## simulated on a log scale; normalize() not called
cfg0$impute_method <- "none"
cfg0$n_features_min <- 5             ## the default of 1000 would filter out everything
cfg0$log_file <- log_file
cfg0$frm <- ~grp + sex
cfg0$test_term <- "grp"

out <- suppressMessages(initialize(list(expression=exprs, features=feats,
  samples=samps), cfg0, minimal=TRUE))
state0 <- out$state
config0 <- out$config

## the four deliberately untestable features. x_all_na is measured nowhere, which is
##   what used to error out of the apply(); x_one_group is measured in one level of grp
##   only, so the term under test has no estimable degrees of freedom left;
##   x_saturated is measured in exactly three observations spanning the rank three
##   design, so there is no residual variance to test against; x_one_resid_df is
##   measured in four such observations, which is the boundary case that must still be
##   tested rather than dropped alongside them:

i_ctl <- which(state0$samples$grp %in% "ctl")
i_trt <- which(state0$samples$grp %in% "trt")

state0$expression["x_all_na", ] <- NA
state0$expression["x_one_group", i_trt] <- NA

keep_sat <- c(i_ctl[state0$samples$sex[i_ctl] %in% "F"][1],
  i_ctl[state0$samples$sex[i_ctl] %in% "M"][1],
  i_trt[state0$samples$sex[i_trt] %in% "F"][1])
state0$expression["x_saturated", -keep_sat] <- NA

keep_one <- c(keep_sat, i_trt[state0$samples$sex[i_trt] %in% "M"][1])
state0$expression["x_one_resid_df", -keep_one] <- NA

design0 <- f.design_test_cols(state0, config0)
report(ncol(design0$X) %in% 3 && f.design_rank(design0$X) %in% 3,
  "fixture: the design for ~grp + sex has three columns and rank three")
report(sum(!is.na(state0$expression["x_saturated", ])) %in% 3 &&
  f.design_rank(design0$X[!is.na(state0$expression["x_saturated", ]), ,
    drop=FALSE]) %in% 3,
  "fixture: x_saturated is measured in exactly as many observations as the rank")

ref0 <- ref_all(state0, config0)
untestable <- c("x_all_na", "x_one_group", "x_saturated")
report(all(is.na(ref0[untestable, "pval"])),
  "fixture: the reference cannot test the three features meant to be untestable")
report(!is.na(ref0["x_one_resid_df", "pval"]),
  "fixture: the reference does test the one residual df feature")

###############################################################################
section("the exact F test replaces the likelihood ratio chi-square")

## in a try(), because the failure this guards against is not a wrong number: a feature
##   measured in no observation used to stop stats::lm() and take the whole call with it:

res <- try(suppressMessages(test_lm(state0, config0)), silent=TRUE)
report(!inherits(res, "try-error") && nrow(res$hits) %in% nfeat,
  "test_lm() returns one row per feature of the state it was given")

if(inherits(res, "try-error")) {
  cat("\n## test_lm() did not return; the remaining assertions cannot be made\n")
  cat("\n## passes:", n_pass, "; failures:", n_fail + 1, "; elapsed:",
    round(as.numeric(difftime(Sys.time(), t0, units="secs"))), "s\n")
  quit(save="no", status=1)
}

hits <- res$hits
rownames(hits) <- hits$feature
hits <- hits[rownames(ref0), , drop=FALSE]

report("stat" %in% names(hits), "the hit table carries a stat column")

report(close_enough(hits$pval[!is.na(ref0[, "pval"])],
  ref0[!is.na(ref0[, "pval"]), "pval"]),
  "the p-value is that of stats::anova() on the same two fits")
report(close_enough(hits$stat[!is.na(ref0[, "stat"])],
  ref0[!is.na(ref0[, "stat"]), "stat"]),
  "and stat is that comparison's F statistic")

## the change is not cosmetic: on this fixture the chi-square is smaller for every
##   feature, which is the anti-conservatism it was replaced for. Checked as a strict
##   inequality rather than as inequality of the vectors, so that a reference that
##   accidentally recomputed the F test could not pass it:

ok <- !is.na(ref0[, "pval"]) & !is.na(ref0[, "chisq"])
report(sum(ok) > 0 && all(ref0[ok, "chisq"] < ref0[ok, "pval"]),
  "the chi-square this replaced is smaller for every feature, so anti-conservative")
report(sum(ok) > 0 && !close_enough(hits$pval[ok], ref0[ok, "chisq"], tol=1e-6),
  "and the reported p-value is not the chi-square one")

## the F statistic against its own p-value, so that stat and pval cannot come from
##   different comparisons: with df_test and df_resid known, one determines the other:

d_test <- f.design_rank(design0$X) - f.design_rank(design0$X_red)
n_ok <- rownames(hits) %in% paste0("f", sprintf("%02d", 1:20))
d_resid <- nobs - f.design_rank(design0$X)
report(close_enough(hits$pval[n_ok],
  stats::pf(hits$stat[n_ok], df1=d_test, df2=d_resid, lower.tail=FALSE)),
  "stat and pval describe the same comparison")

###############################################################################
section("a feature the design cannot support is left NA, not reported as 1.0")

report(all(is.na(hits[untestable, "pval"])),
  "the three untestable features get an NA p-value")
report(all(is.na(hits[untestable, "p.adj"])),
  "and an NA adjusted p-value, rather than the 1.0 they used to be reported as")
report(all(is.na(hits[untestable, "stat"])), "and an NA statistic")
report(!is.na(hits["x_one_resid_df", "pval"]),
  "a feature with one residual df is still tested")

## the whole point of the guard: the feature measured nowhere used to error out of the
##   apply() in test_lm() with '0 (non-NA) cases', losing all 24 results, so a run
##   returning at all is the assertion here:

report(sum(!is.na(hits$pval)) %in% sum(!is.na(ref0[, "pval"])),
  "a feature measured in no observation no longer costs every other feature its result")

## exactly the features the estimability filter drops, since one screens what the other
##   tests: with estimability "test" and df_resid_min 1 the filter keeps a feature when
##   df_test >= 1 and df_resid >= 1, which is the condition the engine now checks:

st_f <- suppressMessages(filter_features_by_estimability(state0, config0,
  estimability="test", df_resid_min=1))
kept <- as.character(st_f$features[[config0$feat_col]])
report(setequal(kept, hits$feature[!is.na(hits$pval)]),
  "the features left NA are exactly those filter_features_by_estimability() drops")

## and the untested features do not count against the tested ones. stats::p.adjust()
##   takes n from the p-values that are not NA, so this is the arithmetic that reporting
##   them as 1.0 got wrong:

ok2 <- !is.na(hits$pval)
report(close_enough(hits$p.adj[ok2], stats::p.adjust(hits$pval[ok2], method="BH")),
  "the adjusted p-values are those of the tested features alone")
report(!close_enough(hits$p.adj[ok2],
  stats::p.adjust(ifelse(is.na(hits$pval), 1.0, hits$pval), method="BH")[ok2]),
  "which is not what filling the untested ones with 1.0 gave")

## the same features, run through test(), reach the standardized table as NA rather
##   than as a hit:

res_t <- suppressMessages(test(state0, config0, method="lm"))
std <- res_t$standard
rownames(std) <- std$feature
report(all(is.na(std[untestable, "pval"])) && all(is.na(std[untestable, "adj_pval"])),
  "test(method='lm') carries the NA p-values into the standardized table")
report(!any(untestable %in% std$feature[!is.na(std$adj_pval) & std$adj_pval < 0.05]),
  "so an untested feature is not counted as a hit")

###############################################################################
section("the statistic reaches the standardized table")

report(!all(is.na(std$stat)),
  "test(method='lm') no longer reports an empty stat column")
report(close_enough(std$stat[!is.na(std$stat)],
  hits[std$feature[!is.na(std$stat)], "stat"]),
  "and the value is the F statistic the engine reported")
report(all(is.na(std$stat) == is.na(std$pval)),
  "stat is present for exactly the features that were tested")

## f.format_lm() now requires the column, rather than filling it with NA, so a table
##   without it is refused instead of silently producing an empty column:

bad <- hits[, setdiff(names(hits), "stat"), drop=FALSE]
report(inherits(try(suppressMessages(f.format_lm(bad, "feature", config0)),
  silent=TRUE), "try-error"),
  "f.format_lm() refuses a table with no stat column")

###############################################################################
section("fdr.method defaults to BH")

report(close_enough(hits$p.adj[ok2], stats::p.adjust(hits$pval[ok2], method="BH")),
  "the default adjustment is BH, as every other test method here uses")
report(!close_enough(hits$p.adj[ok2], stats::p.adjust(hits$pval[ok2], method="BY")),
  "and not the BY it used to default to")

res_by <- suppressMessages(test_lm(state0, config0, fdr.method="BY"))
h_by <- res_by$hits
ok3 <- !is.na(h_by$pval)
report(close_enough(h_by$p.adj[ok3], stats::p.adjust(h_by$pval[ok3], method="BY")),
  "fdr.method='BY' is still available on request")
report(close_enough(h_by$pval[ok3], hits$pval[ok2]),
  "and changes only the adjustment, not the p-values")

m <- mark()
invisible(suppressMessages(test_lm(state0, config0)))
report(logged("fdr.method: BH", since=m),
  "the log says which adjustment was used")

###############################################################################
cat("\n## passes:", n_pass, "; failures:", n_fail, "; elapsed:",
  round(as.numeric(difftime(Sys.time(), t0, units="secs"))), "s\n")

quit(save="no", status=if(n_fail > 0) 1 else 0)
