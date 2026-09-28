## Tests for the mixed path of h0testr::test_prolfqua(), which config$test_method
##   "prolfqua_lmer" selects: one model per gene over the rows of its features, with
##   the feature and the observation as random effects, giving gene level inference
##   from feature level input without aggregating.

usage <- function(msg=NULL) {
  if(!is.null(msg)) cat("ERROR:", msg, "\n\n", file=stderr())
  cat(
    "Test the mixed path of h0testr::test_prolfqua(): that it returns one row per",
    "gene from feature level input, that the random observation effect is present by",
    "default and that dropping it via config$test_random_obs raises the",
    "Satterthwaite denominator degrees of freedom and the rejection rate, that",
    "config$test_term and config$contrast keep the same meaning as on the least",
    "squares path, that a gene with a single observed feature falls back to a least",
    "squares fit rather than being dropped, that already aggregated input is refused,",
    "that the moderation keys are reported as unused, that test_h0() and tune() carry the",
    "results at gene level, and that a continuous config$test_term is recovered as a",
    "slope per unit of the covariate, through test_h0() since the engine reports no",
    "effect size of its own.",
    "",
    "Usage: Rscript test_prolfqua_lmer.R <r_dir>",
    "",
    "Required positional arguments:",
    "  <r_dir>  Path to the h0testr package R/ source directory; all .R files",
    "             there are sourced (the installed package is not used).",
    "",
    "Optional named arguments: none.",
    "",
    "Requires the prolfqua, lmerTest and lme4 packages to be installed.",
    "",
    "Output: one PASS/FAIL line per assertion to stdout, then a count of passes and",
    "  failures and the elapsed time. Messages from expected errors are written to a",
    "  temporary log file, whose path is reported at the end.",
    "",
    "Exit codes: 0 all assertions passed; 1 one or more failed; 2 usage error;",
    "  3 prolfqua not installed; 4 lmerTest not installed; 5 lme4 not installed.",
    "",
    "Examples:",
    "  Rscript test_prolfqua_lmer.R C:/path/to/h0test/h0testr/R",
    "  Rscript test_prolfqua_lmer.R ../../h0test/h0testr/R",
    "  Rscript test_prolfqua_lmer.R C:/path/to/h0testr/R > lmer.out 2>&1",
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

if(!requireNamespace("lmerTest", quietly=TRUE)) {
  cat("ERROR: lmerTest not installed\n", file=stderr())
  quit(status=4)
}

if(!requireNamespace("lme4", quietly=TRUE)) {
  cat("ERROR: lme4 not installed\n", file=stderr())
  quit(status=5)
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

## mark() before the call, logged(pat, since=mark()) after it; f.err() writes its whole
##   message to the log and then stops with "Stopping", so the detail is in the log:
mark <- function() length(readLines(log_file))
logged <- function(pat, since=0) {
  txt <- readLines(log_file)
  if(since >= length(txt)) return(FALSE)
  any(grepl(pat, txt[(since + 1):length(txt)], fixed=TRUE))
}

###############################################################################
## shared data: 20 genes of 3 peptides, 12 observations in two groups of six with a
##   crossed second covariate. Marker columns record which level a metadata column
##   belongs to, so that a gene row carrying pep_tag would be a join to the wrong
##   level:

set.seed(101)
nsamps <- 6
sim <- sim_design(
  sim_samples(factors=list(grp=c("ctl", "trt")), n_per_cell=nsamps),
  frm=~grp, test_term="grp", n_genes=20, n_genes_signif=5,
  effects=2, peps_per_gene=3, reps_per_sample=1,
  p_drop=0.05, mnar_c0=-Inf, mnar_c1=0, mcar_p=0
)
exprs <- log2(sim$state$expression + 1)
gene <- sim$state$features$gene_id
feats <- data.frame(pep=rownames(exprs), gene=gene,
  pep_tag=paste0("T_", rownames(exprs)), gene_tag=paste0("G_", gene),
  stringsAsFactors=FALSE)
samps <- data.frame(
  obs=colnames(exprs),
  grp=c(rep("ctl", nsamps), rep("trt", nsamps)),
  sex=rep(c("M", "F"), nsamps),
  stringsAsFactors=FALSE
)

state0 <- list(expression=exprs, features=feats, samples=samps)
n_gene <- length(unique(gene))

cfg0 <- list(
  obs_id_col="obs", sample_id_col="obs", feat_id_col="pep", gene_id_col="gene",
  frm=~grp+sex, test_term="grp", reference_levels=c(grp="ctl", sex="F"),
  log_file=log_file, is_log_transformed=TRUE
)

init <- function(state=state0, cfg=cfg0) init_state(state, cfg, minimal=TRUE)
mixed <- function(o) test_prolfqua(o$state, o$config, mixed=TRUE)

###############################################################################
section("one row per gene, keyed by the gene id")

o <- init()
m0 <- mark()
res <- mixed(o)

report(nrow(res$hits) %in% n_gene,
  "hits has one row per gene, not one per feature")
report(nrow(exprs) > n_gene,
  "the input really was at feature level (more features than genes)")
report(cfg0$gene_id_col %in% names(res$hits),
  "hits is keyed by config$gene_id_col")
report(!(cfg0$feat_id_col %in% names(res$hits)),
  "hits carries no feature id column")
report(setequal(res$hits$gene, unique(gene)),
  "every gene of the input appears in hits")
report(!anyDuplicated(res$hits$gene), "one row per gene, no duplicates")
report(all(res$hits$factor %in% "grp"), "the factor column names config$test_term")
report(is.numeric(res$hits$p.value) && all(res$hits$p.value >= 0 &
  res$hits$p.value <= 1), "p-values are in [0, 1]")
report(logged("(1|`pep`) + (1|sample)", m0),
  "the log records the mixed formula with both random effects")

###############################################################################
section("the columns a Wald F on a mixed fit has no value for")

report("fit_type" %in% names(res$hits), "hits has a fit_type column")
report(all(res$hits$fit_type %in% c("lmer", "lm")),
  "fit_type is one of lmer or lm")
report(all(!res$hits$moderated), "moderated is FALSE on every row")
report(all(!res$hits$trend), "trend is FALSE on every row")
report(all(is.na(res$hits$df.prior)) && all(is.na(res$hits$s2.prior)),
  "the prior columns are NA")
report(all(is.na(res$hits$F.value.unmod)) && all(is.na(res$hits$p.value.unmod)),
  "the unmoderated columns are NA")
report(isTRUE(all.equal(res$hits$s2.denom, res$hits$Mean.Sq / res$hits$F.value)),
  "s2.denom is the effective denominator variance of the reported F")
report(all(res$hits$df.denom > 0) && !any(is.na(res$hits$df.denom)),
  "df.denom is present and positive on every row")
report(isTRUE(all.equal(res$hits$FDR, stats::p.adjust(res$hits$p.value, "BH"))),
  "FDR is BH over the genes tested")
report(all(diff(res$hits$p.value) >= 0), "hits is sorted by p-value")

###############################################################################
section("the returned fits and coefficients")

report(is.matrix(res$coefs) && nrow(res$coefs) %in% n_gene,
  "coefs has one row per tested gene")
report(identical(colnames(res$coefs), colnames(res$design$X)[res$design$cols_test]),
  "coefs is named by the design columns carrying the test")
report(setequal(rownames(res$coefs), res$hits$gene),
  "coefs rows are the genes of hits")
report(!is.null(res$fit) && "Model" %in% class(res$fit),
  "fit is a prolfqua Model built from the mixed strategy")
report(is.null(res$fit_reduced), "fit_reduced is NULL; no reduced model is fitted")
report(!is.null(res$design$X), "the design used is returned")

###############################################################################
section("the random observation effect is what calibrates the test")

cfg <- cfg0
cfg$test_random_obs <- FALSE
o_off <- init(cfg=cfg)
m0 <- mark()
off <- mixed(o_off)

report(logged("config$test_random_obs is FALSE", m0),
  "dropping the random observation effect is warned about")
report(logged("(1|`pep`)", m0) && !logged("(1|sample)", m0),
  "the formula then carries the feature effect only")

j <- merge(res$hits[, c("gene", "df.denom", "p.value", "fit_type")],
  off$hits[, c("gene", "df.denom", "p.value")], by="gene", suffixes=c(".on", ".off"))
j <- j[j$fit_type %in% "lmer", , drop=FALSE]

report(nrow(j) > 0, "some genes were fitted by lmer, so the comparison is meaningful")

## what the observation effect does to the denominator df is not a drift, it is one of two
##   things per gene, and the earlier median-and-most-genes form of these assertions hid that.
##   Either lme4 estimates the observation variance away, at the boundary of zero, and the fit
##   is the same either way; or the term binds, and the denominator falls from about the number
##   of measurements the gene contributes, features by observations, to about the number of
##   observations. That second case is the pseudo-replication the warning is about: without a
##   term for the observation, the features of a gene are asserted to be independent
##   measurements of it. On this fixture the split was 10 genes each way, the binding ones
##   dropping 22 df (11 for a two-feature gene) and the others moving by under 0.005 df of
##   numerical noise. Neither the split nor the exact df is asserted, only the two shapes and
##   the gap between them, both of which are far from the tolerances used:

n_obs <- ncol(o$state$expression)
d <- j$df.denom.off - j$df.denom.on
binds <- d > 1

report(all(binds | abs(d) < 0.01),
  "the observation effect either binds on a gene or costs it nothing, with no middle")
report(any(binds), "and it binds on at least one gene, so the comparison is not vacuous")
report(all(j$df.denom.on[binds] <= n_obs) && all(j$df.denom.off[binds] > n_obs),
  "where it binds, the denominator falls from above the observation count to at or below it")
report(any(abs(j$p.value.on - j$p.value.off) > 1e-8),
  "so the p-values differ between the two structures")

## the key absent must behave as TRUE, since hand built configs do not carry it:
cfg <- cfg0
cfg$test_random_obs <- TRUE
o_on <- init(cfg=cfg)
on2 <- mixed(o_on)
report(isTRUE(all.equal(res$hits$p.value, on2$hits$p.value)),
  "an absent config$test_random_obs behaves as TRUE")

###############################################################################
section("null calibration: the obs effect controls the rejection rate")

## a per observation effect makes the features of a gene correlated within an
##   observation; without a term for it they are treated as independent measurements
##   and the test rejects too often. Small and seeded, so this is a guard on the
##   substantive claim rather than a precise rate:

set.seed(7)
n_g <- 40
n_p <- 4
mk_null <- function() {
  ex <- matrix(as.numeric(NA), nrow=n_g * n_p, ncol=2 * nsamps)
  for(g in seq_len(n_g)) {
    u_pep <- stats::rnorm(n_p, 0, 2)
    u_obs <- stats::rnorm(2 * nsamps, 0, 1)
    for(p in seq_len(n_p)) {
      ex[(g - 1) * n_p + p, ] <- 20 + u_pep[p] + u_obs +
        stats::rnorm(2 * nsamps, 0, 1)
    }
  }
  rownames(ex) <- paste0("g", rep(seq_len(n_g), each=n_p), "_p", rep(seq_len(n_p), n_g))
  colnames(ex) <- samps$obs
  ft <- data.frame(pep=rownames(ex),
    gene=paste0("g", rep(seq_len(n_g), each=n_p)), stringsAsFactors=FALSE)
  return(list(expression=ex, features=ft, samples=samps))
}

st_null <- mk_null()
cfg <- cfg0
cfg$frm <- ~grp
cfg$reference_levels <- c(grp="ctl")
o_n1 <- init(st_null, cfg)
r_on <- mixed(o_n1)
cfg$test_random_obs <- FALSE
o_n2 <- init(st_null, cfg)
r_off <- mixed(o_n2)

rate_on <- mean(r_on$hits$p.value < 0.05)
rate_off <- mean(r_off$hits$p.value < 0.05)
cat("   rejection rate at 0.05: with obs effect", signif(rate_on, 3),
  "; without", signif(rate_off, 3), "\n")

report(rate_off > rate_on, "the feature-only structure rejects more often")
report(rate_on <= 0.2, "the calibrated structure stays near the nominal level")
report(rate_off >= 0.15, "the feature-only structure is visibly anti-conservative")

###############################################################################
section("config$test_term and config$contrast keep their meaning")

report(all(res$hits$Df %in% 1), "a two level factor gives a 1 df test")

cfg <- cfg0
cfg$frm <- ~grp*sex
o2 <- init(cfg=cfg)
r2 <- mixed(o2)
report(all(r2$hits$Df %in% 2),
  "testing grp in ~grp*sex is a joint 2 df test, as on the least squares path")
report(nrow(r2$hits) %in% n_gene, "and still one row per gene")

samps3 <- samps
samps3$grp3 <- rep(c("a", "b", "c"), 4)
cfg <- cfg0
cfg$frm <- ~grp3
cfg$test_term <- "grp3"
cfg$reference_levels <- c(grp3="a")
st3 <- state0
st3$samples <- samps3
o3 <- init(st3, cfg)
r3 <- mixed(o3)
report(all(r3$hits$Df %in% 2), "a three level factor gives a 2 df test")

## a contrast that is the same hypothesis as the term test must give the same test:
cfg <- cfg0
cfg$test_term <- ""
cfg$contrast <- "grptrt"
o4 <- init(cfg=cfg)
r4 <- mixed(o4)
report(all(r4$hits$Df %in% 1), "a contrast is a 1 df test")
j4 <- merge(res$hits[, c("gene", "p.value")], r4$hits[, c("gene", "p.value")],
  by="gene", suffixes=c(".term", ".con"))
report(isTRUE(all.equal(j4$p.value.term, j4$p.value.con)),
  "the contrast grptrt reproduces the term test of grp")

###############################################################################
section("a gene with one observed feature falls back to least squares")

f2 <- feats
f2$gene[1] <- "solo"
f2$gene_tag[1] <- "G_solo"
st_solo <- state0
st_solo$features <- f2
o5 <- init(st_solo)
r5 <- mixed(o5)

report("solo" %in% r5$hits$gene, "the single feature gene is reported, not dropped")
report(r5$hits$fit_type[r5$hits$gene %in% "solo"] %in% "lm",
  "and is marked as fitted by least squares")
report(nrow(r5$hits) %in% (n_gene + 1), "the other genes are unaffected in number")

## the fallback fit is the feature level least squares fit, so its unmoderated F must
##   match what the least squares path reports for that one feature with moderation off:
cfg <- cfg0
cfg$test_moderate <- FALSE
o6 <- init_state(st_solo, cfg, minimal=TRUE)
plain <- test_prolfqua(o6$state, o6$config)
p_plain <- plain$hits$p.value[plain$hits$pep %in% f2$pep[1]]
p_solo <- r5$hits$p.value[r5$hits$gene %in% "solo"]
report(length(p_plain) %in% 1 && isTRUE(all.equal(p_plain, p_solo)),
  "the fallback p-value equals the unmoderated least squares p-value for that feature")

## every gene having exactly one feature is the aggregated case by another name, and
##   init_state() already treats it as one: gene ids with no duplicates mean there is
##   nothing to aggregate, so it sets config$gene_id_col to config$feat_id_col. The
##   mixed path then has no feature level to model and refuses, which is the same
##   refusal as for input that has been through combine_features():
f3 <- data.frame(pep=feats$pep, gene=paste0("solo_", feats$pep),
  stringsAsFactors=FALSE)
st_all1 <- state0
st_all1$features <- f3
o7 <- init(st_all1)
report(identical(o7$config$gene_id_col, o7$config$feat_id_col),
  "init_state() collapses the id columns when every gene has one feature")
m0 <- mark()
report(threw(mixed(o7)), "so the mixed path refuses that state")
report(logged("the mixed path needs feature level input", m0),
  "naming the condition rather than fitting every gene by least squares")

###############################################################################
section("degenerate features and genes")

## a feature observed nowhere cannot contribute a level, so its gene is counted by
##   what is observed; a gene observed nowhere at all is dropped and reported:
ex4 <- exprs
ex4[gene %in% "gene1", ] <- NA
st4 <- state0
st4$expression <- ex4
o8 <- init(st4)
m0 <- mark()
r8 <- try(mixed(o8), silent=TRUE)
report(!inherits(r8, "try-error"), "a gene with no observed value does not stop the run")
report(logged("genes with no observed value", m0),
  "and is reported as dropped, with a count")
if(!inherits(r8, "try-error")) {
  report(!("gene1" %in% r8$hits$gene), "the unobservable gene is absent from hits")
} else report(FALSE, "the unobservable gene is absent from hits")

## a blank gene id gets an id of its own rather than being pooled with every other:
f5 <- feats
f5$gene[1:2] <- ""
st5 <- state0
st5$features <- f5
o9 <- init(st5)
r9 <- try(mixed(o9), silent=TRUE)
report(!inherits(r9, "try-error"), "blank gene ids do not stop the run")
if(!inherits(r9, "try-error")) {
  report(sum(grepl("^unknown_", r9$hits$gene)) %in% 2,
    "each feature with no gene id becomes its own unknown_ gene")
} else report(FALSE, "each feature with no gene id becomes its own unknown_ gene")

###############################################################################
section("already aggregated input is refused")

cfg <- cfg0
cfg$feat_id_col <- cfg$gene_id_col <- "gene"
m0 <- mark()
report(threw(test_prolfqua(o$state, cfg, mixed=TRUE)),
  "one column named as both the feature and the gene id is an error")
report(logged("the mixed path needs feature level input", m0),
  "and the message says what is wrong")
report(logged("test_method 'prolfqua'", m0),
  "and names the method to use instead")

## the same condition by the route a user reaches it: combine_features() sets
##   feat_id_col to gene_id_col when it aggregates
agg <- combine_features(o$state, o$config)
report(agg$config$feat_id_col %in% agg$config$gene_id_col,
  "combine_features() leaves feat_id_col and gene_id_col the same column")
report(threw(test_prolfqua(agg$state, agg$config, mixed=TRUE)),
  "so the mixed path refuses an aggregated state")

## and the least squares path still accepts it, being what it reduces to:
report(!threw(test_prolfqua(agg$state, agg$config)),
  "while the least squares path accepts that state")

###############################################################################
section("the moderation keys are reported as unused")

cfg <- cfg0
cfg$test_moderate <- TRUE
o10 <- init(cfg=cfg)
m0 <- mark()
invisible(mixed(o10))
report(logged("does not", m0) && logged("config$test_moderate", m0),
  "config$test_moderate=TRUE is reported as not consulted")
## as a note rather than a warning: new_config() sets test_moderate=TRUE, so this
##   applies to a default configuration and describes nothing wrong. Matched on the
##   prefix of the line naming the key, since f.quiet_fits() may log a genuine warning
##   about convergence in the same window:
hit <- grep("test_prolfqua: config$test_moderate", readLines(log_file)[-seq_len(m0)],
  fixed=TRUE, value=TRUE)
report(length(hit) %in% 1 && grepl("^NOTE:", hit),
  "and reported as a note, not a warning: new_config() sets that key TRUE")

cfg$test_moderate <- FALSE
cfg$test_trend <- TRUE
o11 <- init(cfg=cfg)
m0 <- mark()
invisible(mixed(o11))
report(logged("config$test_trend", m0),
  "config$test_trend=TRUE is reported as not consulted")

m0 <- mark()
invisible(mixed(o))                      ## cfg0 carries neither key
report(!logged("does not consult", m0),
  "a config carrying neither key gets no note")

###############################################################################
section("registration and the predicates")

report("prolfqua_lmer" %in% test_methods(),
  "prolfqua_lmer is in test_methods()")
report(f.gene_level_method("prolfqua_lmer"),
  "f.gene_level_method() calls it a gene level method")
report(identical(f.test_id_col("prolfqua_lmer", cfg0), cfg0$gene_id_col),
  "f.test_id_col() keys its rows by the gene id")
report(!exists("f.test_max_cols"),
  "it has no cap on the number of coefficients tested at once, and neither has any")

###############################################################################
section("the helpers")

d_term <- f.design_test_cols(o$state, o$config)
L1 <- f.test_L(d_term)
report(nrow(L1) %in% length(d_term$cols_test) &&
  ncol(L1) %in% ncol(d_term$X), "f.test_L() gives one row per tested column")
report(all(L1[cbind(seq_len(nrow(L1)), d_term$cols_test)] %in% 1) &&
  sum(L1) %in% nrow(L1), "and those rows are the identity rows of the tested columns")

d_con <- f.design_test_cols(o4$state, o4$config)
L2 <- f.test_L(d_con)
report(nrow(L2) %in% 1, "a contrast gives a single row")
report(isTRUE(all.equal(unname(L2[1, ]), unname(d_con$contrast))),
  "and that row is the contrast weights")

## f.wald_f() on a least squares fit must agree with stats::anova() of the nested pair
dat <- data.frame(y=exprs[1, ], samps)
X <- stats::model.matrix(~grp+sex, data=samps)
colnames(X) <- make.names(colnames(X), unique=TRUE)
dd <- data.frame(y=exprs[1, ], X)
fit_a <- stats::lm(stats::as.formula(paste("y ~ 0 +",
  paste(colnames(X), collapse=" + "))), data=dd)
Lw <- matrix(0, nrow=1, ncol=ncol(X), dimnames=list(NULL, colnames(X)))
Lw[1, "grptrt"] <- 1
w <- f.wald_f(fit_a, Lw)
fit_b <- stats::lm(stats::as.formula(paste("y ~ 0 +",
  paste(setdiff(colnames(X), "grptrt"), collapse=" + "))), data=dd)
av <- stats::anova(fit_b, fit_a)
report(isTRUE(all.equal(w$F.value, av[["F"]][2])),
  "f.wald_f() on a least squares fit is the nested model F")
report(isTRUE(all.equal(w$p.value, av[["Pr(>F)"]][2])),
  "and gives the same p-value")
report(w$df.denom %in% stats::df.residual(fit_a),
  "with the residual degrees of freedom of the full fit")

Lbad <- Lw
colnames(Lbad)[colnames(Lbad) %in% "grptrt"] <- "not_a_coefficient"
report(is.null(f.wald_f(fit_a, Lbad)),
  "f.wald_f() returns NULL when a tested coefficient is not in the fit")

###############################################################################
section("through test_h0()")

o12 <- init()
tst <- test_h0(o12$state, o12$config, method="prolfqua_lmer")

report(nrow(tst$standard) %in% n_gene, "test_h0() returns one row per gene")
report(setequal(tst$standard$feature, unique(gene)),
  "the standardized table is keyed by gene id")
report(all(c("expr", "logfc", "stat", "pval", "adj_pval") %in% names(tst$standard)),
  "the standardized table has the usual columns")
report(!any(is.na(tst$standard$logfc)),
  "logfc is filled from the coefficients of the fits")
report(!any(is.na(tst$standard$expr)), "expr is filled at gene level")

## expr is the mean over the peptides of the gene as well as over the observations:
mn <- rowMeans(exprs, na.rm=TRUE)
want <- tapply(mn, gene, mean, na.rm=TRUE)
got <- tst$standard$expr[match(names(want), tst$standard$feature)]
report(isTRUE(all.equal(as.numeric(want), as.numeric(got))),
  "and equals the mean over the features of the gene")

report("gene_tag" %in% names(tst$original),
  "gene level metadata is joined to the results")
report(!("pep_tag" %in% names(tst$original)),
  "and feature level metadata, which no gene row has a value for, is not")
report(all(tst$original$gene_tag %in% paste0("G_", tst$original$gene[1:nrow(tst$original)]) |
  grepl("^G_", tst$original$gene_tag)), "the joined metadata belongs to that gene")
report("n_feats" %in% names(tst$original),
  "the features-per-gene count is reported alongside")

## tune() must not aggregate before this method, or there would be no feature level
report(!f.gene_level_method("prolfqua"),
  "the least squares path is not a gene level method, so tune() aggregates for it")

###############################################################################
section("f.tune2() reports rather than stops")

## a full config, since f.tune2() runs filter_state() and impute() before the test and those
##   read many more keys than the tests above need:
cfg <- new_config()
cfg$log_file <- log_file
cfg$save_state <- FALSE
cfg$obs_id_col <- cfg$sample_id_col <- cfg$obs_col <- "obs"
cfg$feat_id_col <- cfg$gene_id_col <- cfg$feat_col <- "gene"
cfg$frm <- ~grp+sex
cfg$test_term <- "grp"
cfg$reference_levels <- c(grp="ctl", sex="F")
cfg$test_method <- "prolfqua_lmer"
cfg$normalization_method <- "none"
cfg$impute_method <- "none"
cfg$permute_var <- ""
cfg$n_features_min <- 1
cfg$n_samples_min <- 1
cfg$is_log_transformed <- TRUE

st6 <- list(expression=exprs, features=data.frame(gene=rownames(exprs),
  stringsAsFactors=FALSE), samples=samps)
o13 <- init_state(st6, cfg, minimal=TRUE)
m0 <- mark()
row <- try(f.tune2(o13$state, o13$config), silent=TRUE)

report(!inherits(row, "try-error"), "f.tune2() does not stop on this combination")
if(!inherits(row, "try-error")) {
  report(is.data.frame(row) && nrow(row) %in% 1, "and returns one row")
  report(is.na(row$nhits) && is.na(row$ntests),
    "with nhits and ntests NA, which tune_check() reads as untested")
} else {
  report(FALSE, "and returns one row")
  report(FALSE, "with nhits and ntests NA, which tune_check() reads as untested")
}
report(logged("needs feature level input", m0),
  "and the log says why the combination was skipped")

###############################################################################
section("a continuous config$test_term")

## every case above tests a factor, whose coefficient is a difference between two of its
##   levels. A continuous covariate's coefficient is a slope, per unit of the covariate as
##   it happens to be recorded, so the planted effect and the reported one are in different
##   units and the conversion between them is what is worth asserting. sim_design() centers
##   and scales its continuous covariates before building the model matrix, so effects=1
##   plants one log2 unit per SD of age, while the engine fits on age as recorded and should
##   report 1/sd(age) per year:

set.seed(101)
samps_ct <- sim_samples(factors=list(sex=c("F", "M")), covariates=list(age=c(20, 60)),
  n_per_cell=9)
sim_ct <- sim_design(samps_ct, frm=~sex + age, test_term="age", n_genes=24,
  n_genes_signif=c(age=12), effects=1, peps_per_gene=6, p_drop=0.3,
  mnar_c0=-Inf, mnar_c1=0, mcar_p=0, log_cv_mean=-2, log_cv_sd=0.2)

state_ct <- sim_ct$state
state_ct$expression <- log2(state_ct$expression + 1)
cfg_ct <- sim_ct$config
cfg_ct$is_log_transformed <- TRUE          ## simulated raw; normalize() is not called here
cfg_ct$log_file <- log_file
out_ct <- suppressMessages(init_state(state_ct, cfg_ct, minimal=TRUE))

tru_ct <- sim_ct$truth[, "age"]            ## dropout can cost a gene all of its peptides
slope_ct <- 1 / stats::sd(out_ct$state$samples$age)      ## the expected logFC, per year

res_ct <- suppressMessages(mixed(out_ct))

report(nrow(res_ct$hits) %in% nrow(sim_ct$truth) && !any(is.na(res_ct$hits$p.value)),
  "a continuous test_term is tested: one row per gene, no missing p-value")
report(all(res_ct$hits$Df %in% 1),
  "a continuous term is one design column, so the F test is on 1 numerator df")
report(!any(c("logFC", "diff") %in% names(res_ct$hits)),
  "the engine reports no effect size of its own here either, only the F test")

## so the slope has to come through test_h0(), which fills logfc from the fitted coefficients.
##   One column carries the test, so that is the signed coefficient rather than the total
##   swing a joint test would report:

cfg_std_ct <- out_ct$config
cfg_std_ct$test_method <- "prolfqua_lmer"
std_ct <- suppressMessages(test_h0(out_ct$state, cfg_std_ct))$standard
fc_ct <- std_ct$logfc
t_ct <- tru_ct[as.character(std_ct$feature)]
q_ct <- stats::p.adjust(std_ct$pval, method="BH")

report(nrow(std_ct) %in% nrow(sim_ct$truth) && !any(is.na(fc_ct)),
  "test_h0() fills a logfc for every gene from those coefficients")
report(all(sign(fc_ct[t_ct != 0]) == sign(t_ct[t_ct != 0])),
  "the sign of every planted slope is recovered")
report(abs(mean(abs(fc_ct[t_ct != 0])) - slope_ct) < 0.05 * slope_ct,
  "and logfc is the slope per year of age, the planted per SD effect over sd(age)")
report(mean(abs(fc_ct[t_ct == 0])) < 0.1 * slope_ct,
  "while a gene with no age effect is left with a slope near zero")
report(sum(q_ct[t_ct != 0] < 0.05) >= 11 && sum(q_ct[t_ct == 0] < 0.05) <= 2,
  "and the planted genes are the ones that reject at q < 0.05")

###############################################################################
cat("\n## passes:", n_pass, "; failures:", n_fail, "; elapsed:",
  round(as.numeric(difftime(Sys.time(), t0, units="secs"))), "s\n")
cat("## log file:", log_file, "\n")

quit(save="no", status=if(n_fail > 0) 1 else 0)
