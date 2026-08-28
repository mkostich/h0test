## Tests for the aggregate path of h0testr::test_msqrob(), which config$test_method
##   "msqrob_agg" selects: msqrob2::msqrobAggregate() fits one mixed model per gene over
##   the rows of its features, with the feature and the observation as random effects,
##   giving gene level inference from feature level input. Also covers the gene level
##   nNonZero both paths now report, which used to be read from the aggregated rowData()
##   and was dropped there whenever a gene's features were observed unequally often.

usage <- function(msg=NULL) {
  if(!is.null(msg)) cat("ERROR:", msg, "\n\n", file=stderr())
  cat(
    "Test the aggregate path of h0testr::test_msqrob(): that it returns one row per",
    "gene from feature level input, that the random observation effect is present by",
    "default and that dropping it via config$test_random_obs changes the fit, warns,",
    "and understates the standard error when there is an observation effect to absorb,",
    "that config$test_ridge renames the fitted parameters and is translated,",
    "that config$test_term and config$contrast keep the same meaning as on the",
    "aggregate-then-fit path, that a gene with a single observed feature falls back to",
    "a gene level fit rather than being reported as NA, that already aggregated input",
    "is refused, that nNonZero is the gene's count of observations for both paths,",
    "that test() and tune() carry the results at gene level, and that a continuous",
    "config$test_term is recovered as a slope per unit of the covariate.",
    "",
    "Usage: Rscript test_msqrob_agg.R <r_dir>",
    "",
    "Required positional arguments:",
    "  <r_dir>  Path to the h0testr package R/ source directory; all .R files",
    "             there are sourced (the installed package is not used).",
    "",
    "Optional named arguments: none.",
    "",
    "Requires the msqrob2, QFeatures and lme4 packages to be installed.",
    "",
    "Output: one PASS/FAIL line per assertion to stdout, then a count of passes and",
    "  failures and the elapsed time. Messages from expected errors are written to a",
    "  temporary log file, whose path is reported at the end.",
    "",
    "Exit codes: 0 all assertions passed; 1 one or more failed; 2 usage error;",
    "  3 msqrob2 not installed; 4 QFeatures not installed; 5 lme4 not installed.",
    "",
    "Examples:",
    "  Rscript test_msqrob_agg.R C:/path/to/h0test/h0testr/R",
    "  Rscript test_msqrob_agg.R ../../h0test/h0testr/R",
    "  Rscript test_msqrob_agg.R C:/path/to/h0testr/R > agg.out 2>&1",
    sep="\n", file=stderr()
  )
  quit(status=2)
}

args <- commandArgs(trailingOnly=TRUE)
if(length(args) != 1) usage("wrong number of arguments")

r_dir <- args[1]
if(!dir.exists(r_dir)) usage(paste("r_dir not a directory:", r_dir))

if(!requireNamespace("msqrob2", quietly=TRUE)) {
  cat("ERROR: msqrob2 not installed\n", file=stderr())
  quit(status=3)
}

if(!requireNamespace("QFeatures", quietly=TRUE)) {
  cat("ERROR: QFeatures not installed\n", file=stderr())
  quit(status=4)
}

if(!requireNamespace("lme4", quietly=TRUE)) {
  cat("ERROR: lme4 not installed\n", file=stderr())
  quit(status=5)
}

for(f in list.files(r_dir, pattern="[.]R$", full.names=TRUE)) source(f)

## msqrob2 fits through BiocParallel::bplapply(); kept serial so that a suite run does
##   not compete with itself, two concurrent processes in this stack having segfaulted:
BiocParallel::register(BiocParallel::SerialParam())

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

init <- function(state=state0, cfg=cfg0) initialize(state, cfg, minimal=TRUE)
agg <- function(o) test_msqrob(o$state, o$config, aggregate=TRUE)
plain <- function(o) test_msqrob(o$state, o$config)

o <- init()
res <- suppressMessages(agg(o))

###############################################################################
section("one row per gene, keyed by the gene id")

report(is.list(res) && all(c("hits", "fit") %in% names(res)),
  "returns a list with hits and fit")
report(nrow(res$hits) %in% n_gene, "one row per gene")
report("gene" %in% names(res$hits), "hits is keyed by config$gene_id_col")
report(!("pep" %in% names(res$hits)), "hits carries no feature id column")
report(setequal(as.character(res$hits$gene), unique(gene)),
  "every gene of the input appears exactly once")
report(!any(duplicated(res$hits$gene)), "no gene appears twice")
report(all(c("nNonZero", ".n", "fit_type") %in% names(res$hits)),
  "hits carries nNonZero, .n and fit_type")
report(all(c("logFC", "se", "df", "t", "pval", "adjPval") %in% names(res$hits)),
  "a one column test reports the moderated t shape msqrob2 returns")
report(!any(is.na(res$hits$pval)), "every gene got a p-value")
report(isTRUE(all.equal(res$hits$adjPval,
  stats::p.adjust(res$hits$pval, method="BH")[order(order(res$hits$pval))])) ||
  all(res$hits$adjPval >= res$hits$pval),
  "adjPval is at least pval on every row")

###############################################################################
section("nNonZero is the gene's count of observations, for both paths")

## the count of observations in which the gene was seen at all, which is what a gene
##   level row can report. Computed from the input rather than trusted:

want <- tapply(1:nrow(exprs), gene, function(idxs) {
  sum(apply(!is.na(exprs[idxs, , drop=FALSE]), 2, any))
})
got <- res$hits$nNonZero[match(names(want), res$hits$gene)]
report(isTRUE(all.equal(as.numeric(want), as.numeric(got))),
  "nNonZero counts the observations in which the gene was measured at all")
## .n is what QFeatures::aggregateFeatures() counted, the features of the gene that
##   reached the fit; initialize() prefilters some away, so it is 2 or 3 here and sums
##   to the number of rows the test was handed:
report(all(res$hits$.n %in% 2:3) &&
  sum(res$hits$.n) %in% nrow(o$state$expression),
  ".n is the number of features behind the gene, and they account for every row")

## the regression this replaced: a gene whose features were observed unequally often
##   made QFeatures::aggregateFeatures() drop the per feature count, and reading it
##   errored. Both paths must now survive that:

st <- state0
st$expression[1, 1:4] <- NA                ## one peptide of the first gene, only
o2 <- init(st)
r_agg <- try(suppressMessages(agg(o2)), silent=TRUE)
r_pln <- try(suppressMessages(plain(o2)), silent=TRUE)
report(!inherits(r_agg, "try-error") && nrow(r_agg$hits) %in% n_gene,
  "the aggregate path runs when a gene's features are observed unequally often")
report(!inherits(r_pln, "try-error") && nrow(r_pln$hits) %in% n_gene,
  "so does the aggregate-then-fit path, which used to error there")
report(!inherits(r_pln, "try-error") &&
  all(r_pln$hits$nNonZero <= ncol(exprs)) && all(r_pln$hits$nNonZero > 0),
  "and its nNonZero is a count of observations, in range")

###############################################################################
section("the fitted models and the returned object")

mods <- SummarizedExperiment::rowData(res$fit[["genes"]])$msqrobModels
report(!is.null(mods) && length(mods) %in% n_gene,
  "fit carries one msqrob2 model per gene in rowData(fit[['genes']])")
report(all(vapply(mods, msqrob2::getFitMethod, character(1)) %in%
  c("lmer", "lm", "rlm")), "every model was fit, none left a fitError")
report("features" %in% names(res$fit) && "genes" %in% names(res$fit),
  "the un-aggregated assay is kept alongside the aggregated one")
report(nrow(SummarizedExperiment::assay(res$fit[["features"]])) %in% nrow(exprs),
  "the fitted assay is the feature level one")

## the coefficients the fit carries include the random effect estimates, which is what
##   distinguishes this from the gene level fit:

b <- msqrob2::getCoef(mods[[1]])
report(all(c("(Intercept)", "grptrt", "sexM") %in% names(b)),
  "the fixed effects are named by the design columns")
report(any(grepl("^\\(Intercept\\)pep", names(b))),
  "and the per feature random effect estimates sit alongside them")

###############################################################################
section("the random observation effect")

cfg <- cfg0
cfg$test_random_obs <- FALSE
o3 <- init(cfg=cfg)
m0 <- mark()
res_off <- suppressMessages(agg(o3))

report(logged("(1 | pep)", m0) && !logged("(1 | obs)", m0),
  "test_random_obs=FALSE fits the feature effect alone")
report(logged("WARNING: test_msqrob: config$test_random_obs is FALSE", m0),
  "and warns that the structure is anti-conservative")

m0 <- mark()
invisible(suppressMessages(agg(o)))
report(logged("(1 | pep) + (1 | obs)", m0),
  "the default fits both random effects")
report(!logged("config$test_random_obs is FALSE", m0),
  "and does not warn")

j <- merge(res$hits[, c("gene", "pval", "df", "se")],
  res_off$hits[, c("gene", "pval", "df", "se")], by="gene", suffixes=c(".on", ".off"))
report(nrow(j) %in% n_gene && !isTRUE(all.equal(j$pval.on, j$pval.off)),
  "the two structures give different p-values, so the switch does something")

## but nothing stronger can be claimed on the shared fixture, and the df in particular
##   cannot: sim_design() draws every feature independently within an observation, so there
##   is no shared observation effect for (1 | obs) to absorb, its variance is estimated at
##   essentially zero, and it costs about 0.002 of effective residual df. What is left of the
##   difference is msqrob2's moderation prior, one number fitted from the spread of variances
##   across the whole dataset, which can fall either way. The standard error says the same
##   thing: it moves in both directions across the genes, and never by anything like the
##   margin below, the largest move on the shared fixture being 11%. The 0.8 is the yardstick
##   the structured fixture below is held to, so the two assertions read against each other:

report(any(j$se.off > j$se.on) && all(j$se.off > 0.8 * j$se.on),
  "with no observation effect in the data, the switch alone does not shrink the se")

## so what the warning claims has to be tested against data that has the structure. Give
##   (1 | obs) something to hold: a per-observation offset on the log scale, of the kind a
##   loading difference leaves behind, over six features per gene rather than three. Six is
##   what makes the claims below hold for every gene rather than most: omitting (1 | obs)
##   treats the features of an observation as independent replicates of the group mean, which
##   is the pseudo-replication the warning is about, and the standard error it understates
##   falls by a factor approaching 1/sqrt(n_feat), so three features leave too little room.
##   With six, both claims held for 20 of 20 genes at every offset seed tried, the largest se
##   ratio being 0.68 and the df margin 12 or more. The p-value follows for 18 to 20 of the
##   20, not all, because the point estimate moves between the two fits as well, so it is not
##   asserted here; the se is the mechanism and the se is what is checked:

set.seed(101)
sim_obs <- sim_design(
  sim_samples(factors=list(grp=c("ctl", "trt")), n_per_cell=nsamps),
  frm=~grp, test_term="grp", n_genes=20, n_genes_signif=5,
  effects=2, peps_per_gene=6, reps_per_sample=1,
  p_drop=0.05, mnar_c0=-Inf, mnar_c1=0, mcar_p=0
)
e_obs <- log2(sim_obs$state$expression + 1)
set.seed(7)
e_obs <- e_obs + rep(stats::rnorm(ncol(e_obs), mean=0, sd=1), each=nrow(e_obs))
state_obs <- list(
  expression=e_obs,
  features=data.frame(pep=rownames(e_obs), gene=sim_obs$state$features$gene_id,
    stringsAsFactors=FALSE),
  samples=samps[match(colnames(e_obs), samps$obs), , drop=FALSE]
)

cfg <- cfg0
cfg$test_random_obs <- FALSE
j2 <- merge(
  suppressMessages(agg(init(state=state_obs)))$hits[, c("gene", "se", "df")],
  suppressMessages(agg(init(state=state_obs, cfg=cfg)))$hits[, c("gene", "se", "df")],
  by="gene", suffixes=c(".on", ".off")
)
report(nrow(j2) %in% n_gene && all(j2$se.off < 0.8 * j2$se.on),
  "with one in the data, omitting it understates every gene's se by a fifth or more")
report(all(j2$df.on < j2$df.off),
  "and inflates every gene's denominator df, which is the anti-conservatism warned about")

## an absent key behaves as TRUE, since new_config() sets it TRUE:
cfg <- cfg0
cfg$test_random_obs <- NULL
o4 <- init(cfg=cfg)
m0 <- mark()
res_abs <- suppressMessages(agg(o4))
report(logged("(1 | pep) + (1 | obs)", m0) &&
  isTRUE(all.equal(res_abs$hits$pval, res$hits$pval)),
  "an absent config$test_random_obs behaves as TRUE")

###############################################################################
section("config$test_ridge")

cfg <- cfg0
cfg$test_ridge <- TRUE
o5 <- init(cfg=cfg)
m0 <- mark()
res_rdg <- suppressMessages(agg(o5))

report(nrow(res_rdg$hits) %in% n_gene, "the ridge fit reports one row per gene")
report(logged("with the fixed effects penalized", m0),
  "and the log says the fixed effects were penalized")

mods_rdg <- SummarizedExperiment::rowData(res_rdg$fit[["genes"]])$msqrobModels
b_rdg <- msqrob2::getCoef(mods_rdg[[1]])
report(all(c("ridgegrptrt", "ridgesexM") %in% names(b_rdg)),
  "msqrob2 renames the penalized fixed effects ridge<column>")
report(!("grptrt" %in% names(b_rdg)),
  "so the design column name is not among the fitted names")
report(!isTRUE(all.equal(res_rdg$hits$logFC, res$hits$logFC)),
  "the penalty changes the coefficients, as it is meant to")
report(!any(is.na(res_rdg$hits$pval)),
  "the renaming is translated, so no gene comes back NA for want of a name")

report(isTRUE(all.equal(f.msqrob_parms(c("(Intercept)", "grptrt"), TRUE),
  c("(Intercept)", "ridgegrptrt"))),
  "f.msqrob_parms() prefixes every name but the intercept when ridge is TRUE")
report(isTRUE(all.equal(f.msqrob_parms(c("(Intercept)", "grptrt"), FALSE),
  c("(Intercept)", "grptrt"))),
  "and leaves every name alone when it is FALSE")

## the ridge fit refuses a mean model with fewer than two non-intercept columns, which
##   is msqrob2's own restriction and reaches the caller as an error rather than as NAs:
cfg <- cfg0
cfg$test_ridge <- TRUE
cfg$frm <- ~grp
cfg$reference_levels <- c(grp="ctl")
o6 <- init(cfg=cfg)
report(threw(suppressMessages(agg(o6))),
  "a two level factor alone is refused under ridge, as msqrob2 requires")

###############################################################################
section("what config$test_term and config$contrast mean")

## a three level factor is a joint test of two columns, which msqrob2::hypothesisTest()
##   cannot express, so it goes through f.msqrob_wald() and reports an F:

st <- state0
st$samples$grp3 <- rep(c("a", "b", "c"), length.out=ncol(exprs))
cfg <- cfg0
cfg$frm <- ~grp3+sex
cfg$test_term <- "grp3"
cfg$reference_levels <- c(grp3="a", sex="F")
o7 <- init(st, cfg)
res3 <- suppressMessages(agg(o7))

report(all(c("f_statistic", "df1", "df2", "pval", "adjPval") %in% names(res3$hits)),
  "a joint test reports the F shape")
report(all(res3$hits$df1 %in% 2), "with two numerator degrees of freedom")
report(!any(is.na(res3$hits$pval)),
  "and every gene gets one: the mixed covariance is a Matrix, which is coerced")
report(!("logFC" %in% names(res3$hits)),
  "and no single log fold change, there being no single difference")

## a contrast is one degree of freedom however many coefficients it weights, so msqrob2
##   answers it directly; the contrast naming the tested column must reproduce the term:

cfg <- cfg0
cfg$test_term <- ""
cfg$contrast <- "grptrt"
o8 <- init(cfg=cfg)
res_con <- suppressMessages(agg(o8))
j <- merge(res$hits[, c("gene", "pval", "logFC")],
  res_con$hits[, c("gene", "pval", "logFC")], by="gene", suffixes=c(".t", ".c"))
report(nrow(j) %in% n_gene && isTRUE(all.equal(j$pval.t, j$pval.c)),
  "the contrast grptrt reproduces the test of term grp exactly")
report(isTRUE(all.equal(j$logFC.t, j$logFC.c)),
  "and reports the same log fold change")

###############################################################################
section("the gene level fallback")

## a gene left with one observed feature cannot carry a random feature effect: that
##   effect would be a single unknown confounded with the intercept. Such a gene is
##   fitted at the gene level instead of being reported as NA:

st <- state0
i <- which(gene %in% unique(gene)[1])
st$expression[i[-1], ] <- NA               ## first gene keeps one observed peptide
o9 <- init(st)
m0 <- mark()
res_fb <- suppressMessages(agg(o9))

report(nrow(res_fb$hits) %in% n_gene, "the gene is kept, not dropped")
report(logged("fewer than two observed", m0),
  "and the log says why it took another route")
row <- res_fb$hits[res_fb$hits$gene %in% unique(gene)[1], ]
report(nrow(row) %in% 1 && !is.na(row$pval),
  "it gets a p-value rather than NA")
report(nrow(row) %in% 1 && row$fit_type %in% c("lm", "rlm"),
  "and fit_type records that it was fitted at the gene level")
report(sum(res_fb$hits$fit_type %in% "lmer") %in% (n_gene - 1),
  "every other gene still took the mixed route")

## and that fallback is the fit test_method="msqrob" performs for it:
r_pln <- suppressMessages(plain(o9))
p_fb <- row$pval
p_pl <- r_pln$hits$pval[r_pln$hits$gene %in% unique(gene)[1]]
report(length(p_pl) %in% 1 && length(p_fb) %in% 1 &&
  isTRUE(all.equal(p_fb, p_pl, tolerance=1e-6)),
  "the fallback p-value is what the gene level path reports for that gene")

###############################################################################
section("degenerate genes")

## a gene with no observed value at all: msqrob2 reports NA rather than dropping it,
##   and it is not mistaken for a fallback candidate that could be fitted:

st <- state0
i <- which(gene %in% unique(gene)[2])
st$expression[i, ] <- NA
o10 <- init(st)
res_na <- try(suppressMessages(agg(o10)), silent=TRUE)
report(!inherits(res_na, "try-error"), "a gene with no observed value does not stop it")
if(!inherits(res_na, "try-error")) {
  report(nrow(res_na$hits) <= n_gene && nrow(res_na$hits) >= n_gene - 1,
    "and the other genes are reported")
}

## blank gene ids become their own genes rather than pooling, via f.gene_ids():
st <- state0
st$features$gene[1] <- ""
o11 <- init(st)
res_bl <- try(suppressMessages(agg(o11)), silent=TRUE)
report(!inherits(res_bl, "try-error"), "a blank gene id does not stop it")
if(!inherits(res_bl, "try-error")) {
  report(any(grepl("^unknown_", as.character(res_bl$hits$gene))),
    "and becomes its own gene, named unknown_<feature>")
}

###############################################################################
section("aggregated input is refused")

cfg <- cfg0
cfg$gene_id_col <- "pep"
o12 <- init(cfg=cfg)
m0 <- mark()
report(threw(suppressMessages(agg(o12))),
  "one column named as both the feature and the gene id is refused")
report(logged("the aggregate path needs feature level input", m0),
  "with a message saying why")
report(logged("use that instead", m0),
  "and naming test_method 'msqrob' as the thing to use")

## the same state runs on the aggregate-then-fit path, which is the recommendation:
report(!inherits(try(suppressMessages(plain(o12)), silent=TRUE), "try-error"),
  "and that recommendation works on the same state")

## the same refusal by way of test():
cfg$test_method <- "msqrob_agg"
o13 <- init(cfg=cfg)
report(threw(suppressMessages(test(o13$state, o13$config))),
  "test() refuses it too")

###############################################################################
section("registration and the predicates")

report("msqrob_agg" %in% test_methods(), "msqrob_agg is in test_methods()")
report(f.gene_level_method("msqrob_agg"),
  "f.gene_level_method() calls it a gene level method")
report(!exists("f.test_max_cols"),
  "there is no cap helper to consult, no engine being capped any more")
report(!("msqrob_agg" %in% eval(formals(tune)$test_methods)),
  "and it is deliberately not in tune()'s default test_methods vector")

###############################################################################
section("test() integration")

cfg <- cfg0
cfg$test_method <- "msqrob_agg"
o14 <- init(cfg=cfg)
out <- suppressMessages(test(o14$state, o14$config))

report(all(c("original", "standard", "fit") %in% names(out)),
  "test() returns original, standard and fit")
report(nrow(out$standard) %in% n_gene, "one standard row per gene")
report(setequal(as.character(out$standard$feature), unique(gene)),
  "keyed by gene id, not by feature id")
report(all(c("feature", "expr", "logfc", "stat", "lod", "pval", "adj_pval") %in%
  names(out$standard)), "with the standard columns")
report(!any(is.na(out$standard$expr)),
  "expr is filled in, at gene level")
report("gene_tag" %in% names(out$original) && !("pep_tag" %in% names(out$original)),
  "the metadata joined into original is the gene level table")
report(all(out$original$gene_tag %in% paste0("G_", unique(gene))),
  "and it carries the gene level marker, so the join was at the right level")

## the gene level expr is the mean over the features of a gene and the observations
##   where they were seen, which is what f.feature_means() computes for a gene level
##   method:
want <- tapply(1:nrow(exprs), gene, function(idxs) {
  mean(exprs[idxs, , drop=FALSE], na.rm=TRUE)
})
got <- out$standard$expr[match(names(want), out$standard$feature)]
report(isTRUE(all.equal(as.numeric(want), as.numeric(got))),
  "expr is the mean over the gene's features and observations")

###############################################################################
section("f.tune2() records an unusable combination rather than failing")

## a full config, since f.tune2() runs filter() and impute() before the test and those
##   read many more keys than the tests above need; the count minima are lowered from
##   the defaults, which are meant for real data and would filter this fixture away:

cfg <- new_config()
cfg$feat_id_col <- cfg$gene_id_col <- cfg$feat_col <- "gene"
cfg$obs_id_col <- cfg$sample_id_col <- cfg$obs_col <- "obs"
cfg$frm <- ~grp+sex
cfg$test_term <- "grp"
cfg$reference_levels <- c(grp="ctl", sex="F")
cfg$test_method <- "msqrob_agg"
cfg$normalization_method <- "none"
cfg$impute_method <- "none"
cfg$permute_var <- ""
cfg$n_features_min <- 1
cfg$n_samples_min <- 1
cfg$log_file <- log_file
cfg$save_state <- FALSE
cfg$is_log_transformed <- TRUE

st <- list(expression=exprs, features=data.frame(gene=rownames(exprs),
  stringsAsFactors=FALSE), samples=samps)
o15 <- init(st, cfg)
m0 <- mark()
row <- try(suppressMessages(f.tune2(o15$state, o15$config)), silent=TRUE)

report(!inherits(row, "try-error"), "f.tune2() returns rather than stopping")
if(!inherits(row, "try-error")) {
  report(nrow(row) %in% 1 && is.na(row$nhits) && is.na(row$ntests),
    "and records the combination as untested")
  report(row$test %in% "msqrob_agg", "under the method that could not run")
}
report(logged("needs feature level", m0),
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
out_ct <- suppressMessages(initialize(state_ct, cfg_ct, minimal=TRUE))

tru_ct <- sim_ct$truth[, "age"]            ## dropout can cost a gene all of its peptides
slope_ct <- 1 / stats::sd(out_ct$state$samples$age)      ## the expected logFC, per year

res_ct <- suppressMessages(test_msqrob(out_ct$state, out_ct$config, aggregate=TRUE))
fc_ct <- res_ct$hits$logFC
t_ct <- tru_ct[as.character(res_ct$hits$gene_id)]
q_ct <- stats::p.adjust(res_ct$hits$pval, method="BH")

report(nrow(res_ct$hits) %in% nrow(sim_ct$truth) && !any(is.na(fc_ct)),
  "a continuous test_term is tested on the aggregate path too: one row per gene, no NA")
report(all(sign(fc_ct[t_ct != 0]) == sign(t_ct[t_ct != 0])),
  "the sign of every planted slope is recovered")
report(abs(mean(abs(fc_ct[t_ct != 0])) - slope_ct) < 0.05 * slope_ct,
  "and logFC is the slope per year of age, the planted per SD effect over sd(age)")
report(mean(abs(fc_ct[t_ct == 0])) < 0.1 * slope_ct,
  "while a gene with no age effect is left with a slope near zero")
report(sum(q_ct[t_ct != 0] < 0.05) >= 11 && sum(q_ct[t_ct == 0] < 0.05) <= 2,
  "and the planted genes are the ones that reject at q < 0.05")

## test() is not called on this fixture, since one msqrobAggregate() pass over 24 genes is
##   half a minute and the standardized logfc of a single column continuous term is asserted
##   in 1/test_msqrob_joint.R, 1/test_deqms.R and 1/test_prolfqua_lmer.R.

###############################################################################
cat("\n## passes:", n_pass, "; failures:", n_fail, "; elapsed:",
  round(as.numeric(difftime(Sys.time(), t0, units="secs"))), "s\n")
cat("## log file:", log_file, "\n")

quit(save="no", status=if(n_fail > 0) 1 else 0)
