## Tests that every test_*() function tests the hypothesis config$test_term names,
##   including the marginality rule that testing a variable tests every term
##   containing it. Covers the coefficient selection shared with test_lm() and
##   filter_features_by_estimability(), the refusal by engines that can only test one
##   coefficient at a time, the refusal by test_prolfqua(), which is bounded by terms
##   rather than coefficients, the likelihood ratio path by which test_proda() runs
##   the joint test instead of refusing it, and the skip that keeps tune() going past
##   the engines that cannot.
##   Also checks that the contrast an engine reports is the one
##   config$reference_levels declares, rather than one the engine re-derived by
##   sorting the levels, that tune_check() does not rank a skipped combination
##   above the combinations that actually ran, and that test_proda() tests the
##   intercept under the name proDA gives that column and refuses a reduced model with
##   no parameters in its own words.

usage <- function(msg=NULL) {
  if(!is.null(msg)) cat("ERROR:", msg, "\n\n", file=stderr())
  cat(
    "Test that test_trend(), test_voom(), test_deqms(), test_msqrob(),",
    "test_proda() and test_prolfqua() test the hypothesis named by",
    "config$test_term, agreeing with test_lm() and",
    "filter_features_by_estimability() on which design matrix columns carry that",
    "test. Checks that a joint test over several coefficients is performed by the",
    "limma-family methods, run as a likelihood ratio test against the reduced",
    "model by test_proda(), refused with an informative error by the methods that",
    "can only test one coefficient or one term, and skipped rather than fatal",
    "inside the tune() loop. Also checks that the p-value test_prolfqua() reports",
    "does not depend on where config$test_term sits in config$frm, that the",
    "contrast each engine reports follows config$reference_levels, that",
    "tune_check() ranks a skipped combination below every combination that ran,",
    "and that test_proda() handles the two formulas the other engines accept: a",
    "test of the intercept, and a test_term leaving no parameters in the reduced",
    "model, which it refuses in its own words.",
    "",
    "Requires the prolfqua package to be installed.",
    "",
    "Usage: Rscript test_marginality.R <r_dir>",
    "",
    "Required positional arguments:",
    "  <r_dir>  Path to the h0testr package R/ source directory; all .R files",
    "             there are sourced (the installed package is not used).",
    "",
    "Output: one PASS/FAIL line per assertion to stdout, then a count of",
    "  passes and failures. Log output, including messages from expected",
    "  errors, is written to a temporary file whose path is reported at the end.",
    "",
    "Exit codes: 0 all assertions passed; 1 one or more failed; 2 usage error;",
    "  3 prolfqua not installed.",
    "",
    "Examples:",
    "  Rscript test_marginality.R C:/path/to/h0test/h0testr/R",
    "  Rscript test_marginality.R ../../h0test/h0testr/R",
    "  Rscript test_marginality.R C:/path/to/h0testr/R > test_marginality.out 2>&1",
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

###############################################################################
## Twenty-four observations; grp has three levels, sex and batch two each, age is
##   continuous. Twenty genes with two peptides each, so DEqMS sees varying peptide
##   counts once the first peptide of every other gene is removed.
##   Genes g01-g10 carry a sex effect confined to batch b2, that is an effect that
##   the sexM coefficient alone does not describe. Testing 'sex' in ~sex*batch is a
##   joint test of sexM and sexM:batchb2 by marginality, and only a method that
##   performs that joint test recovers these genes.

set.seed(101)

nobs <- 24
ngene <- 20
npep <- 2

samps <- data.frame(
  obs=paste0("o", sprintf("%02d", 1:nobs)),
  grp=rep(c("a", "b", "c"), each=nobs / 3),
  sex=rep(c("F", "M"), times=nobs / 2),
  batch=rep(c("b1", "b2"), each=nobs / 2),
  age=round(rnorm(nobs, 50, 10), 1),
  stringsAsFactors=FALSE
)

gene_eff <- matrix(0, nrow=ngene, ncol=nobs)
gene_eff[1:10, samps$sex %in% "M" & samps$batch %in% "b2"] <- 4

exprs <- matrix(rnorm(ngene * npep * nobs, 20, 1.5), nrow=ngene * npep)
exprs <- exprs + gene_eff[rep(1:ngene, each=npep), , drop=FALSE]
gene <- paste0("g", sprintf("%02d", rep(1:ngene, each=npep)))
rownames(exprs) <- paste0(gene, ".p", rep(1:npep, times=ngene))
colnames(exprs) <- samps$obs

## make peptides per gene vary, so test_deqms() has a variance prior to fit, and spread
##   rather than merely varied: DEqMS fits that prior as a loess of residual variance on
##   the count, and dropping two rows overall left eighteen of the twenty genes at one
##   count, where the loess comes back NaN and those genes lose their moderated statistic
##   altogether. Dropping the first peptide of every other gene leaves ten genes with one
##   peptide and ten with two, which it fits:

i_drop <- seq(from=1, by=2 * npep, length.out=ngene %/% 2)
exprs <- exprs[-i_drop, , drop=FALSE]
gene <- gene[-i_drop]

feats <- data.frame(pep=rownames(exprs), gene=gene)

## based on new_config() so that the tune() result rows have every column they are
##   built from, and so the filter thresholds are the documented defaults:

cfg0 <- new_config()
cfg0$obs_id_col <- cfg0$obs_col <- "obs"
cfg0$sample_id_col <- "obs"
cfg0$feat_id_col <- cfg0$feat_col <- "pep"
cfg0$gene_id_col <- "gene"
cfg0$frm <- ~grp + sex + batch + age
cfg0$test_term <- "sex"
cfg0$reference_levels <- c(grp="a", sex="F", batch="b1")
cfg0$estimability <- "test"
cfg0$df_resid_min <- 2
cfg0$save_state <- FALSE
cfg0$permute_var <- ""
cfg0$normalization_method <- "RLE"
## the simulated values are already on a log scale, with the gene effects added
##   rather than multiplied, and normalize() is never called below; declared here
##   so that the is_log_transformed=TRUE passed to the engines agrees with it:
cfg0$is_log_transformed <- TRUE
cfg0$impute_method <- "none"
cfg0$n_features_min <- 5    ## default of 1000 would filter out every observation
cfg0$log_file <- log_file

out <- init_state(list(expression=exprs, features=feats, samples=samps), cfg0,
  minimal=TRUE)
state0 <- out$state
config0 <- out$config

true_ids <- paste0("g", sprintf("%02d", 1:10))

## min p-value per gene, from whichever p-value column the engine reports:

recovered <- function(tbl, ids=true_ids, cutoff=0.05) {
  pcol <- intersect(c("P.Value", "pval", "PValue", "p.value"), names(tbl))[1]
  if(is.na(pcol)) return(NA_integer_)
  key <- if("gene" %in% names(tbl)) tbl$gene else
    if("feature" %in% names(tbl)) tbl$feature else
    if("name" %in% names(tbl)) tbl$name else rownames(tbl)   ## proDA::test_diff()
  key <- sub("[.]p[0-9]+$", "", key)     ## peptide ids back to gene ids
  p <- tapply(tbl[[pcol]], key, min, na.rm=TRUE)
  return(sum(p[intersect(ids, names(p))] < cutoff, na.rm=TRUE))
}

## deqms and msqrob take peptide-level input; the rest take gene-level:

## the methods that report gene level results from feature level input are handed the
##   un-aggregated state, which is what tune() does; the predicate is the package's own,
##   so this cannot drift from it. "prolfqua_lmer" is test_prolfqua(mixed=TRUE) and
##   "msqrob_agg" is test_msqrob(aggregate=TRUE):
run_test <- function(nm, cfg) {
  fn <- get(paste0("test_", sub("_(lmer|agg)$", "", nm)))
  if(f.gene_level_method(nm)) {
    st <- state0
    cf <- cfg
  } else {
    agg <- suppressMessages(combine_features(state0, cfg))
    st <- agg$state
    cf <- agg$config
  }
  args <- list(st, cf)
  if(nm %in% "proda") args <- list(st, cf, is_log_transformed=TRUE, prior_df=5)
  if(nm %in% "prolfqua") args <- list(st, cf, is_log_transformed=TRUE)
  if(nm %in% "prolfqua_lmer") {
    args <- list(st, cf, is_log_transformed=TRUE, mixed=TRUE)
  }
  if(nm %in% "msqrob_agg") args <- list(st, cf, aggregate=TRUE)
  return(try(suppressMessages(do.call(fn, args)), silent=TRUE))
}

###############################################################################
section("joint test: frm ~sex*batch, test_term 'sex'")

cfg1 <- config0
cfg1$frm <- ~sex * batch
cfg1$test_term <- "sex"

agg1 <- suppressMessages(combine_features(state0, cfg1))
d1 <- f.design_test_cols(agg1$state, agg1$config)

report(d1$df_intend %in% 2 && length(d1$cols_test) %in% 2,
  "f.design_test_cols(): testing 'sex' in ~sex*batch is a 2 df joint test")
report(identical(sort(colnames(d1$X)[d1$cols_test]), c("sexM", "sexM:batchb2")),
  "f.design_test_cols(): the joint test covers sexM and sexM:batchb2")

## the interaction genes are recoverable only by a method that runs the joint test;
##   before the design-column selection replaced coefficient name matching, trend and
##   voom tested sexM alone and recovered none of them:

for(nm in c("lm", "trend")) {
  res <- run_test(nm, cfg1)
  report(!inherits(res, "try-error") && recovered(res$hits) >= 9,
    paste0("test_", nm, "() recovers the interaction genes"))
}

## voom models a count mean-variance relationship, so it is a poor fit for this
##   log-normal fixture; held to a lower bar than trend, but well above the zero it
##   recovered when it tested sexM alone:

res <- run_test("voom", cfg1)
report(!inherits(res, "try-error") && recovered(res$hits) >= 5,
  "test_voom() recovers most of the interaction genes")

for(nm in c("trend", "voom")) {
  res <- run_test(nm, cfg1)
  report(!inherits(res, "try-error") && "F" %in% names(res$hits),
    paste0("test_", nm, "() reports F, not logFC, for a multi-coefficient test"))
}

## no engine reports one coefficient alone in place of the joint test any more, which
##   would be a different and narrower hypothesis than the one requested. Both engines
##   that used to were capped at their APIs rather than by their fits: msqrob2's cap was
##   in hypothesisTest(), and DEqMS's is in which statistic
##   DEqMS::spectraCounteBayes() forms from the variance prior it fits, not in the prior:

n0 <- length(readLines(log_file, warn=FALSE))

res <- run_test("deqms", cfg1)
report(!inherits(res, "try-error") && "sca.F" %in% names(res$hits) &&
  all(res$hits$sca.df.num %in% 2),
  "test_deqms() puts both sexM and sexM:batchb2 under test, on 2 df")

res <- run_test("msqrob", cfg1)
report(!inherits(res, "try-error") && "f_statistic" %in% names(res$hits),
  "test_msqrob() reports a joint F rather than refusing")
report(!inherits(res, "try-error") && all(res$hits$df1 %in% 2),
  "test_msqrob() puts both sexM and sexM:batchb2 under test")

res <- run_test("msqrob_agg", cfg1)
report(!inherits(res, "try-error") && "f_statistic" %in% names(res$hits),
  "test_msqrob(aggregate=TRUE) reports the same joint F from its mixed fits")
report(!inherits(res, "try-error") && all(res$hits$df1 %in% 2) &&
  !any(is.na(res$hits$pval)),
  "and answers every gene, the mixed covariance being coerced before subsetting")

## there used to be a refusal here, from deqms, and three assertions on the wording of
##   its message: the number of coefficients, their names, and the list of methods it
##   recommended instead. Nothing refuses this test now, so what is checked is that
##   silence: no engine reports a cap, and every one of them reports 2 numerator degrees
##   of freedom rather than quietly testing sexM alone:

txt <- readLines(log_file, warn=FALSE)[-seq_len(n0)]
report(!any(grepl("can test at most", txt, fixed=TRUE)),
  "no engine reports a cap on the number of coefficients")

## every method runs this test. A recommendation that is refused in turn would be worse
##   than none, and there is nothing left to recommend anything else to:

for(nm in c("lm", "trend", "voom", "deqms", "msqrob", "msqrob_agg", "prolfqua",
    "prolfqua_lmer", "proda")) {
  res <- run_test(nm, cfg1)
  report(!inherits(res, "try-error"),
    paste0("test_", nm, "() runs the joint test"))
}

## proDA::test_diff() takes a reduced model as well as a
##   contrast, so test_proda() compares the full model against the model with both
##   tested columns dropped, which is the joint test config$test_term asks for. It is
##   reported as an F statistic with no fold change, there being no single difference
##   to report over two coefficients:

res <- run_test("proda", cfg1)
report(!inherits(res, "try-error"),
  "test_proda() runs the joint test by likelihood ratio instead of refusing it")
report(!inherits(res, "try-error") && "f_statistic" %in% names(res$hits) &&
  !("diff" %in% names(res$hits)),
  "test_proda() reports f_statistic, not diff, for a multi-coefficient test")
report(!inherits(res, "try-error") && all(res$hits$df1 %in% 2),
  "test_proda() reports the joint test as 2 numerator df")
report(!inherits(res, "try-error") && recovered(res$hits) >= 9,
  "test_proda() recovers the interaction genes")

## the reduced model is the full design with the tested columns dropped, so it must
##   be exactly what f.design_test_cols() derived, not a second derivation from
##   formula text that could order factor levels or name interaction columns
##   differently:

report(!inherits(res, "try-error") &&
  identical(colnames(proDA::design(res$fit)),
    sub("^\\(Intercept\\)$", "Intercept", colnames(d1$X))),
  "proDA::proDA() builds the same design f.design_test_cols() did")

## the standardized table carries the F statistic through with an absent logfc,
##   the same convention f.format_limma() uses for limma's F test:

if(!inherits(res, "try-error")) {
  std <- f.format_proda(res$hits, agg1$config)
  report(all(is.na(std$logfc)) && !any(is.na(std$stat)) &&
    identical(names(std), c("feature", "expr", "logfc", "stat", "lod", "pval",
      "adj_pval")),
    "f.format_proda() reports the F statistic with logfc NA")
} else report(FALSE, "f.format_proda() reports the F statistic with logfc NA")

###############################################################################
section("continuous covariate in an interaction: ~age*sex, test_term 'age'")

cfg2 <- config0
cfg2$frm <- ~age * sex
cfg2$test_term <- "age"

agg2 <- suppressMessages(combine_features(state0, cfg2))
d2 <- f.design_test_cols(agg2$state, agg2$config)

report(d2$df_intend %in% 2 &&
  identical(sort(colnames(d2$X)[d2$cols_test]), c("age", "age:sexM")),
  "testing a continuous covariate in an interaction is a 2 df joint test")

res <- run_test("trend", cfg2)
report(!inherits(res, "try-error") && "F" %in% names(res$hits),
  "test_trend() F-tests a continuous covariate and its interaction")

res <- run_test("deqms", cfg2)
report(!inherits(res, "try-error") && all(res$hits$sca.df.num %in% 2),
  "test_deqms() F-tests a continuous covariate and its interaction")

###############################################################################
section("single coefficient: frm ~sex + batch, test_term 'sex'")

cfg3 <- config0
cfg3$frm <- ~sex + batch
cfg3$test_term <- "sex"

agg3 <- suppressMessages(combine_features(state0, cfg3))
d3 <- f.design_test_cols(agg3$state, agg3$config)

report(d3$df_intend %in% 1 && identical(colnames(d3$X)[d3$cols_test], "sexM"),
  "testing 'sex' in ~sex + batch is a 1 df test of sexM")

## every method runs, and for a single coefficient every one of them reports the
##   statistic its own engine reports:

for(nm in c("lm", "trend", "voom", "deqms", "msqrob", "msqrob_agg")) {
  res <- run_test(nm, cfg3)
  report(!inherits(res, "try-error"),
    paste0("test_", nm, "() runs a single-coefficient test"))
}

res <- run_test("trend", cfg3)
report(!inherits(res, "try-error") && "logFC" %in% names(res$hits),
  "test_trend() reports logFC for a single-coefficient test")

###############################################################################
section("multi-level factor: frm ~grp, test_term 'grp'")

cfg4 <- config0
cfg4$frm <- ~grp
cfg4$test_term <- "grp"

res <- run_test("trend", cfg4)
report(!inherits(res, "try-error") && "F" %in% names(res$hits),
  "test_trend() F-tests a three-level factor")

res <- run_test("deqms", cfg4)
report(!inherits(res, "try-error") && all(res$hits$sca.df.num %in% 2),
  "test_deqms() F-tests a three-level factor")

res <- run_test("msqrob", cfg4)
report(!inherits(res, "try-error") && "f_statistic" %in% names(res$hits) &&
  all(res$hits$df1 %in% 2), "test_msqrob() F-tests a three-level factor")

res <- run_test("msqrob_agg", cfg4)
report(!inherits(res, "try-error") && "f_statistic" %in% names(res$hits) &&
  all(res$hits$df1 %in% 2),
  "test_msqrob(aggregate=TRUE) F-tests a three-level factor too")

## a three-level factor is two contrasts, so proda takes the likelihood ratio path
##   here too; the reduced model is the intercept alone:

res <- run_test("proda", cfg4)
report(!inherits(res, "try-error") && "f_statistic" %in% names(res$hits),
  "test_proda() F-tests a three-level factor")
report(!inherits(res, "try-error") && all(res$hits$df1 %in% 2),
  "test_proda() reports a three-level factor as 2 numerator df")

## while a single coefficient still takes the contrast path, reporting a difference
##   whose direction follows config$reference_levels:

res <- run_test("proda", cfg3)
report(!inherits(res, "try-error") && "diff" %in% names(res$hits) &&
  !("f_statistic" %in% names(res$hits)),
  "test_proda() still reports diff for a single-coefficient test")

###############################################################################
section("test_proda(): the intercept, and a reduced model with no parameters")

## the two formulas that every other engine ran but test_proda() did not. Testing the
##   intercept selects the design column stats::model.matrix() names '(Intercept)',
##   which proDA::proDA() renames 'Intercept'; looking the contrast up under the
##   original name found nothing and refused a test that is perfectly well defined:

cfg5 <- config0
cfg5$frm <- ~sex
cfg5$test_term <- "1"

agg5 <- suppressMessages(combine_features(state0, cfg5))
d5 <- f.design_test_cols(agg5$state, agg5$config)

report(d5$df_intend %in% 1 &&
  identical(colnames(d5$X)[d5$cols_test], "(Intercept)"),
  "testing '1' in ~sex is a 1 df test of the intercept column")

res <- run_test("proda", cfg5)
report(!inherits(res, "try-error"), "test_proda() tests the intercept")
report(!inherits(res, "try-error") && "diff" %in% names(res$hits) &&
  !any(is.na(res$hits$pval)),
  "test_proda() runs the intercept test by contrast, with no missing p-values")

## and it is the intercept that was tested, not the sexM coefficient the design also
##   carries: the fixture is log intensities around 20, so the intercept is far from
##   zero for every feature, while sexM is near zero for the genes without an effect:

report(!inherits(res, "try-error") && all(res$hits$diff > 10) &&
  all(res$hits$adj_pval < 0.05),
  "the reported contrast is the intercept, which every feature rejects")

## a test_term that names every term of an intercept-free formula leaves nothing in
##   the reduced model. proDA cannot fit a model with no parameters, and reported the
##   empty matrix from inside proDA::test_diff() as "'d' must be a nonempty numeric
##   vector", which names neither config$frm nor the reduced model:

cfg6 <- config0
cfg6$frm <- ~0 + sex
cfg6$test_term <- "sex"

agg6 <- suppressMessages(combine_features(state0, cfg6))
d6 <- suppressMessages(f.design_test_cols(agg6$state, agg6$config))

report(d6$df_intend %in% 2 &&
  identical(colnames(d6$X)[d6$cols_test], c("sexF", "sexM")) &&
  f.design_rank(d6$X[, -d6$cols_test, drop=FALSE]) %in% 0,
  "dropping 'sex' from ~0 + sex leaves a reduced model with no columns")

n0 <- length(readLines(log_file, warn=FALSE))

res <- run_test("proda", cfg6)
report(inherits(res, "try-error"),
  "test_proda() refuses a reduced model with no parameters")

txt <- readLines(log_file, warn=FALSE)[-seq_len(n0)]
report(any(grepl("leaves a reduced model with no parameters, which proDA cannot fit",
  txt, fixed=TRUE)),
  "the refusal names the empty reduced model rather than a colinear covariate")
report(any(grepl("test_method 'lm', 'trend', 'voom', 'deqms', 'msqrob', 'msqrob_agg',",
  txt, fixed=TRUE)) &&
  any(grepl("'prolfqua' or 'prolfqua_lmer' to test against zero", txt, fixed=TRUE)),
  "the refusal names the methods that can run the test instead")

## which they do: the test is of whether the sex means are zero rather than whether
##   they differ, which f.design_test_cols() warns about but is well defined. deqms is
##   among them now, this being a 2 df joint test of sexF and sexM:

for(nm in c("lm", "trend", "voom", "deqms", "msqrob", "msqrob_agg", "prolfqua",
    "prolfqua_lmer")) {
  res <- run_test(nm, cfg6)
  report(!inherits(res, "try-error"),
    paste0("test_", nm, "() runs the test against zero that proDA cannot fit"))
}

###############################################################################
section("test_prolfqua(): no longer bounded by terms")

## test_prolfqua() used to read the rows of a per-term anova table, which cannot
##   express a joint test over several terms, so a variable that also appears in an
##   interaction was refused. It now compares an explicit full and reduced design,
##   so the joint test runs and agrees with stats::anova() on the same pair:

res <- run_test("prolfqua", cfg1)          ## ~sex*batch, test_term 'sex'
report(!inherits(res, "try-error"),
  "test_prolfqua() runs a test_term that spans two terms")
report(!inherits(res, "try-error") && nrow(res$hits) %in% ngene,
  "the joint test gives one row per gene")
report(!inherits(res, "try-error") && all(res$hits$Df %in% 2),
  "sex in ~sex*batch is a 2 df joint test of sex and sex:batch")

if(!inherits(res, "try-error")) {

  ## the columns dropped are sex and its interaction, so the reduced model is the
  ##   intercept and batch; ~sex+batch would be a 1 df test of the interaction
  ##   alone, which is a different hypothesis:

  agg <- suppressMessages(combine_features(state0, cfg1))
  y <- agg$state$expression["g01", ]
  sm <- agg$state$samples
  ## against p.value.unmod, since the reported p.value has the error variance
  ##   moderated across features and stats::anova() uses the feature's own:

  a <- stats::anova(stats::lm(y ~ batch, data=sm),
    stats::lm(y ~ sex * batch, data=sm))
  report(isTRUE(all.equal(res$hits$p.value.unmod[res$hits$gene %in% "g01"],
    a[2, "Pr(>F)"])),
    "the joint test matches stats::anova() on the same nested pair")
}

## a factor with more than two levels is a single term, and still runs:

res <- run_test("prolfqua", cfg4)          ## ~grp, test_term 'grp'
report(!inherits(res, "try-error") && nrow(res$hits) %in% ngene,
  "test_prolfqua() runs a three-level factor, one row per gene")
report(!inherits(res, "try-error") && all(res$hits$Df %in% 2),
  "test_prolfqua() reports the three-level factor as a single 2 df F-test")

## config$test_term naming the interaction itself is one term, and has to match
##   whichever variable order it was written in; the raw string comparison this
##   replaced returned an empty table without complaint for 'batch:sex':

cfg5 <- config0
cfg5$frm <- ~sex * batch
cfg5$test_term <- "batch:sex"

res <- run_test("prolfqua", cfg5)
report(!inherits(res, "try-error") && nrow(res$hits) %in% ngene,
  "test_prolfqua() matches an interaction test_term written in either order")
## the reported term label is config$test_term verbatim, so it reads the way it was
##   written rather than canonicalized; the raw string comparison this replaced
##   returned an empty table without complaint for 'batch:sex':

report(!inherits(res, "try-error") && all(res$hits$factor %in% "batch:sex"),
  "the rows are labelled with config$test_term as written")

###############################################################################
section("test_prolfqua(): reported p-value does not depend on term order")

## grp and batch are unbalanced against each other in this fixture, so a
##   sequential (Type I) decomposition gives a different p-value for grp depending
##   on whether batch was entered first. Confirmed here with stats::anova() so
##   that the invariance checks below cannot pass vacuously:

cfg6 <- config0
cfg6$frm <- ~grp + batch
cfg6$test_term <- "grp"

cfg7 <- cfg6
cfg7$frm <- ~batch + grp

agg6 <- suppressMessages(combine_features(state0, cfg6))
y <- agg6$state$expression["g01", ]
sm <- agg6$state$samples
p_first <- stats::anova(stats::lm(y ~ grp + batch, data=sm))["grp", "Pr(>F)"]
p_last <- stats::anova(stats::lm(y ~ batch + grp, data=sm))["grp", "Pr(>F)"]

report(!isTRUE(all.equal(p_first, p_last)),
  "fixture is unbalanced: sequential SS for 'grp' depends on term order")

r6 <- run_test("prolfqua", cfg6)
r7 <- run_test("prolfqua", cfg7)
ok <- !inherits(r6, "try-error") && !inherits(r7, "try-error")
report(ok, "test_prolfqua() runs both term orders")

if(ok) {
  k6 <- r6$hits[order(r6$hits$gene), , drop=FALSE]
  k7 <- r7$hits[order(r7$hits$gene), , drop=FALSE]
  report(isTRUE(all.equal(k6$p.value, k7$p.value)),
    "test_prolfqua() reports the same p-value for either term order")

  ## and it is the test adjusted for the other terms, which is what
  ##   config$test_term means; matches the fit with 'grp' entered last, comparing the
  ##   unmoderated column since stats::anova() does not borrow variance:

  report(isTRUE(all.equal(k6$p.value.unmod[k6$gene %in% "g01"], p_last,
    tolerance=1e-6)),
    "the reported p-value is the one adjusted for the other terms")
}

###############################################################################
section("f.tune2() skips instead of aborting")

n0 <- length(readLines(log_file, warn=FALSE))

## a combination a method cannot run yields an NA row, so that one unusable
##   test_method does not end a sweep over the others. The skips for an engine that
##   could not express the hypothesis config$test_term implied are gone, no engine being
##   capped that way any more, so the one provoked here is the live deqms skip, which is
##   a property of the data rather than of config$frm: with one feature per gene the
##   count DEqMS fits its variance prior against does not vary and there is nothing to
##   fit. One peptide of each gene rather than a subset of the genes, so that the
##   twenty-gene guard above it does not fire first and take the credit:

first_pep <- !duplicated(state0$features[[config0$gene_id_col]])
state_u <- state0
state_u$expression <- state0$expression[first_pep, , drop=FALSE]
state_u$features <- state0$features[first_pep, , drop=FALSE]

cfg <- cfg1
cfg$test_method <- "deqms"
row <- try(f.tune2(state_u, cfg), silent=TRUE)
ok <- !inherits(row, "try-error") && is.data.frame(row) && nrow(row) %in% 1 &&
  is.na(row$nhits) && is.na(row$ntests) && row$test %in% "deqms"
report(ok, "f.tune2() returns an NA row for test_method 'deqms'")

txt <- readLines(log_file, warn=FALSE)[-seq_len(n0)]
report(any(grepl("needs the number of features per gene to vary", txt, fixed=TRUE)),
  "f.tune2() logs why the combination was skipped")

## and neither msqrob nor deqms is skipped on the joint test itself any more, so the
##   sweep gets a real row for both:

n0 <- length(readLines(log_file, warn=FALSE))

for(nm in c("msqrob", "deqms")) {
  cfg <- cfg1
  cfg$test_method <- nm
  row <- try(f.tune2(state0, cfg), silent=TRUE)
  report(!inherits(row, "try-error") && is.data.frame(row) && nrow(row) %in% 1 &&
    !is.na(row$ntests) && row$test %in% nm,
    paste("f.tune2() runs test_method", nm, "on a joint test rather than skipping it"))
}

txt <- readLines(log_file, warn=FALSE)[-seq_len(n0)]
report(!any(grepl("can test at most", txt, fixed=TRUE)),
  "f.tune2() has no cap left to report")

## prolfqua is no longer skipped on that joint test: it was the only term-bounded
##   engine, and the term count guard has been removed from f.tune2() along with
##   f.test_max_terms(), which returned Inf for every method, so f.tune2() runs it
##   and returns a real row rather than an NA one:

cfg <- cfg1
cfg$test_method <- "prolfqua"
row <- try(f.tune2(state0, cfg), silent=TRUE)
report(!inherits(row, "try-error") && is.data.frame(row) && nrow(row) %in% 1 &&
  !is.na(row$ntests) && row$test %in% "prolfqua",
  "f.tune2() runs test_method 'prolfqua' on a joint test over two terms")
report(!exists("f.test_max_terms"),
  "f.test_max_terms() is gone rather than returning Inf for every method")
src <- unlist(lapply(list.files(r_dir, pattern="[.]R$", full.names=TRUE), readLines,
  warn=FALSE))
report(!any(grepl("f.test_max_terms", src, fixed=TRUE) &
  !grepl("^ *##", src)),
  "no term count bound remains in the package sources")

## prolfqua still runs when the test is a single term, including the three-level
##   factor that the coefficient-bounded engines cannot handle:

cfg <- cfg4
cfg$test_method <- "prolfqua"
row <- try(f.tune2(state0, cfg), silent=TRUE)
report(!inherits(row, "try-error") && !is.na(row$nhits),
  "f.tune2() runs test_method 'prolfqua' on a three-level factor")

## and a method that can express it still runs on the same joint test:

cfg <- cfg1
cfg$test_method <- "trend"
row <- try(f.tune2(state0, cfg), silent=TRUE)
report(!inherits(row, "try-error") && !is.na(row$nhits),
  "f.tune2() still runs test_method 'trend' on the same joint test")
report(!inherits(row, "try-error") && row$nhits >= 9,
  "f.tune2() reports the interaction genes as hits under 'trend'")

## proda is no longer skipped on the joint test, so tune() sweeps combinations that
##   used to come back NA; this is the sweep-level consequence of the likelihood
##   ratio path in test_proda():

cfg <- cfg1
cfg$test_method <- "proda"
row <- try(f.tune2(state0, cfg), silent=TRUE)
report(!inherits(row, "try-error") && is.data.frame(row) && !is.na(row$nhits) &&
  !is.na(row$ntests), "f.tune2() runs test_method 'proda' on the joint test")
report(!inherits(row, "try-error") && !is.na(row$nhits) && row$nhits >= 9,
  "f.tune2() reports the interaction genes as hits under 'proda'")

###############################################################################
section("reported contrast follows config$reference_levels")

## An engine that builds its own design from the sample annotations re-derives the
##   factor levels by sorting unless it is handed factors, which renames the
##   coefficients: with reference level 'M' h0testr names the tested column sexF,
##   an alphabetical re-derivation names it sexM. test_msqrob() asks
##   msqrob2::hypothesisTest() for its contrast by name, and a name the fit does not
##   have used to come back as a column of NAs with no error at all.
##   Declaring the other level of sex as the reference asks for the opposite
##   contrast, so every estimate should negate exactly while the p-values stay put:

mk_ref <- function(ref) {
  cfg <- cfg0
  cfg$frm <- ~sex + batch
  cfg$test_term <- "sex"
  cfg$reference_levels <- c(sex=ref, batch="b1")
  return(init_state(list(expression=exprs, features=feats, samples=samps), cfg,
    minimal=TRUE))
}

ref_f <- mk_ref("F")
ref_m <- mk_ref("M")

report(identical(levels(ref_m$state$samples$sex), c("M", "F")),
  "reference level 'M' puts M first, so the tested column is sexF")

d_f <- f.design_test_cols(ref_f$state, ref_f$config)
d_m <- f.design_test_cols(ref_m$state, ref_m$config)
report(identical(colnames(d_f$X)[d_f$cols_test], "sexM") &&
  identical(colnames(d_m$X)[d_m$cols_test], "sexF"),
  "f.design_test_cols() names the tested column from the declared reference")

## msqrob and deqms take peptide-level input, proda gene-level; all three are
##   ordered by the engine, so compare on a common key:

by_gene <- function(r, key, est, pv) {
  h <- r$hits
  k <- if(key %in% names(h)) h[[key]] else rownames(h)
  o <- order(k)
  return(list(k=k[o], est=h[[est]][o], p=h[[pv]][o]))
}

run_ref <- function(nm, out) {
  st <- out$state
  cf <- out$config
  if(!f.gene_level_method(nm)) {
    agg <- suppressMessages(combine_features(st, cf))
    st <- agg$state
    cf <- agg$config
  }
  args <- list(st, cf)
  if(nm %in% "proda") args <- list(st, cf, is_log_transformed=TRUE, prior_df=5)
  if(nm %in% "msqrob_agg") args <- list(st, cf, aggregate=TRUE)
  return(try(suppressMessages(do.call(get(paste0("test_",
    sub("_(lmer|agg)$", "", nm))), args)), silent=TRUE))
}

## msqrob first, since it is the one that was broken:

m_f <- run_ref("msqrob", ref_f)
m_m <- run_ref("msqrob", ref_m)
report(!inherits(m_m, "try-error"), "test_msqrob() runs with reference level 'M'")
report(!inherits(m_m, "try-error") && !any(is.na(m_m$hits$logFC)),
  "test_msqrob() reports estimates, not NAs, with reference level 'M'")

## swapping the reference level relabels the fit rather than changing it, so this
##   is exact up to floating point: measured at 1e-14 on 4297 protein groups of
##   rdtc_seer2, holding the random imputation fixed:

if(!inherits(m_f, "try-error") && !inherits(m_m, "try-error")) {
  a <- by_gene(m_f, "gene", "logFC", "pval")
  b <- by_gene(m_m, "gene", "logFC", "pval")
  report(isTRUE(all.equal(a$est, -b$est)),
    "test_msqrob() logFC negates when the reference level is swapped")
  report(isTRUE(all.equal(a$p, b$p)),
    "test_msqrob() p-values are unchanged when the reference level is swapped")
}

## deqms builds its design from state$samples directly, so it was never affected;
##   asserted so that a future change to that path does not go unnoticed:

q_f <- run_ref("deqms", ref_f)
q_m <- run_ref("deqms", ref_m)
if(!inherits(q_f, "try-error") && !inherits(q_m, "try-error")) {
  a <- by_gene(q_f, "gene", "logFC", "P.Value")
  b <- by_gene(q_m, "gene", "logFC", "P.Value")
  report(isTRUE(all.equal(a$est, -b$est)),
    "test_deqms() logFC negates when the reference level is swapped")
  report(isTRUE(all.equal(a$p, b$p)),
    "test_deqms() p-values are unchanged when the reference level is swapped")
} else {
  report(FALSE, "test_deqms() runs under both reference levels")
}

## proda fits its dropout model by iterative optimization, so reparameterizing
##   perturbs the optimizer rather than just relabelling the answer: the estimates
##   negate and the p-values agree, but only to a few significant figures. Measured
##   at max |dp| 1.3e-2 on this fixture, against max |de| of 5.3e-9 for the estimates,
##   so the tolerances below are loose on purpose; a genuine mismatch of contrast is a
##   sign flip or an NA, not a wobble the estimates do not share. It was 1.5e-4 before
##   the peptide counts above were spread for DEqMS's benefit, the optimizer being
##   sensitive to the fixture rather than to anything under test here:

p_f <- run_ref("proda", ref_f)
p_m <- run_ref("proda", ref_m)
if(!inherits(p_f, "try-error") && !inherits(p_m, "try-error")) {
  a <- by_gene(p_f, "name", "diff", "pval")
  b <- by_gene(p_m, "name", "diff", "pval")
  report(isTRUE(all.equal(a$est, -b$est, tolerance=1e-4)),
    "test_proda() diff negates when the reference level is swapped")
  report(max(abs(a$p - b$p)) < 3e-2,
    "test_proda() p-values agree when the reference level is swapped")
} else {
  report(FALSE, "test_proda() runs under both reference levels")
}

## the design test_msqrob() selects from and the design msqrob2::msqrob() fits are
##   built from the same formula, so an interaction coefficient, which is a single
##   coefficient and so within what this engine can test, is found by name:

cfg_int <- cfg0
cfg_int$frm <- ~sex * batch
cfg_int$test_term <- "batch:sex"
cfg_int$reference_levels <- c(sex="F", batch="b1")
out_int <- init_state(list(expression=exprs, features=feats, samples=samps),
  cfg_int, minimal=TRUE)
res <- try(suppressMessages(test_msqrob(out_int$state, out_int$config)), silent=TRUE)
report(!inherits(res, "try-error") && !any(is.na(res$hits$logFC)),
  "test_msqrob() tests an interaction coefficient without name mismatch")

###############################################################################
section("tune_check(): combinations that never ran")

## f.tune2() returns nhits NA and ntests NA for a combination it skipped, and the
##   permuted runs skip identically, so before the fix the NA to 0 substitution made
##   it 0 hits with max1 0, hence fdr 0: the best score in the table, awarded to a
##   combination that was never tested. It belongs below every combination that ran:

tune_dir <- file.path(tempdir(), "tune_check_fixture")
dir.create(tune_dir, showWarnings=FALSE)
suffix <- ".tune.tsv"

mk_tune_row <- function(test, perm, nhits, ntests) {
  data.frame(norm="RLE", nquant=0.5, impute="none", iquant=0, scale=1, span=0.5,
    npcs=3, k=5, test=test, perm=perm, nhits=nhits, ntests=ntests,
    time="00:00:00", stringsAsFactors=FALSE)
}

## 'trend' ran and found 30 of 500, with at most 1 hit in the permutations;
##   'deqms' was skipped, so it has no counts at all:

for(prm in c("0", "1", "2")) {
  is_perm <- !(prm %in% "0")
  rows <- rbind(
    mk_tune_row("trend", prm, if(is_perm) 1 else 30, 500),
    mk_tune_row("deqms", prm, NA, NA)
  )
  utils::write.table(rows, file=file.path(tune_dir, paste0(prm, suffix)),
    sep="\t", quote=FALSE, row.names=FALSE)
}

chk <- try(suppressMessages(tune_check(tune_dir, "", suffix,
  list(log_file=log_file))), silent=TRUE)

report(!inherits(chk, "try-error"), "tune_check() runs on a table with a skipped row")

if(!inherits(chk, "try-error")) {
  report("ntests" %in% names(chk),
    "tune_check() reports ntests, so a skipped row is distinguishable")
  i_skip <- chk$test %in% "deqms"
  i_ran <- chk$test %in% "trend"
  report(sum(i_skip) %in% 1 && is.na(chk$fdr[i_skip]),
    "tune_check() gives a skipped combination no fdr")
  report(sum(i_skip) %in% 1 && chk$ntests[i_skip] %in% 0,
    "tune_check() reports ntests 0 for a skipped combination")
  report(which(i_ran) < which(i_skip),
    "tune_check() ranks the combination that ran above the one that did not")
  report(isTRUE(all.equal(chk$fdr[i_ran], 1 / 30)),
    "tune_check() still computes fdr as max1 / nhits for a combination that ran")
}

unlink(tune_dir, recursive=TRUE)

###############################################################################
cat("\n", strrep("=", 70), "\n", sep="")
cat("\n## passes:", n_pass, "; failures:", n_fail, "; elapsed:",
  round(as.numeric(difftime(Sys.time(), t0, units="secs"))), "s\n")
cat("## log file:", log_file, "\n")

quit(save="no", status=if(n_fail > 0) 1 else 0)
