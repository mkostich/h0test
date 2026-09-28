## Tests for sim_samples() and sim_design(), added to simulate.R on Aug 19. sim1() has no design
##   and sim2() has exactly one two-level factor, so nothing in the suites exercised a third
##   group, a continuous term, an adjustment covariate, an interaction, or a design that cannot
##   be fitted. sim_design() takes the mean structure from a model matrix instead:
##   log2(mean[feat, obs]) = log2(feat_mean[feat]) + x[obs, ] %*% beta[feat, ], leaving sim2()'s
##   intensity, noise and missingness layers and their empirical defaults alone. The assertions
##   below cover what the returned truth claims (that the planted coefficients are what the data
##   actually carries, per SD for a continuous term and per level for a factor), the per-term
##   effect specification, the shape and alignment of the state it returns, the configuration it
##   returns with it, and the argument checking. f.sim0() was generalized from a mean vector to a
##   mean matrix to do this, so one assertion pins the vector path's draws to what they were.

usage <- function(msg=NULL) {
  if(!is.null(msg)) cat("ERROR:", msg, "\n\n", file=stderr())
  cat(
    "Test that h0testr::sim_samples() builds the samples table it documents, balanced or not,",
    "and that h0testr::sim_design() plants the coefficients it reports: that a factor effect",
    "appears as a log2 fold change between one non-reference level and the reference, that a",
    "continuous effect appears as a slope per SD of that covariate, that continuous variables",
    "are scaled for the model matrix but kept raw in the samples table, that effects and",
    "n_genes_signif honor their per-term specification, that the expression matrix, features",
    "and samples it returns line up with each other and run through init_state(), that a",
    "design which cannot be fitted is passed through rather than refused, and that every",
    "documented argument constraint is enforced with a message naming the argument and the",
    "offending value.",
    "",
    "Usage: Rscript test_sim_design.R <r_dir>",
    "",
    "Required positional arguments:",
    "  <r_dir>  Path to the h0testr package R/ source directory; all .R files",
    "             there are sourced (the installed package is not used).",
    "",
    "Optional named arguments: none.",
    "",
    "Needs the packages h0testr itself needs, one assertion running a simulation through",
    "  init_state(), normalize(), filter_state(), impute() and test_h0() to confirm the returned",
    "  configuration matches the returned state.",
    "",
    "Output: one PASS/FAIL line per assertion to stdout, then a count of passes and",
    "  failures and the elapsed time.",
    "",
    "Exit codes: 0 all assertions passed; 1 one or more failed; 2 usage error.",
    "",
    "Examples:",
    "  Rscript test_sim_design.R C:/path/to/h0test/h0testr/R",
    "  Rscript test_sim_design.R ../../h0test/h0testr/R",
    "  Rscript test_sim_design.R C:/path/to/h0testr/R > sim_design.out 2>&1",
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

why <- function(expr) {                  ## the message, for asserting on its content
  out <- try(expr, silent=TRUE)
  if(!inherits(out, "try-error")) return("")
  return(conditionMessage(attr(out, "condition")))
}

quiet <- function(expr) {                ## the pipeline steps narrate; keep the suite readable
  invisible(utils::capture.output(out <- expr))
  return(out)
}

warned <- function(expr) {               ## the value and the warning, for asserting on both
  msg <- ""
  val <- withCallingHandlers(expr, warning=function(w) {
    msg <<- conditionMessage(w)
    invokeRestart("muffleWarning")
  })
  return(list(val=val, msg=msg))
}

## a design-only simulation: no missing values and no dropout, so what comes back is exactly
##   what the mean structure produced:

clean <- function(...) {
  return(sim_design(..., mnar_c0=-Inf, mnar_c1=0, mcar_p=0, p_drop=0))
}

###############################################################################
section("sim_samples() builds the table it documents")

set.seed(1)
samps <- sim_samples(factors=list(sex=c("F", "M"), geno=c("WT", "KO", "HET")))
report(nrow(samps) == 18 && ncol(samps) == 3,
  "a two by three crossing at the default of three per cell is 18 rows")
report(identical(names(samps), c("sample_id", "sex", "geno")),
  "sample_id comes first, then one column per factor")
report(identical(levels(samps$sex), c("F", "M")) && identical(levels(samps$geno),
  c("WT", "KO", "HET")),
  "factor levels keep the order given, so the reference level is the one named first")
report(!any(duplicated(samps$sample_id)) && identical(samps$sample_id[1], "samp1"),
  "sample ids are distinct and start at samp1")
cells <- unique(samps[, c("sex", "geno")])
report(identical(as.character(cells$sex), rep(c("F", "M"), 3)) &&
  identical(as.character(cells$geno), rep(c("WT", "KO", "HET"), each=2)),
  "cells come in expand.grid() order, the first factor fastest, which n_per_cell is read in")
report(identical(as.character(samps$sex[1:4]), c("F", "F", "F", "M")),
  "and a cell's samples are consecutive, so the ordering is systematic, not interleaved")

samps <- sim_samples(factors=list(sex=c("F", "M")), n_per_cell=c(5, 2))
report(sum(samps$sex == "F") == 5 && sum(samps$sex == "M") == 2,
  "n_per_cell as a vector unbalances the design deliberately")

samps <- sim_samples(factors=list(sex=c("F", "M"), geno=c("WT", "KO")),
  n_per_cell=c(3, 3, 3, 0))
report(nrow(samps) == 9 && nrow(unique(samps[, c("sex", "geno")])) == 3,
  "a cell given 0 is left out entirely, which is one way to reach an empty cell")
report(identical(levels(samps$geno), c("WT", "KO")),
  "and the level of the empty cell is still declared, so the design stays rank deficient")

set.seed(2)
samps <- sim_samples(factors=list(sex=c("F", "M")), covariates=list(age=c(20, 60)), n_per_cell=8)
report(is.numeric(samps$age) && all(samps$age >= 20 & samps$age <= 60),
  "a length 2 covariate range is drawn uniformly inside that range")

samps <- sim_samples(covariates=c("age", "intake"), n=200)
report(nrow(samps) == 200 && all(c("age", "intake") %in% names(samps)),
  "a character vector of covariate names is shorthand for the default draw")
report(abs(mean(samps$age)) < 0.3 && abs(stats::sd(samps$age) - 1) < 0.2,
  "and that default is standard normal, which unlike a uniform draw has tails")

samps <- sim_samples(covariates=list(age=function(n) rep(seq(20, 60, length.out=5), n / 5)), n=10)
report(nrow(samps) == 10 && length(unique(samps$age)) == 5,
  "a function of n draws a covariate however the caller likes")

report(!("sample_id" %in% c()) && identical(rownames(sim_samples(covariates="age", n=3)),
  c("1", "2", "3")),
  "rownames are reset, so deleting rows to unbalance leaves no gaps behind")

###############################################################################
section("sim_samples() argument checking")

report(grepl("both were NULL", why(sim_samples())),
  "factors and covariates cannot both be NULL")
report(grepl("n is required", why(sim_samples(covariates="age"))),
  "n is required when there are no cells to count")
report(grepl("not both", why(sim_samples(factors=list(sex=c("F", "M")), n=6))),
  "and refused when there are, n_per_cell determining the count")
report(grepl("distinct name per factor", why(sim_samples(factors=list(c("F", "M"))))),
  "factors must be named")
report(grepl("distinct non-missing levels", why(sim_samples(factors=list(sex=c("F", "F"))))),
  "and its levels must be distinct")
report(grepl("whole numbers >= 0", why(sim_samples(factors=list(sex=c("F", "M")),
  n_per_cell=-1))),
  "n_per_cell cannot be negative")
report(grepl("whole numbers >= 0", why(sim_samples(factors=list(sex=c("F", "M")),
  n_per_cell=2.5))),
  "nor fractional")
report(grepl("length 1 or one per cell", why(sim_samples(factors=list(sex=c("F", "M")),
  n_per_cell=c(1, 2, 3)))),
  "and a length that is neither 1 nor one per cell is refused, naming both counts")
report(grepl("no samples", why(sim_samples(factors=list(sex=c("F", "M")), n_per_cell=0))),
  "every cell at 0 would leave no samples at all")
report(grepl("already in use", why(sim_samples(factors=list(sex=c("F", "M")),
  covariates=list(sex=c(0, 1))))),
  "a covariate cannot take a name a factor already has")
report(grepl("already in use", why(sim_samples(covariates=list(sample_id=c(0, 1)), n=4))),
  "nor the name sample_id")
report(grepl("increasing length 2 numeric range", why(sim_samples(covariates=list(age=c(60, 20)),
  n=4))),
  "a covariate range must increase, a constant covariate carrying no per-SD effect")
report(grepl("must return 4 finite numeric", why(sim_samples(covariates=list(age=function(n) 1),
  n=4))),
  "a covariate function returning the wrong length is refused, naming the length wanted")

###############################################################################
section("sim_design() returns a state whose three pieces line up")

set.seed(3)
samps <- sim_samples(factors=list(sex=c("F", "M")), covariates=list(age=c(20, 60)), n_per_cell=4)
sim <- clean(samps, frm=~sex + age, test_term="sex", n_genes=10, n_genes_signif=c(sex=3, age=5))
report(is.matrix(sim$state$expression) && all(dim(sim$state$expression) == c(10, 8)),
  "one row per gene and one column per observation")
report(identical(rownames(sim$state$expression), sim$state$features$feature_id),
  "the features table is in the same order as the matrix rows")
report(identical(colnames(sim$state$expression), sim$state$samples$observation_id),
  "and the samples table is in the same order as the matrix columns")
report(all(sim$state$expression >= 1) && all(sim$state$expression %% 1 == 0),
  "the matrix is whole numbers of at least 1, ceiling() keeping a positive draw positive")
report(identical(names(sim$state$samples)[1:2], c("observation_id", "sample_id")),
  "the identifier columns come first in the samples table")
report(identical(sim$state$samples$observation_id, sim$state$samples$sample_id),
  "with one observation per sample when there are no technical replicates")
report(identical(sim$state$features$gene_id, as.character(sim$feat_gene)),
  "gene_id agrees with the feat_gene map")
report(identical(dim(sim$truth), c(10L, 2L)) &&
  identical(colnames(sim$truth), c("sexM", "age")),
  "truth is genes by model matrix columns, the intercept excluded")
report(identical(rownames(sim$x), samps$sample_id) &&
  identical(colnames(sim$x), colnames(sim$truth)),
  "and x is samples by those same columns")

set.seed(4)
samps <- sim_samples(factors=list(geno=c("WT", "KO", "HET")), n_per_cell=2)
sim <- clean(samps, frm=~geno, test_term="geno", n_genes=4, peps_per_gene=3, reps_per_sample=2)
report(all(dim(sim$state$expression) == c(12, 12)),
  "peps_per_gene multiplies the rows and reps_per_sample the columns")
report(identical(sim$state$samples$observation_id[1:2], c("samp1_rep1", "samp1_rep2")) &&
  identical(sim$state$samples$sample_id[1:2], c("samp1", "samp1")),
  "a technical replicate gets its own observation id and repeats its sample id")
report(all(table(sim$state$samples$sample_id) == 2) && nrow(sim$state$samples) == 12,
  "so sample_id is what identifies the biological sample, as config$sample_id_col says")
report(all(table(sim$feat_gene) == 3) && length(unique(sim$feat_gene)) == 4,
  "every gene keeps its three features when p_drop is 0")
report(identical(names(sim$feat_gene), rownames(sim$state$expression)),
  "and feat_gene is named by feature")
report(identical(names(sim$feat_mean), rownames(sim$state$expression)) &&
  identical(names(sim$feat_cv), rownames(sim$state$expression)),
  "as are the returned parameter means and CVs")

###############################################################################
section("the planted coefficients are what the data actually carries")

## a large simulation with a small CV, so the observed effect is the planted one to two decimals
##   rather than merely in the right direction:

set.seed(5)
samps <- sim_samples(factors=list(sex=c("F", "M")), n_per_cell=40)
sim <- clean(samps, frm=~sex, test_term="sex", n_genes=400, n_genes_signif=200, effects=1,
  log_cv_mean=-2, log_cv_sd=0.2)
mat <- log2(sim$state$expression)
i_m <- sim$state$samples$sex == "M"
obs <- rowMeans(mat[, i_m]) - rowMeans(mat[, !i_m])
tru <- sim$truth[rownames(mat), "sexM"]
report(abs(mean(obs[tru == 1]) - 1) < 0.05 && abs(mean(obs[tru == -1]) + 1) < 0.05,
  "a factor effect is a log2 fold change of the non-reference level against the reference")
report(abs(mean(obs[tru == 0])) < 0.05,
  "and a gene with no planted effect has no fold change")
report(stats::cor(obs, tru) > 0.95,
  "so the observed fold changes track the reported truth gene by gene")
report(sum(tru > 0) > 60 && sum(tru < 0) > 60,
  "increases and decreases are drawn with roughly equal chance")

set.seed(6)
samps <- sim_samples(covariates=list(age=c(20, 60)), n=200)
sim <- clean(samps, frm=~age, test_term="age", n_genes=400, n_genes_signif=200, effects=1,
  log_cv_mean=-2, log_cv_sd=0.2)
mat <- log2(sim$state$expression)
age_z <- as.numeric(scale(sim$state$samples$age))
slope <- apply(mat, 1, function(v) stats::coef(stats::lm(v ~ age_z))[2])
tru <- sim$truth[rownames(mat), "age"]
report(abs(mean(slope[tru == 1]) - 1) < 0.05 && abs(mean(slope[tru == -1]) + 1) < 0.05,
  "a continuous effect is a log2 fold change per SD of that covariate")
report(abs(mean(slope[tru == 0])) < 0.05,
  "and a covariate with no planted effect has no slope")

report(all(sim$state$samples$age >= 20) && mean(sim$state$samples$age) > 30,
  "the samples table keeps the covariate on the scale it was drawn on")
report(abs(mean(sim$x[, "age"])) < 1e-8 && abs(stats::sd(sim$x[, "age"]) - 1) < 1e-8,
  "while the model matrix holds it centered and scaled, which is what per SD means")

## location matters: an uncentered covariate would carry the intensity distribution off
##   log_m_mean, and with it the global missingness rate:

set.seed(7)
s1 <- sim_samples(covariates=list(age=c(20, 60)), n=30)
s2 <- s1
s2$age <- s2$age + 1000
set.seed(8)
a <- sim_design(s1, frm=~age, test_term="age", n_genes=300, n_genes_signif=150, effects=1)
set.seed(8)
b <- sim_design(s2, frm=~age, test_term="age", n_genes=300, n_genes_signif=150, effects=1)
report(identical(sum(is.na(a$state$expression)), sum(is.na(b$state$expression))),
  "so shifting a covariate by 1000 leaves the missingness rate exactly where it was")

###############################################################################
section("effects and n_genes_signif, per term")

set.seed(9)
samps <- sim_samples(factors=list(sex=c("F", "M")), covariates=list(age=c(20, 60)), n_per_cell=4)
sim <- clean(samps, frm=~sex + age, test_term="sex", n_genes=20, n_genes_signif=c(sex=3, age=7))
report(sum(sim$truth[, "sexM"] != 0) == 3 && sum(sim$truth[, "age"] != 0) == 7,
  "each term gets the number of significant genes it was given, not a shared number")

sim <- clean(samps, frm=~sex + age, test_term="sex", n_genes=20, n_genes_signif=c(sex=3),
  effects=c(sex=2))
report(all(sim$truth[, "age"] == 0) && all(abs(sim$truth[, "sexM"])[sim$truth[, "sexM"] != 0] == 2),
  "a term the caller did not name gets no effect, and the named one gets the size given")

sim <- clean(samps, frm=~sex + age, test_term="sex", n_genes=20, n_genes_signif=4, effects=1.5)
report(sum(sim$truth[, "sexM"] != 0) == 4 && sum(sim$truth[, "age"] != 0) == 4 &&
  all(abs(sim$truth[sim$truth != 0]) == 1.5),
  "a scalar applies to every term of frm")
report(any(rowSums(sim$truth != 0) == 2),
  "and significant genes are drawn per term, so one gene can carry two terms' effects")

## the two per-term specifications are read the same way, so naming one term in either zeroes
##   the rest. Genes without an effect size is the asymmetry worth a word, the reverse being how
##   a term is deliberately left alone:

out <- warned(clean(samps, frm=~sex + age, test_term="sex", n_genes=20,
  n_genes_signif=c(sex=3, age=5), effects=c(age=1)))
report(all(out$val$truth[, "sexM"] == 0) && sum(out$val$truth[, "age"] != 0) == 5 &&
  grepl("^sim_design: no effect planted on term sex", out$msg),
  "a term given genes but an effect size of 0 gets no effect, and says so")
report(grepl("n_genes_signif asks for 3", out$msg) && grepl("effects: sex=0, age=1", out$msg),
  "the warning naming the count asked for and what effects actually holds for every term")

out <- warned(clean(samps, frm=~sex + age, test_term="sex", n_genes=20,
  n_genes_signif=c(sex=3), effects=1))
report(sum(out$val$truth[, "sexM"] != 0) == 3 && all(out$val$truth[, "age"] == 0) &&
  !nzchar(out$msg),
  "while an effect size with no genes to spend it on is quiet, that being the usual idiom")

out <- warned(clean(samps, frm=~sex + age, test_term="sex", n_genes=20, n_genes_signif=4,
  effects=0))
report(all(out$val$truth == 0) && grepl("no effect planted on term", out$msg),
  "and effects of 0 for every term warns per term rather than returning an empty truth quietly")

report(grepl("sexx", why(clean(samps, frm=~sex + age, test_term="sex", n_genes=8,
  effects=c(sexx=1)))) &&
  grepl("terms of frm: sex, age", why(clean(samps, frm=~sex + age, test_term="sex", n_genes=8,
  effects=c(sexx=1)))),
  "a misspelled term is refused, naming it and the terms of frm, not silently ignored")
report(grepl("more than once", why(clean(samps, frm=~sex + age, test_term="sex", n_genes=8,
  effects=c(sex=1, sex=2)))),
  "and a term named twice is refused")
report(grepl("scalar or a vector named", why(clean(samps, frm=~sex + age, test_term="sex",
  n_genes=8, effects=c(1, 2)))),
  "an unnamed vector of length above 1 is refused, there being no order to read it in")
report(grepl("n_genes_signif", why(clean(samps, frm=~sex, test_term="sex", n_genes=8,
  n_genes_signif=9))) &&
  grepl("\\[0, 8\\]", why(clean(samps, frm=~sex, test_term="sex", n_genes=8, n_genes_signif=9))),
  "more significant genes than genes is refused, naming the bound")

set.seed(10)
samps <- sim_samples(factors=list(geno=c("WT", "KO", "HET")), n_per_cell=3)
sim <- clean(samps, frm=~geno, test_term="geno", n_genes=30, n_genes_signif=12, effects=1)
report(identical(colnames(sim$truth), c("genoKO", "genoHET")),
  "a three level factor gives one column per non-reference level")
report(all(rowSums(sim$truth != 0) <= 1) && sum(rowSums(sim$truth != 0)) == 12,
  "and each significant gene is changed in exactly one of those levels")
report(sum(sim$truth[, "genoKO"] != 0) > 0 && sum(sim$truth[, "genoHET"] != 0) > 0,
  "with the level chosen at random, so both levels are used across 12 genes")

mat <- log2(sim$state$expression)
g <- sim$state$samples$geno
i_gene <- which(sim$truth[, "genoKO"] != 0)[1]
lfc_ko <- mean(mat[i_gene, g == "KO"]) - mean(mat[i_gene, g == "WT"])
lfc_het <- mean(mat[i_gene, g == "HET"]) - mean(mat[i_gene, g == "WT"])
report(abs(lfc_ko) > 0.5 && abs(lfc_het) < abs(lfc_ko),
  "so a gene changed in KO moves against WT and its untouched level does not")

set.seed(11)
samps <- sim_samples(factors=list(sex=c("F", "M"), geno=c("WT", "KO")), n_per_cell=4)
sim <- clean(samps, frm=~sex * geno, test_term="sex:geno", n_genes=20,
  n_genes_signif=c("sex:geno"=5), effects=1)
report(identical(colnames(sim$truth), c("sexM", "genoKO", "sexM:genoKO")),
  "an interaction gets its own model matrix column")
report(sum(sim$truth[, "sexM:genoKO"] != 0) == 5 && all(sim$truth[, "sexM"] == 0),
  "and an effect can be planted on the interaction alone")

## a continuous by factor interaction is the case where per SD scaling could have gone wrong:
##   covariates are scaled on the samples table, before model.matrix() builds any term from them,
##   so the interaction column is sexM * age_z and its coefficient is the difference of the two
##   groups' per SD slopes. Large n and a small CV, so the recovered value is the planted one:

set.seed(11)
samps <- sim_samples(factors=list(sex=c("F", "M")), covariates=list(age=c(20, 60)), n_per_cell=100)
sim <- clean(samps, frm=~sex * age, test_term="sex:age", n_genes=400,
  n_genes_signif=c("sex:age"=200), effects=1, log_cv_mean=-2, log_cv_sd=0.2)
mat <- log2(sim$state$expression)
age_z <- as.numeric(scale(sim$state$samples$age))
i_m <- sim$state$samples$sex == "M"
f.slope <- function(v, i) return(stats::coef(stats::lm(v[i] ~ age_z[i]))[2])
d_slope <- apply(mat, 1, function(v) f.slope(v, i_m) - f.slope(v, !i_m))
tru <- sim$truth[rownames(mat), "sexM:age"]
report(all(sim$truth[, "sexM"] == 0) && all(sim$truth[, "age"] == 0) && sum(tru != 0) == 200,
  "an effect on a continuous by factor interaction leaves both main effect columns alone")
report(abs(mean(d_slope[tru == 1]) - 1) < 0.05 && abs(mean(d_slope[tru == -1]) + 1) < 0.05,
  "and it recovers as the difference of the two groups' per SD slopes, in SD units not raw units")
report(abs(mean(d_slope[tru == 0])) < 0.05,
  "while a gene with no interaction effect has the same slope in both groups")

###############################################################################
section("the configuration returned with the state")

set.seed(12)
samps <- sim_samples(factors=list(sex=c("F", "M")), covariates=list(age=c(20, 60)), n_per_cell=6)
sim <- sim_design(samps, frm=~sex + age, test_term="sex", n_genes=300,
  n_genes_signif=c(sex=40, age=80), effects=1)
cfg <- sim$config
report(identical(cfg$test_term, "sex") && identical(deparse(cfg$frm), deparse(~sex + age)),
  "frm and test_term are set to what was simulated")
report(identical(cfg$feat_id_col, "feature_id") && identical(cfg$gene_id_col, "gene_id") &&
  identical(cfg$obs_id_col, "observation_id") && identical(cfg$sample_id_col, "sample_id"),
  "the four id columns name the columns actually present")
report(identical(cfg$reference_levels, c(sex="F")),
  "reference_levels names the factor's first level, and only the factor")
report(isFALSE(cfg$save_state),
  "save_state is FALSE, so a simulation leaves no files behind")
report(identical(cfg$n_features_min, 1),
  "n_features_min is 1, its default of 1000 being meant for a real dataset")
report(is.null(cfg$covariate_types) || !length(cfg$covariate_types),
  "covariate_types is left unset, being init_state()'s to derive")

out <- quiet(init_state(sim$state, cfg))
report(!inherits(out, "try-error") && is.matrix(out$state$expression),
  "the state and the configuration go through init_state() as they come")
report(identical(out$config$covariate_types[["sex"]], "factor") &&
  identical(out$config$covariate_types[["age"]], "numeric"),
  "and init_state() reads sex as a factor and age as continuous, unprompted")

out$state <- quiet(add_filter_stats(out$state, out$config))
out$state <- quiet(prefilter(out$state, out$config))
out <- quiet(normalize(out$state, out$config))
out <- quiet(filter_state(out$state, out$config))
out <- quiet(impute(out$state, out$config))
res <- quiet(test_h0(out$state, out$config))
tbl <- res$standard
i_true <- tbl$feature %in% rownames(sim$truth)[sim$truth[, "sexM"] != 0]
report(nrow(tbl) > 100 && all(c("logfc", "pval", "adj_pval") %in% names(tbl)),
  "and all the way through test_h0(), so the returned config really does match the state")
report(stats::median(tbl$pval[i_true]) < stats::median(tbl$pval[!i_true]),
  "the planted genes come out with smaller p values than the untouched ones")
report(mean(sign(tbl$logfc[i_true]) == sign(sim$truth[tbl$feature[i_true], "sexM"])) > 0.8,
  "and test_h0()'s logfc agrees in sign with the planted coefficient")

set.seed(13)
samps <- sim_samples(covariates=list(age=c(20, 60)), n=8)
sim <- clean(samps, frm=~age, test_term="age", n_genes=6, n_genes_signif=2)
report(is.character(sim$config$reference_levels) && !length(sim$config$reference_levels),
  "a design of continuous terms only leaves reference_levels empty rather than absent")
report(!threw(quiet(init_state(sim$state, sim$config))),
  "and init_state() accepts that")

###############################################################################
section("designs that cannot be fitted are passed through, not refused")

set.seed(14)
samps <- sim_samples(factors=list(sex=c("F", "M"), geno=c("WT", "KO")),
  n_per_cell=c(3, 3, 3, 0))
out <- warned(clean(samps, frm=~sex * geno, test_term="sex:geno", n_genes=8, n_genes_signif=2))
sim <- out$val
report(is.matrix(sim$state$expression) && nrow(sim$state$expression) == 8,
  "an empty cell makes the design rank deficient, which sim_design() simulates anyway")
report(qr(sim$x)$rank < ncol(sim$x),
  "the rank deficiency being real, and left for test_h0() to have an opinion about")
report(all(sim$x[, "sexM:genoKO"] == 0) && all(sim$truth[, "sexM:genoKO"] == 0) &&
  sum(sim$truth[, "sexM"] != 0) == 2 && grepl("no effect planted on term sex:geno", out$msg),
  "an empty cell being also why the interaction column is constant, so nothing is planted on it")

samps <- sim_samples(factors=list(sex=c("F", "M")), n_per_cell=4)
samps <- samps[samps$sex == "F", ]
out <- warned(clean(samps, frm=~sex, test_term="sex", n_genes=8, n_genes_signif=2))
report(is.list(out$val) && nrow(out$val$state$expression) == 8,
  "a factor with only one level present is likewise simulated, not refused")

## but the effect goes on a column the data can express, never on the all zero dummy of a level
##   that was declared and then not sampled, which would be a truth the data does not carry:

report(all(out$val$truth == 0) && grepl("^sim_design: no effect planted on term sex", out$msg) &&
  grepl("sexM", out$msg),
  "no effect is planted when every column of the term is constant, and it says so")

set.seed(141)
samps <- sim_samples(factors=list(geno=c("WT", "KO", "HET")), n_per_cell=c(6, 6, 0))
out <- warned(clean(samps, frm=~geno, test_term="geno", n_genes=40, n_genes_signif=20, effects=2,
  log_cv_mean=-2, log_cv_sd=0.2))
sim <- out$val
report(identical(colnames(sim$truth), c("genoKO", "genoHET")) && all(sim$x[, "genoHET"] == 0),
  "an unsampled level of a three level factor still gets its own all zero column")
report(sum(sim$truth[, "genoKO"] != 0) == 20 && all(sim$truth[, "genoHET"] == 0) &&
  !nzchar(out$msg),
  "so all 20 significant genes take their effect on the level that varies, with no warning")

mat <- log2(sim$state$expression)
g <- sim$state$samples$geno
lfc <- rowMeans(mat[, g == "KO"]) - rowMeans(mat[, g == "WT"])
tru <- sim$truth[rownames(mat), "genoKO"]
report(abs(mean(abs(lfc[tru != 0])) - 2) < 0.1 && abs(mean(lfc[tru == 0])) < 0.1,
  "and the data carries every one of them, which is what the constant column check buys")

samps <- sim_samples(factors=list(sex=c("F", "M")), covariates=list(age=c(20, 60)), n_per_cell=4)
samps <- samps[-c(2, 5), ]
report(nrow(clean(samps, frm=~sex + age, test_term="sex", n_genes=8)$state$samples) == 6,
  "deleting rows to unbalance a design needs nothing else to be told about it")

###############################################################################
section("sim_design() argument checking")

set.seed(15)
samps <- sim_samples(factors=list(sex=c("F", "M")), covariates=list(age=c(20, 60)), n_per_cell=3)

report(grepl("samps must be a data.frame", why(clean(as.matrix(1:4), frm=~sex,
  test_term="sex", n_genes=8))),
  "samps must be a data.frame, its class named")
report(grepl("at least one row", why(clean(samps[0, ], frm=~sex, test_term="sex", n_genes=8))),
  "and must have rows")
report(grepl("frm must be a formula", why(clean(samps, frm="~sex", test_term="sex", n_genes=8))),
  "frm must be a formula, not a string that looks like one")
report(grepl("at least one term", why(clean(samps, frm=~1, test_term="sex", n_genes=8))),
  "and must carry a term")
report(grepl("not a term of frm", why(clean(samps, frm=~sex, test_term="age", n_genes=8))),
  "test_term must be one of frm's terms")
report(grepl("missing from samps", why(clean(samps, frm=~sex + weight, test_term="sex",
  n_genes=8))) &&
  grepl("weight", why(clean(samps, frm=~sex + weight, test_term="sex", n_genes=8))),
  "and every variable of frm must be a column of samps, the missing one named")

report(grepl("n_genes", why(clean(samps, frm=~sex, test_term="sex", n_genes=1))),
  "n_genes below 2 is refused")
report(grepl("peps_per_gene", why(clean(samps, frm=~sex, test_term="sex", n_genes=8,
  peps_per_gene=0))),
  "as is a gene with no features")
report(grepl("reps_per_sample", why(clean(samps, frm=~sex, test_term="sex", n_genes=8,
  reps_per_sample=1.5))),
  "and a fractional number of technical replicates")
## ported from test_simulate.R, where these four f.chk_num() and f.chk_term_vec() call sites
##   were probed through sim2(). sim_design() has the same four, and nothing else here
##   reaches them; the other six of that file's checks are already covered above:

report(grepl("cv_reps", why(clean(samps, frm=~sex, test_term="sex", n_genes=8, cv_reps=-1))),
  "cv_reps cannot be negative")
report(grepl("p_drop", why(sim_design(samps, frm=~sex, test_term="sex", n_genes=8,
  p_drop=1.5))) &&
  grepl("p_drop", why(sim_design(samps, frm=~sex, test_term="sex", n_genes=8, p_drop=-0.1))),
  "and p_drop stays a probability at both ends")
report(!threw(sim_design(samps, frm=~sex, test_term="sex", n_genes=8, peps_per_gene=2,
  p_drop=1)),
  "including p_drop of exactly 1, which keeps one feature per gene rather than none")
report(grepl("effects", why(clean(samps, frm=~sex, test_term="sex", n_genes=8, effects=NA))),
  "effects must be present and finite, a missing effect size being no effect size")

report(!threw(clean(samps, frm=~sex, test_term="sex", n_genes=8, log_m_mean=0)) &&
  !threw(clean(samps, frm=~sex, test_term="sex", n_genes=8, log_m_sd=0)),
  "either of log_m_mean and log_m_sd may be 0 on its own")
report(grepl("cannot both be 0", why(clean(samps, frm=~sex, test_term="sex", n_genes=8,
  log_m_mean=0, log_m_sd=0))) &&
  grepl("log_m_mean: 0; log_m_sd: 0", why(clean(samps, frm=~sex, test_term="sex", n_genes=8,
  log_m_mean=0, log_m_sd=0))),
  "but not both, and sim_design refuses that itself rather than by way of f.sim_rnorm_pos()")

report(grepl("mnar_off must lie in \\(0, 1\\)", why(sim_design(samps, frm=~sex, test_term="sex",
  n_genes=8, mnar_off=1))),
  "mnar_off keeps the open interval its own documentation gives it")
report(grepl("mcar_p", why(sim_design(samps, frm=~sex, test_term="sex", n_genes=8, mcar_p=1.5))),
  "and mcar_p stays a probability")
report(!threw(sim_design(samps, frm=~sex, test_term="sex", n_genes=8, mnar_c0=-Inf, mnar_c1=0)),
  "mnar_c0 of -Inf is allowed, that being the documented way to turn MNAR off")
report(grepl("must be finite", why(sim_design(samps, frm=~sex, test_term="sex", n_genes=8,
  mnar_c0=Inf))),
  "while +Inf is refused, p(mnar) being NaN there")

bad <- samps
bad$age <- 50
report(grepl("no variance", why(clean(bad, frm=~sex + age, test_term="sex", n_genes=8))) &&
  grepl("age", why(clean(bad, frm=~sex + age, test_term="sex", n_genes=8))),
  "a covariate of one repeated value cannot carry a per-SD effect, and is named")

bad <- samps
bad$age[1] <- NA
report(grepl("must be finite", why(clean(bad, frm=~age, test_term="age", n_genes=8))),
  "a missing covariate value is refused here rather than dropping a row in model.matrix()")

bad <- samps
bad$sex[1] <- NA
report(grepl("non-missing", why(clean(bad, frm=~sex, test_term="sex", n_genes=8))),
  "as is a missing factor value")

bad <- samps
bad$z <- complex(real=1:6, imaginary=1)
report(grepl("must be numeric, logical, character or factor",
  why(clean(bad, frm=~z, test_term="z", n_genes=8))),
  "and a variable of a type no design can use is named with its class")

bad <- samps
bad$sample_id <- "same"
report(grepl("must be unique", why(clean(bad, frm=~sex, test_term="sex", n_genes=8))),
  "sample ids that repeat are refused, an observation needing to be identifiable")

###############################################################################
section("treatment contrasts, without which truth means something else")

## a factor coefficient is a log2 fold change against the reference level only under treatment
##   contrasts. An ordered factor gets contr.poly from model.matrix() whatever options("contrasts")
##   says, a contrasts attribute on the variable is honored the same way, and the option itself
##   can be set session-wide; each would leave truth's column names, and the reference level
##   config$reference_levels names, describing something the data does not carry. Refusing is
##   deliberate rather than passing contrasts.arg: test_h0() reads the same option, so pinning it in
##   the simulator alone would score truth in one parameterization against a fit in another:

samps_c <- sim_samples(factors=list(geno=c("WT", "KO", "HET")), n_per_cell=3)

report(!threw(clean(samps_c, frm=~geno, test_term="geno", n_genes=6)),
  "an unordered factor under the default contrasts is simulated as before")

bad <- samps_c
bad$geno <- factor(bad$geno, levels=c("WT", "KO", "HET"), ordered=TRUE)
report(grepl("^sim_design: variable geno is an ordered factor",
  why(clean(bad, frm=~geno, test_term="geno", n_genes=6))),
  "an ordered factor is refused up front, naming the variable")
report(grepl("polynomial contrasts", why(clean(bad, frm=~geno, test_term="geno", n_genes=6))) &&
  grepl("factor\\(levels=c\\(\"WT\", \"KO\", \"HET\"\\)\\)",
  why(clean(bad, frm=~geno, test_term="geno", n_genes=6))),
  "and says why, with the declaration to use instead")

bad <- samps_c
bad$geno <- stats::C(bad$geno, "contr.sum")
report(grepl("carries a contrasts attribute",
  why(clean(bad, frm=~geno, test_term="geno", n_genes=6))),
  "a contrasts attribute on the variable is refused too, model.matrix() honoring it")

old_contr <- getOption("contrasts")
options(contrasts=c("contr.sum", "contr.poly"))
msg_f <- why(clean(samps_c, frm=~geno, test_term="geno", n_genes=6))
samps_n <- sim_samples(covariates=list(age=c(20, 60)), n=8)
threw_n <- threw(clean(samps_n, frm=~age, test_term="age", n_genes=6))
options(contrasts=old_contr)

report(grepl("options\\(\"contrasts\"\\) gives contr.sum", msg_f) &&
  grepl("contr.treatment", msg_f),
  "a session-wide contrasts option other than contr.treatment is refused, and named")
report(!threw_n,
  "but only where there is a factor to interpret: an all-continuous design is unaffected")
report(identical(getOption("contrasts"), old_contr),
  "and the option is put back, so the assertions after this one see the default")

###############################################################################
section("dropout leaves every gene a feature")

## f.pep_drop() draws per feature and floors per gene, so the heaviest dropout there is leaves
##   one feature of each gene rather than one feature altogether. The floor used to be on the
##   whole matrix, and a gene could lose every feature and disappear from truth with it:

set.seed(16)
samps <- sim_samples(factors=list(sex=c("F", "M")), n_per_cell=3)
sim <- sim_design(samps, frm=~sex, test_term="sex", n_genes=6, n_genes_signif=3,
  peps_per_gene=2, p_drop=1, mnar_c0=-Inf, mnar_c1=0, mcar_p=0)
report(nrow(sim$state$expression) == 6,
  "p_drop of 1 leaves one feature of every gene rather than one feature overall")
report(nrow(sim$truth) == 6 && all(rownames(sim$truth) %in% sim$feat_gene) &&
  all(sim$feat_gene %in% rownames(sim$truth)),
  "and truth keeps every gene, none of them having lost all of its features")
report(length(sim$feat_mean) == 6 && length(sim$feat_cv) == 6 &&
  nrow(sim$state$features) == 6 && all(table(sim$feat_gene) == 1),
  "as do the parameters and the features table, at one feature each")

set.seed(17)
sim <- sim_design(samps, frm=~sex, test_term="sex", n_genes=40, n_genes_signif=10,
  peps_per_gene=4, p_drop=0.5, mnar_c0=-Inf, mnar_c1=0, mcar_p=0)
report(nrow(sim$state$expression) < 160 && nrow(sim$state$expression) > 40,
  "a partial dropout leaves a heterogeneous number of features per gene")
report(all(rownames(sim$truth) %in% sim$feat_gene) &&
  all(sim$feat_gene %in% rownames(sim$truth)),
  "and truth holds exactly the genes that still have a feature")
report(nrow(sim$truth) == 40 && all(table(sim$feat_gene) >= 1) &&
  all(table(sim$feat_gene) <= 4) && length(unique(table(sim$feat_gene))) > 1,
  "which is all forty of them, each with one to peps_per_gene of them, and not all alike")

###############################################################################
section("f.sim0()'s vector path draws what it always drew")

## f.sim0() was generalized to accept a matrix of per-observation means, which is the whole of
##   how a design is imposed. It has always gone observation by observation, so a matrix whose
##   columns are all the same vector has to consume the RNG exactly as the vector call does; if
##   it does not, sim1() and sim2() no longer reproduce their seeded fixtures:

m <- c(10, 20, 30, 40)
s <- c(1, 2, 3, 4)
set.seed(101)
a <- f.sim0(n_obs=5, feat_means=m, feat_sds=s)
set.seed(101)
b <- f.sim0(n_obs=5, feat_means=matrix(m, nrow=4, ncol=5),
  feat_sds=matrix(s, nrow=4, ncol=5))
report(identical(dim(a), dim(b)) && isTRUE(all.equal(unname(a), unname(b))),
  "the matrix path with constant columns draws the same values as the vector path")

set.seed(102)
a <- sim1(n_obs=4, n_feats=6)
set.seed(102)
b <- sim1(n_obs=4, n_feats=6)
report(isTRUE(all.equal(a, b)),
  "and sim1() is still reproducible from a seed")

###############################################################################
cat("\n## passes:", n_pass, "; failures:", n_fail, "; elapsed:",
  round(as.numeric(difftime(Sys.time(), t0, units="secs"))), "s\n")

quit(save="no", status=if(n_fail > 0) 1 else 0)
