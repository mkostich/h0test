## Tests the joint Wald test test_msqrob() now runs when config$test_term resolves to
##   more than one design matrix column. msqrob2::hypothesisTest() returns one table
##   per contrast, so those tests used to be refused; the fitted msqrob2 StatModel
##   carries the coefficients, the unscaled covariance and the moderated variance and
##   posterior degrees of freedom, so the statistic is computed from the fit instead.
##   The single column path is unchanged and still goes through
##   msqrob2::hypothesisTest(); the central assertion here is that the new code
##   reproduces that function exactly when handed one column, which is what justifies
##   referring the statistic to an F on dfPosterior rather than to a chi-square.

usage <- function(msg=NULL) {
  if(!is.null(msg)) cat("ERROR:", msg, "\n\n", file=stderr())
  cat(
    "Test the joint Wald test h0testr::test_msqrob() runs when config$test_term",
    "spans more than one design matrix column. Checks that the new code reproduces",
    "msqrob2::hypothesisTest() when handed a single column, that the joint statistic",
    "is an F on the posterior degrees of freedom and is more conservative than the",
    "chi-square it would otherwise be referred to, that the previously refused tests",
    "(a factor with more than two levels, a variable inside an interaction) now run",
    "and rank the known effects first, that logfc holds the total swing of the",
    "msqrob2 coefficients, and that a coefficient the fit does not carry gives NA",
    "rather than an error. Also checks that a continuous config$test_term resolves to",
    "a single column and is recovered as a slope per unit of the covariate, and that",
    "no engine is capped at one coefficient",
    "any more, test_deqms() having been the last; see test_deqms.R.",
    "",
    "Usage: Rscript test_msqrob_joint.R <r_dir>",
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
    "  3 msqrob2 not installed.",
    "",
    "Examples:",
    "  Rscript test_msqrob_joint.R C:/path/to/h0test/h0testr/R",
    "  Rscript test_msqrob_joint.R ../../h0test/h0testr/R",
    "  Rscript test_msqrob_joint.R C:/path/to/h0testr/R > test_msqrob_joint.out 2>&1",
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

log_len <- function() {
  if(!file.exists(log_file)) return(0L)
  return(length(readLines(log_file, warn=FALSE)))
}

log_has <- function(n0, ...) {
  txt <- readLines(log_file, warn=FALSE)
  if(length(txt) <= n0) return(FALSE)
  txt <- paste(txt[-seq_len(n0)], collapse=" ")
  for(pat in c(...)) if(!grepl(pat, txt, fixed=TRUE)) return(FALSE)
  return(TRUE)
}

close_enough <- function(a, b, tol=1e-8) {
  if(length(a) != length(b)) return(FALSE)
  i <- !is.na(a) & !is.na(b)
  if(!any(i)) return(FALSE)
  if(any(is.na(a) != is.na(b))) return(FALSE)
  return(max(abs(a[i] - b[i])) < tol)
}

###############################################################################
## Twenty-four observations, fully crossed: grp has three levels, sex and batch two
##   each, age is continuous. Sixty genes with two peptides each, and a per gene
##   variance rather than a common one, so that limma's prior degrees of freedom are
##   finite: with a common variance msqrob2 returns dfPosterior = Inf, which makes the
##   F and the chi-square coincide and the choice between them untestable.

set.seed(101)

nobs <- 24
ngene <- 60
npep <- 2

samps <- expand.grid(rep=1:2, batch=c("b1", "b2"), sex=c("F", "M"),
  grp=c("a", "b", "c"), stringsAsFactors=FALSE)
samps$obs <- paste0("o", sprintf("%02d", 1:nobs))
samps$age <- round(rnorm(nobs, 50, 10), 1)
samps <- samps[, c("obs", "grp", "sex", "batch", "age")]

## a known effect on the first ten genes: level c is up, level b is not:

gene_eff <- matrix(0, nrow=ngene, ncol=nobs)
gene_eff[1:10, samps$grp %in% "c"] <- 3

sds <- rep(sqrt(1 / rgamma(ngene, shape=3, rate=3)), each=npep)
exprs <- matrix(rnorm(ngene * npep * nobs, 20, rep(sds, times=nobs)),
  nrow=ngene * npep)
exprs <- exprs + gene_eff[rep(1:ngene, each=npep), , drop=FALSE]
gene <- paste0("g", sprintf("%02d", rep(1:ngene, each=npep)))
rownames(exprs) <- paste0(gene, ".p", rep(1:npep, times=ngene))
colnames(exprs) <- samps$obs

feats <- data.frame(pep=rownames(exprs), gene=gene)

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
cfg0$is_log_transformed <- TRUE
cfg0$impute_method <- "none"
cfg0$n_features_min <- 5
cfg0$test_prior_df <- 5
cfg0$log_file <- log_file

out <- initialize(list(expression=exprs, features=feats, samples=samps), cfg0,
  minimal=TRUE)
state0 <- out$state
config0 <- out$config

cfg_for <- function(frm, test_term) {
  cfg <- config0
  cfg$frm <- frm
  cfg$test_term <- test_term
  cfg$reference_levels <-
    config0$reference_levels[names(config0$reference_levels) %in% all.vars(frm)]
  return(cfg)
}

## the truth, for ranking checks:

is_effect <- function(ids) as.integer(sub("^g", "", ids)) <= 10

###############################################################################
section("one column: the package's own path, unchanged")

cfg_1 <- cfg_for(~grp + sex + batch + age, "sex")
res_1 <- suppressMessages(test_msqrob(state0, cfg_1))

report(is.data.frame(res_1$hits) &&
  all(c("logFC", "se", "df", "t", "pval", "adjPval") %in% names(res_1$hits)),
  "one column: hits has the msqrob2::hypothesisTest() columns")

report(nrow(res_1$hits) == ngene && !any(is.na(res_1$hits$pval)),
  "one column: one row per gene, no missing p-values")

report(all(is.finite(res_1$hits$df)) && all(res_1$hits$df > 0),
  "one column: dfPosterior is finite here, so the F and the chi-square differ")

###############################################################################
section("the joint code reproduces msqrob2::hypothesisTest() at one column")

## the same fitted models the branch above used, asked for the same single column
##   through the new code path:

design_1 <- f.design_test_cols(state0, cfg_1)
cols_1 <- colnames(design_1$X)[design_1$cols_test]
report(length(cols_1) == 1 && cols_1 == "sexM", "one column: the column is sexM")

wald_1 <- f.msqrob_wald(res_1$fit, cols_1, cfg_1)
hits_1 <- res_1$hits[match(rownames(wald_1), res_1$hits$gene), ]

report(close_enough(wald_1$pval, hits_1$pval, tol=1e-12),
  "one column: f.msqrob_wald() p-value equals msqrob2::hypothesisTest()")

report(close_enough(wald_1$f_statistic, hits_1$t^2, tol=1e-8),
  "one column: the Wald statistic is the square of the moderated t")

report(all(wald_1$df1 == 1) && close_enough(wald_1$df2, hits_1$df, tol=1e-10),
  "one column: df1 is 1 and df2 is the dfPosterior msqrob2 reports")

report(close_enough(wald_1$adjPval, hits_1$adjPval, tol=1e-12),
  "one column: adjusted p-values agree")

## the decision recorded in the summary: referring the same statistic to a chi-square
##   would disagree with msqrob2's own answer, always downward:

p_chisq <- stats::pchisq(wald_1$f_statistic, df=1, lower.tail=FALSE)
report(all(p_chisq <= wald_1$pval + 1e-12) && any(p_chisq < wald_1$pval - 1e-6),
  "one column: a chi-square would be anticonservative relative to msqrob2")

###############################################################################
section("joint test over a three level factor")

cfg_g <- cfg_for(~grp + sex + batch + age, "grp")
res_g <- try(suppressMessages(test_msqrob(state0, cfg_g)), silent=TRUE)

report(!inherits(res_g, "try-error"),
  "joint: a three level factor no longer refused")

report(is.data.frame(res_g$hits) &&
  all(c("f_statistic", "df1", "df2", "pval", "adjPval") %in% names(res_g$hits)),
  "joint: hits carries the F statistic and both degrees of freedom")

report(all(res_g$hits$df1 == 2), "joint: df1 is the number of columns tested")

report(all(is.finite(res_g$hits$f_statistic)) && all(res_g$hits$f_statistic >= 0),
  "joint: every gene has a finite non-negative statistic")

report(close_enough(res_g$hits$pval,
  stats::pf(res_g$hits$f_statistic, res_g$hits$df1, res_g$hits$df2,
    lower.tail=FALSE)),
  "joint: the p-value is the F tail its own degrees of freedom imply")

report(all(is_effect(res_g$hits$gene[1:10])),
  "joint: the ten genes with the simulated effect are the ten most significant")

###############################################################################
section("joint test agrees with an independent joint test of the same design")

agg <- suppressMessages(combine_features(state0, cfg_g))
tr <- suppressMessages(test(agg$state, agg$config, method="trend"))$standard
cmp <- merge(res_g$hits[, c("gene", "f_statistic", "pval")],
  tr[, c("feature", "pval")], by.x="gene", by.y="feature", suffixes=c(".m", ".t"))

report(stats::cor(cmp$pval.m, cmp$pval.t, method="spearman") > 0.8,
  "joint: ranks agree with limma trend's F test on the same design")

###############################################################################
section("the standardized table")

std <- suppressMessages(test(state0, cfg_g, method="msqrob"))$standard

report(all(c("feature", "expr", "logfc", "stat", "pval", "adj_pval") %in% names(std)),
  "standard: the usual seven columns")

report(close_enough(std$stat,
  res_g$hits$f_statistic[match(std$feature, res_g$hits$gene)]),
  "standard: stat is the F statistic")

## logfc: the total swing of the msqrob2 coefficients over the tested columns,
##   computed here from the fit rather than taken from the package:

design_g <- f.design_test_cols(state0, cfg_g)
cols_g <- colnames(design_g$X)[design_g$cols_test]
models <- SummarizedExperiment::rowData(res_g$fit[["genes"]])$msqrobModels
ids <- as.character(
  SummarizedExperiment::rowData(res_g$fit[["genes"]])[[cfg_g$gene_id_col]])

swing <- rep(as.numeric(NA), length(models))
names(swing) <- ids
for(idx in seq_along(models)) {
  beta <- msqrob2::getCoef(models[[idx]])[cols_g]
  swing[idx] <- diff(range(design_g$X[, design_g$cols_test, drop=FALSE] %*% beta))
}

report(close_enough(std$logfc, swing[std$feature]),
  "standard: logfc is the swing of the msqrob2 coefficients, not of a refit")

report(all(std$logfc >= 0, na.rm=TRUE), "standard: the swing is unsigned")

report(mean(std$logfc[is_effect(std$feature)]) >
  3 * mean(std$logfc[!is_effect(std$feature)]),
  "standard: the swing is much larger for the genes with the simulated effect")

pep_means <- rowMeans(state0$expression, na.rm=TRUE)
gene_means <- tapply(pep_means, as.character(state0$features$gene), mean)
report(close_enough(std$expr, as.numeric(gene_means[std$feature])),
  "standard: expr is the mean over the peptides of each gene")

###############################################################################
section("joint test spanning an interaction")

cfg_i <- cfg_for(~grp * sex + batch, "sex")
res_i <- try(suppressMessages(test_msqrob(state0, cfg_i)), silent=TRUE)

report(!inherits(res_i, "try-error"),
  "interaction: a variable inside an interaction no longer refused")

report(!inherits(res_i, "try-error") && all(res_i$hits$df1 == 3),
  "interaction: marginality puts sexM and both grp:sex columns under test")

report(!inherits(res_i, "try-error") && all(is.finite(res_i$hits$pval)),
  "interaction: every gene has a p-value")

###############################################################################
section("degenerate cases")

## a column the fit does not carry: NA rather than an error, matching what
##   msqrob2::hypothesisTest() does for a fit that failed:

wald_bad <- try(f.msqrob_wald(res_g$fit, c("grpb", "nosuchcolumn"), cfg_g),
  silent=TRUE)
report(!inherits(wald_bad, "try-error") && all(is.na(wald_bad$f_statistic)) &&
  all(is.na(wald_bad$pval)), "degenerate: an unknown column gives NA, not an error")

## no fitted models at all is a package level mistake rather than a data problem,
##   so it is an error:

n0 <- log_len()
obj_bad <- res_g$fit
SummarizedExperiment::rowData(obj_bad[["genes"]])$msqrobModels <- NULL
res_bad <- try(f.msqrob_wald(obj_bad, cols_g, cfg_g), silent=TRUE)
report(inherits(res_bad, "try-error") && log_has(n0, "msqrobModels"),
  "degenerate: a missing model column is an error naming the column")

###############################################################################
section("the caps elsewhere")

## there are none left. f.test_max_cols() and f.design_test_cols_max() are gone, deqms
##   having been the last capped engine: DEqMS::spectraCounteBayes() reports one
##   coefficient's moderated t, but the variance prior it fits is per gene and mentions
##   no coefficient, so test_deqms() forms the joint test from it the way this file's
##   subject forms one from the msqrob2 models. See 1/test_deqms.R, which tests it, and
##   f.deqms_moderated_f():

report(!exists("f.test_max_cols") && !exists("f.design_test_cols_max"),
  "caps: no engine is capped, so both cap helpers are gone")

## this fixture gives every gene the same two peptides, and test_deqms() refuses that
##   whatever is being tested: the count is the covariate its variance prior is fitted
##   against. One peptide dropped from a third of the genes so that the refusal under
##   test here is the one that used to be about the cap:

state_v <- state0
drop <- state0$features$pep %in%
  paste0("g", sprintf("%02d", seq(1, 60, by=3)), ".p2")
state_v$expression <- state0$expression[!drop, , drop=FALSE]
state_v$features <- state0$features[!drop, , drop=FALSE]

res_d <- try(suppressMessages(test_deqms(state_v, cfg_g)), silent=TRUE)
report(!inherits(res_d, "try-error") &&
  all(c("sca.F", "sca.P.Value") %in% names(res_d$hits)) &&
  all(res_d$hits$sca.df.num %in% 2),
  "caps: test_deqms() runs the same joint test and reports a moderated F on 2 df")

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

## a continuous term is one column whatever its spread, so this is the single column path
##   the top of this file pins to msqrob2::hypothesisTest(), not the joint one:

des_ct <- f.design_test_cols(out_ct$state, out_ct$config)
report(length(des_ct$cols_test) %in% 1 &&
  identical(colnames(des_ct$X)[des_ct$cols_test], "age"),
  "a continuous test_term is a single design column, whatever its range")

res_ct <- suppressMessages(test_msqrob(out_ct$state, out_ct$config))
fc_ct <- res_ct$hits$logFC
t_ct <- tru_ct[as.character(res_ct$hits$gene_id)]
q_ct <- stats::p.adjust(res_ct$hits$pval, method="BH")

report(nrow(res_ct$hits) %in% nrow(sim_ct$truth) && !any(is.na(fc_ct)),
  "it is tested rather than refused: one row per gene, no missing logFC")
report(all(sign(fc_ct[t_ct != 0]) == sign(t_ct[t_ct != 0])),
  "the sign of every planted slope is recovered")
report(abs(mean(abs(fc_ct[t_ct != 0])) - slope_ct) < 0.05 * slope_ct,
  "and logFC is the slope per year of age, the planted per SD effect over sd(age)")
report(mean(abs(fc_ct[t_ct == 0])) < 0.1 * slope_ct,
  "while a gene with no age effect is left with a slope near zero")
report(sum(q_ct[t_ct != 0] < 0.05) >= 11 && sum(q_ct[t_ct == 0] < 0.05) <= 2,
  "and the planted genes are the ones that reject at q < 0.05")

## one column carries the test, so f.logfc_effect() passes the signed coefficient through
##   rather than reporting the total swing it reports for the joint cases above:

cfg_std_ct <- out_ct$config
cfg_std_ct$test_method <- "msqrob"
std_ct <- suppressMessages(test(out_ct$state, cfg_std_ct))$standard
report(close_enough(std_ct$logfc, fc_ct[match(std_ct$feature, res_ct$hits$gene_id)]),
  "test() reports that same slope as logfc, not a total swing")

###############################################################################
cat("\n#############################################\n")
cat("\n## passes:", n_pass, "; failures:", n_fail, "; elapsed:",
  round(as.numeric(difftime(Sys.time(), t0, units="secs"))), "s\n")
cat("## log file:", log_file, "\n")

quit(save="no", status=if(n_fail > 0) 1 else 0)
