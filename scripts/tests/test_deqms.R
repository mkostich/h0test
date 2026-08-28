## Tests the moderated test test_deqms() runs when config$test_term resolves to more
##   than one design matrix column, and the statistic it reports in either case.
##   DEqMS::spectraCounteBayes() moderates the t-statistic of one coefficient and the
##   package has no F-analogue, so those tests used to be refused. Its moderation has no
##   such limit: that function fits a variance prior against the number of features
##   behind each gene and returns a per-gene posterior variance and a prior degrees of
##   freedom, neither of which mentions a coefficient, so f.deqms_moderated_f() forms the
##   joint test from them.
##   The central assertion is that the new code reproduces DEqMS's own moderated t and
##   p-value exactly when handed one column, which is what justifies referring the joint
##   statistic to an F on the same posterior degrees of freedom. The joint case is
##   checked against an independent nested-model computation that does not go through
##   limma's covariance matrix at all.
##   Also covers the standardized table, which used to receive limma's t and P.Value
##   from f.format_limma() because DEqMS's table carries limma's columns too, so that
##   test_method "deqms" reported limma::eBayes(trend=FALSE) in every column anything
##   downstream reads.

usage <- function(msg=NULL) {
  if(!is.null(msg)) cat("ERROR:", msg, "\n\n", file=stderr())
  cat(
    "Test the moderated test h0testr::test_deqms() runs when config$test_term spans",
    "one or more design matrix columns. Checks that the single-column case reproduces",
    "DEqMS::spectraCounteBayes()'s own sca.t and sca.p exactly, that the joint case",
    "agrees with an independent nested-model F with the same posterior variance, that",
    "the standardized table reports DEqMS's statistic rather than limma's, that logfc",
    "holds the total swing for a joint test, that a config$contrast reaches the engine",
    "as a single coefficient, that a joint test refuses a matrix with missing values,",
    "that uniform feature counts are still refused, that a continuous config$test_term",
    "is recovered as a slope per unit of the covariate, that the variance prior is",
    "fitted from the genes that have residual degrees of freedom and reported as NA for",
    "the rest rather than misaligned across all of them, and the degenerate cases of",
    "f.deqms_moderated_f() and f.format_deqms().",
    "",
    "Usage: Rscript test_deqms.R <r_dir>",
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
    "  3 DEqMS not installed.",
    "",
    "Examples:",
    "  Rscript test_deqms.R C:/path/to/h0test/h0testr/R",
    "  Rscript test_deqms.R ../../h0test/h0testr/R",
    "  Rscript test_deqms.R C:/path/to/h0testr/R > test_deqms.out 2>&1",
    sep="\n", file=stderr()
  )
  quit(status=2)
}

args <- commandArgs(trailingOnly=TRUE)
if(length(args) != 1) usage("wrong number of arguments")

r_dir <- args[1]
if(!dir.exists(r_dir)) usage(paste("r_dir not a directory:", r_dir))

if(!requireNamespace("DEqMS", quietly=TRUE)) {
  cat("ERROR: DEqMS not installed\n", file=stderr())
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

errs_with <- function(expr, ...) {
  n0 <- log_len()
  res <- try(suppressMessages(expr), silent=TRUE)
  return(inherits(res, "try-error") && log_has(n0, ...))
}

close_enough <- function(a, b, tol=1e-8) {
  if(length(a) != length(b)) return(FALSE)
  if(any(is.na(a) != is.na(b))) return(FALSE)
  i <- !is.na(a) & !is.na(b)
  if(!any(i)) return(FALSE)
  return(max(abs(a[i] - b[i])) < tol)
}

###############################################################################
## Twenty-four observations, fully crossed: grp has three levels, sex and batch two
##   each, age is continuous, so that even ~grp*sex*batch is of full rank. Sixty genes
##   whose peptide counts run from one to six: DEqMS fits its variance prior as a loess
##   of the log residual variance on the log count, and needs those counts spread rather
##   than merely varied. One distinct count leaves that loess with nothing to fit, and so
##   does a count all but a few genes share; both fail with an error, which the last
##   section covers along with the genes such a prior cannot be fitted for.

set.seed(101)

nobs <- 24
ngene <- 60
npep_max <- 6

samps <- expand.grid(rep=1:2, batch=c("b1", "b2"), sex=c("F", "M"),
  grp=c("a", "b", "c"), stringsAsFactors=FALSE)
samps <- data.frame(
  obs=paste0("o", sprintf("%02d", 1:nobs)),
  grp=samps$grp,
  sex=samps$sex,
  batch=samps$batch,
  age=round(stats::rnorm(nobs, 50, 10), 1),
  stringsAsFactors=FALSE
)

## the first ten genes are up in group c and the second ten up in males of batch b2,
##   so that a three-level factor and an interaction each have something to find:

gene_eff <- matrix(0, nrow=ngene, ncol=nobs)
gene_eff[1:10, samps$grp %in% "c"] <- 3
gene_eff[11:20, samps$sex %in% "M" & samps$batch %in% "b2"] <- 3

npep <- ((seq_len(ngene) - 1) %% npep_max) + 1        ## 1..6 peptides per gene
gene <- rep(paste0("g", sprintf("%02d", seq_len(ngene))), times=npep)

exprs <- matrix(stats::rnorm(length(gene) * nobs, 20, 1.5), nrow=length(gene))
exprs <- exprs + gene_eff[match(gene, paste0("g", sprintf("%02d", seq_len(ngene)))), ,
  drop=FALSE]
rownames(exprs) <- paste0(gene, ".p", unlist(lapply(npep, seq_len)))
colnames(exprs) <- samps$obs
feats <- data.frame(pep=rownames(exprs), gene=gene, stringsAsFactors=FALSE)

cfg0 <- new_config()
cfg0$obs_id_col <- cfg0$obs_col <- "obs"
cfg0$sample_id_col <- "obs"
cfg0$feat_id_col <- cfg0$feat_col <- "pep"
cfg0$gene_id_col <- "gene"
cfg0$reference_levels <- c(grp="a", sex="F", batch="b1")
cfg0$estimability <- "test"
cfg0$df_resid_min <- 2
cfg0$save_state <- FALSE
cfg0$permute_var <- ""
cfg0$is_log_transformed <- TRUE      ## simulated on a log scale; normalize() not called
cfg0$impute_method <- "none"
cfg0$n_features_min <- 5             ## default of 1000 would filter out everything
cfg0$log_file <- log_file
cfg0$frm <- ~grp + sex + batch + age
cfg0$test_term <- "sex"
cfg0$test_method <- "deqms"

out <- init_state(list(expression=exprs, features=feats, samples=samps), cfg0,
  minimal=TRUE)
state0 <- out$state
config0 <- out$config

up_grp <- paste0("g", sprintf("%02d", 1:10))
up_int <- paste0("g", sprintf("%02d", 11:20))

cfg_trm <- function(frm, term) {
  cfg <- config0
  cfg$frm <- frm
  cfg$test_term <- term
  cfg$contrast <- ""
  vars <- sort(unique(f.parse_frm(frm, cfg)$vars))
  cfg$reference_levels <- config0$reference_levels[
    names(config0$reference_levels) %in% vars]
  return(cfg)
}

cfg_con <- function(frm, con) {
  cfg <- cfg_trm(frm, "sex")
  cfg$test_term <- ""
  cfg$contrast <- con
  return(cfg)
}

## the aggregated matrix test_deqms() tests, rebuilt here so that the independent
##   computations below start from the same values it does:

agg_of <- function(cfg) {
  cf <- cfg
  cf$save_state <- FALSE
  return(suppressMessages(combine_features(state0, cf, method="medianPolish",
    rescale=FALSE))$state)
}

###############################################################################
section("one coefficient: the numbers are DEqMS's own")

cfg_1 <- cfg_trm(~grp + sex + batch + age, "sex")
res_1 <- suppressMessages(test_deqms(state0, cfg_1))
des_1 <- f.design_test_cols(agg_of(cfg_1), cfg_1)

report(length(des_1$cols_test) %in% 1 &&
  identical(colnames(des_1$X)[des_1$cols_test], "sexM"),
  "testing 'sex' in ~grp + sex + batch + age is a 1 df test of sexM")

report(all(res_1$hits$sca.df.num %in% 1) && all(res_1$hits$sca.df.den > 0),
  "the reported numerator df is 1 and the denominator df is positive")

## DEqMS::spectraCounteBayes() is called without a coef_col, so the sca.t and sca.p it
##   leaves on the fit cover every coefficient: fit$sca.t[, 'sexM'] is that function's
##   own answer for this test, computed by DEqMS and not by h0testr:

idx <- match("sexM", colnames(res_1$fit$coefficients))
own_t <- res_1$fit$sca.t[res_1$hits$gene, idx]
own_p <- res_1$fit$sca.p[res_1$hits$gene, idx]

report(close_enough(res_1$hits$sca.t, own_t, tol=1e-10),
  "sca.t is exactly what DEqMS::spectraCounteBayes() computed for that coefficient")
report(close_enough(res_1$hits$sca.P.Value, own_p, tol=1e-12),
  "sca.P.Value is exactly what DEqMS::spectraCounteBayes() computed")

report(close_enough(res_1$hits$sca.F, res_1$hits$sca.t^2, tol=1e-10),
  "at 1 df the moderated F is the square of DEqMS's moderated t")
report(close_enough(res_1$hits$sca.P.Value,
  stats::pf(res_1$hits$sca.F, 1, res_1$hits$sca.df.den, lower.tail=FALSE), tol=1e-12),
  "and the F on 1 df gives back the same p-value")

j <- match(res_1$hits$gene, rownames(res_1$fit$coefficients))
report(close_enough(res_1$hits$sca.df.den,
  res_1$fit$sca.dfprior + res_1$fit$df.residual[j], tol=1e-10),
  "the denominator df is DEqMS's prior df plus the residual df of the fit")

report(all(is.finite(res_1$fit$sca.postvar)) &&
  length(unique(res_1$fit$count)) >= 5,
  "the variance prior was fitted for every gene, the feature counts having spread")

###############################################################################
section("the standardized table reports DEqMS, not limma")

std_1 <- suppressMessages(test_h0(state0, cfg_1))$standard
i <- match(std_1$feature, res_1$hits$gene)

report(close_enough(std_1$pval, res_1$hits$sca.P.Value[i], tol=1e-12),
  "standard$pval is DEqMS's sca.P.Value")
report(close_enough(std_1$adj_pval, res_1$hits$sca.adj.pval[i], tol=1e-12),
  "standard$adj_pval is DEqMS's sca.adj.pval")
report(close_enough(std_1$stat, res_1$hits$sca.t[i], tol=1e-10),
  "standard$stat is DEqMS's sca.t")

## and they are not limma's, which is the bug this replaced: f.format_limma() matched on
##   the limma columns of DEqMS's table, all of which are present, so the standardized
##   table used to carry limma's t and P.Value and DEqMS's prior reached nothing:

report(!close_enough(std_1$pval, res_1$hits$P.Value[i], tol=1e-6),
  "standard$pval is not limma's P.Value")
report(!close_enough(std_1$stat, res_1$hits$t[i], tol=1e-6),
  "standard$stat is not limma's t")

report(close_enough(std_1$logfc, res_1$hits$logFC[i], tol=1e-10),
  "standard$logfc is the coefficient under test, which limma and DEqMS share")
report(all(is.na(std_1$lod)),
  "standard$lod is empty rather than limma's B, which comes from limma's own prior")

###############################################################################
section("several coefficients: a three-level factor")

cfg_g <- cfg_trm(~grp + sex + batch + age, "grp")
res_g <- suppressMessages(test_deqms(state0, cfg_g))
agg_g <- agg_of(cfg_g)
des_g <- f.design_test_cols(agg_g, cfg_g)

report(length(des_g$cols_test) %in% 2 &&
  identical(sort(colnames(des_g$X)[des_g$cols_test]), c("grpb", "grpc")),
  "testing 'grp' is a 2 df joint test of grpb and grpc")

report(all(res_g$hits$sca.df.num %in% 2),
  "the reported numerator df is 2")
report(all(is.na(res_g$hits$sca.t)) && !any(is.na(res_g$hits$sca.F)),
  "sca.t is NA for a joint test and sca.F is not")
report(!("logFC" %in% names(res_g$hits)) &&
  all(c("grpb", "grpc") %in% names(res_g$hits)),
  "limma::topTable() reports one coefficient column per tested column, not a logFC")

## the joint statistic, independently: the extra sum of squares of the nested model
##   comparison, per gene, over the numerator df, divided by DEqMS's posterior variance.
##   This goes nowhere near fit$cov.coefficients, which is what f.deqms_moderated_f()
##   uses, so agreement is a check of the algebra and not of the same code twice:

## named xf/xr rather than X/X_red because base::apply() has an argument of its own
##   called X, and a named X= here would be taken as that one:

ess <- function(y, xf, xr) {
  ok <- !is.na(y)
  r_full <- stats::lm.fit(xf[ok, , drop=FALSE], y[ok])$residuals
  r_red <- stats::lm.fit(xr[ok, , drop=FALSE], y[ok])$residuals
  return(sum(r_red^2) - sum(r_full^2))
}

q_g <- apply(agg_g$expression, 1, ess, xf=des_g$X, xr=des_g$X_red)
f_ref <- (q_g / 2) / res_g$fit$sca.postvar[match(names(q_g),
  rownames(res_g$fit$coefficients))]
f_ref <- f_ref[res_g$hits$gene]

report(close_enough(res_g$hits$sca.F, f_ref, tol=1e-8),
  "the joint moderated F equals the nested extra sum of squares over the posterior var")

p_ref <- stats::pf(f_ref, 2, res_g$hits$sca.df.den, lower.tail=FALSE)
report(close_enough(res_g$hits$sca.P.Value, p_ref, tol=1e-12),
  "and its p-value is that F on 2 and the posterior degrees of freedom")

## the same F with limma's own posterior variance is limma's F, which is what
##   test_trend(trend=FALSE) would report: same numerator, different prior:

report(!close_enough(res_g$hits$sca.P.Value, res_g$hits$P.Value, tol=1e-6) &&
  stats::cor(res_g$hits$sca.P.Value, res_g$hits$P.Value) > 0.9,
  "DEqMS's prior gives different p-values from limma's, and correlated ones")

report(sum(utils::head(res_g$hits$gene, 10) %in% up_grp) >= 8,
  "the joint test ranks the genes with the simulated group effect first")

###############################################################################
section("several coefficients: a variable inside an interaction")

cfg_i <- cfg_trm(~sex * batch, "sex")
res_i <- suppressMessages(test_deqms(state0, cfg_i))
des_i <- f.design_test_cols(agg_of(cfg_i), cfg_i)

report(length(des_i$cols_test) %in% 2 &&
  identical(sort(colnames(des_i$X)[des_i$cols_test]), c("sexM", "sexM:batchb2")),
  "testing 'sex' in ~sex * batch is a 2 df joint test, by marginality")
report(all(res_i$hits$sca.df.num %in% 2) && !any(is.na(res_i$hits$sca.F)),
  "test_deqms() runs it rather than testing sexM alone")
report(sum(utils::head(res_i$hits$gene, 10) %in% up_int) >= 8,
  "and ranks the genes with the simulated interaction effect first")

###############################################################################
section("the effect size of a joint test")

std_g <- suppressMessages(test_h0(state0, cfg_g))$standard

## f.logfc_effect()'s total swing: the range, over the observations, of the fitted
##   contribution of the tested columns. For grp that is the largest difference between
##   any two of its three levels, including the b-to-c pair no coefficient names:

b_g <- res_g$fit$coefficients[, des_g$cols_test, drop=FALSE]
fitted_g <- b_g %*% t(des_g$X[, des_g$cols_test, drop=FALSE])
swing <- apply(fitted_g, 1, function(v) diff(range(v)))
swing <- swing[std_g$feature]

report(close_enough(std_g$logfc, swing, tol=1e-8),
  "standard$logfc is the total swing over the tested coefficients")
report(all(std_g$logfc >= 0, na.rm=TRUE),
  "the total swing is unsigned")
report(stats::median(std_g$logfc[std_g$feature %in% up_grp]) > 2 &&
  stats::median(std_g$logfc[!(std_g$feature %in% up_grp)]) < 1.5,
  "and it is larger for the genes with the simulated 3-unit group effect")

report(close_enough(std_g$stat, res_g$hits$sca.F[match(std_g$feature,
  res_g$hits$gene)], tol=1e-10),
  "standard$stat is the moderated F for a joint test")

###############################################################################
section("config$contrast reaches the engine as one coefficient")

cfg_c <- cfg_con(~grp + sex + batch + age, "grpc - grpb")
res_c <- suppressMessages(test_deqms(state0, cfg_c))

report(all(res_c$hits$sca.df.num %in% 1) && !any(is.na(res_c$hits$sca.t)),
  "a contrast within a three-level factor is 1 df, so DEqMS's own t is reported")
report("logFC" %in% names(res_c$hits),
  "and limma::contrasts.fit() leaves a single coefficient with a logFC")

std_c <- suppressMessages(test_h0(state0, cfg_c))$standard
report(close_enough(std_c$stat, res_c$hits$sca.t[match(std_c$feature,
  res_c$hits$gene)], tol=1e-10),
  "the standardized table reports that t")

## and the contrast is a different hypothesis from the joint test of the term:

report(!close_enough(res_c$hits$sca.P.Value[order(res_c$hits$gene)],
  res_g$hits$sca.P.Value[order(res_g$hits$gene)], tol=1e-6),
  "the contrast and the joint test of the same factor are different hypotheses")

###############################################################################
section("a joint test needs a complete matrix")

## limma fits each gene on the observations that gene has, so fit$stdev.unscaled is per
##   gene while fit$cov.coefficients, which the joint test needs, comes from the
##   complete design. Refused rather than reported wrong:

state_na <- state0
i_na <- state_na$features[[config0$gene_id_col]] %in% "g01"
state_na$expression[i_na, 1] <- NA

report(errs_with(test_deqms(state_na, cfg_g),
  "no per-gene covariance for the coefficients under test"),
  "a joint test refuses an aggregated matrix with a missing value, naming the reason")

res_na1 <- try(suppressMessages(test_deqms(state_na, cfg_1)), silent=TRUE)
report(!inherits(res_na1, "try-error") && all(res_na1$hits$sca.df.num %in% 1),
  "a single coefficient still runs on the same matrix, reading stdev.unscaled per gene")

res_nac <- try(suppressMessages(test_deqms(state_na, cfg_c)), silent=TRUE)
report(!inherits(res_nac, "try-error") && all(res_nac$hits$sca.df.num %in% 1),
  "and so does a contrast, which is one coefficient after contrasts.fit()")

###############################################################################
section("what is still refused, and the degenerate cases")

## the variance prior is fitted against the number of features per gene, so it needs
##   that number to vary. Unchanged by the joint test, and checked here as a regression:

state_u <- state0
first_pep <- !duplicated(state0$features[[config0$gene_id_col]])
state_u$expression <- state0$expression[first_pep, , drop=FALSE]
state_u$features <- state0$features[first_pep, , drop=FALSE]

report(errs_with(test_deqms(state_u, cfg_1), "every gene has the same number of"),
  "uniform feature counts are still refused, whatever is being tested")
report(errs_with(test_deqms(state_u, cfg_g), "every gene has the same number of"),
  "including for a joint test")

## f.deqms_moderated_f() on something that is not a DEqMS fit:

fit_bad <- res_g$fit
fit_bad$sca.postvar <- NULL
report(errs_with(f.deqms_moderated_f(fit_bad, des_g$cols_test, cfg_g),
  "did not come from DEqMS::spectraCounteBayes()"),
  "f.deqms_moderated_f() refuses a fit without the moderated variance")

fit_bad <- res_g$fit
fit_bad$cov.coefficients <- fit_bad$cov.coefficients[1, 1, drop=FALSE]
report(errs_with(f.deqms_moderated_f(fit_bad, des_g$cols_test, cfg_g),
  "does not carry an unscaled covariance"),
  "f.deqms_moderated_f() refuses a fit whose covariance lacks a tested coefficient")

## a gene whose posterior variance could not be fitted comes back NA, with a count:

n0 <- log_len()
fit_nan <- res_g$fit
fit_nan$sca.postvar[1:3] <- NaN
mod_nan <- suppressMessages(f.deqms_moderated_f(fit_nan, des_g$cols_test, cfg_g))
report(sum(is.na(mod_nan$p.value)) %in% 3 && log_has(n0, "have no"),
  "a gene with no posterior variance gives NA and is counted in a warning")

## f.format_deqms() on a table that cannot have come from one hypothesis:

tbl_bad <- res_g$hits
tbl_bad$sca.df.num[1] <- 3
report(errs_with(f.format_deqms(tbl_bad, cfg_g), "different numerator"),
  "f.format_deqms() refuses a table reporting more than one numerator df")

tbl_bad <- res_g$hits
tbl_bad$sca.P.Value <- NULL
report(errs_with(f.format_deqms(tbl_bad, cfg_g), "sca.P.Value"),
  "f.format_deqms() names the column it is missing")

## and there is no cap helper left to consult:

report(!exists("f.test_max_cols") && !exists("f.design_test_cols_max"),
  "f.test_max_cols() and f.design_test_cols_max() are gone")

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

res_ct <- suppressMessages(test_deqms(out_ct$state, out_ct$config))
fc_ct <- res_ct$hits$logFC
t_ct <- tru_ct[as.character(res_ct$hits$gene)]
q_mod_ct <- stats::p.adjust(res_ct$hits$sca.P.Value, method="BH")
q_raw_ct <- stats::p.adjust(res_ct$hits$P.Value, method="BH")

report(nrow(res_ct$hits) %in% nrow(sim_ct$truth) && !any(is.na(fc_ct)),
  "a continuous test_term is tested: one row per gene, no missing logFC")
report(all(res_ct$hits$sca.df.num %in% 1),
  "a continuous term is one design column, so the moderated test is on 1 df")
report(all(sign(fc_ct[t_ct != 0]) == sign(t_ct[t_ct != 0])),
  "the sign of every planted slope is recovered")
report(abs(mean(abs(fc_ct[t_ct != 0])) - slope_ct) < 0.05 * slope_ct,
  "and logFC is the slope per year of age, the planted per SD effect over sd(age)")
report(mean(abs(fc_ct[t_ct == 0])) < 0.1 * slope_ct,
  "while a gene with no age effect is left with a slope near zero")
report(sum(q_mod_ct[t_ct != 0] < 0.05) >= 11 && sum(q_mod_ct[t_ct == 0] < 0.05) <= 2,
  "and the planted genes are the ones that reject at q < 0.05")
report(sum(q_raw_ct[t_ct != 0] < 0.05) >= 11,
  "which limma's own unmoderated p-value agrees about")

## the fixture gives its genes 6 peptides and drops them heavily, so the feature counts
##   spread out. That is deliberate and worth asserting rather than assuming: DEqMS fits
##   its variance prior as a loess of the residual variance on the feature count, and
##   fails outright when the counts are lopsided enough to leave that loess singular.
##   Nothing to do with a continuous term, which is why the check is on the fixture
##   rather than on the term; the last section covers the failure itself:

report(length(unique(res_ct$hits$count)) > 2 && !any(is.na(res_ct$hits$sca.t)),
  "the feature counts vary enough for the variance prior to be fitted")

## and the standardized table carries the slope itself: one design column carries the test,
##   so f.logfc_effect() passes the signed coefficient through rather than a total swing:

cfg_std_ct <- out_ct$config
cfg_std_ct$test_method <- "deqms"
std_ct <- suppressMessages(test_h0(out_ct$state, cfg_std_ct))$standard
report(close_enough(std_ct$logfc, fc_ct[match(std_ct$feature, res_ct$hits$gene)]),
  "test_h0() reports that same slope as logfc")

###############################################################################
section("the genes DEqMS's variance prior is fitted from")

## DEqMS::spectraCounteBayes() fits that prior with
##   loess(log(fit$sigma^2) ~ log2(fit$count)), and stats::loess() defaults to
##   na.action=na.omit, so a gene whose residual variance is not finite is dropped from
##   the fit and stats::fitted() comes back shorter than the fit. DEqMS then recycles
##   those predictions against fit$df.residual, which hands every gene at or after the
##   dropped row another gene's prior variance. Silently: the recycling restores the
##   length, so f.deqms_moderated_f()'s own length check cannot see it, and the prior a
##   gene receives depends on nothing but where the untestable genes sit in
##   state$expression. f.deqms_prior() fits the prior on the genes that can carry one
##   instead. A gene with as many observations as the design has columns is the reachable
##   case, unimputed missing values being what leaves one behind:

cfg_z <- cfg_trm(~sex, "sex")            ## two design columns: intercept and sexM
i_f <- which(state0$samples$sex %in% "F")[1]
i_m <- which(state0$samples$sex %in% "M")[1]

## two observations, one per sex: full rank, and no residual degree of freedom left:

cripple <- function(state, genes) {
  i <- state$features$gene %in% genes
  state$expression[i, -c(i_f, i_m)] <- NA
  return(state)
}

drop_genes <- function(state, genes) {
  i <- !(state$features$gene %in% genes)
  state$expression <- state$expression[i, , drop=FALSE]
  state$features <- state$features[i, , drop=FALSE]
  return(state)
}

dead <- c("g05", "g11")                  ## five peptides each, so the counts stay spread
state_z <- cripple(state0, dead)

n0 <- log_len()
res_z <- suppressMessages(test_deqms(state_z, cfg_z))
na_z <- is.na(res_z$hits$sca.t)

report(nrow(res_z$hits) %in% ngene &&
  identical(sort(res_z$hits$gene[na_z]), sort(dead)),
  "a gene with no residual degrees of freedom keeps its row and loses its statistic")
report(all(is.finite(res_z$hits$sca.t[!na_z])),
  "and every other gene is moderated")
## limma still reports for such a gene: limma::eBayes() shrinks it toward its own prior,
##   which needs no residual variance of the gene's own, and tests it on the prior's
##   degrees of freedom alone. DEqMS's prior is a loess against the feature counts and has
##   nothing to say at a gene it was not fitted from, so only the sca. columns go missing:

report(all(is.finite(res_z$hits$P.Value[na_z])),
  "limma's own p-value survives there, coming from a prior of its own")
report(log_has(n0, "left out of DEqMS", "no residual degrees of freedom:"),
  "the log says how many genes were left out of the prior, and why")
report(sum(!is.na(res_z$hits$sca.adj.pval)) %in% (ngene - length(dead)),
  "and they are left out of the multiplicity correction rather than counted in it")

## the prior the rest receive is DEqMS's own, fitted from them: the same run with those
##   genes' features removed from the matrix altogether reports the same moderated
##   p-value for every gene the two share, which is what a correctly aligned prior means.
##   limma's own columns do move, its prior being pooled over the genes present:

res_y <- suppressMessages(test_deqms(drop_genes(state0, dead), cfg_z))
i_zy <- match(res_y$hits$gene, res_z$hits$gene)

report(!any(is.na(i_zy)) &&
  close_enough(res_y$hits$sca.P.Value, res_z$hits$sca.P.Value[i_zy]),
  "the surviving genes get the prior a fit of those genes alone would give them")

## and it no longer matters where in the matrix the untestable genes sit, which is the
##   symptom the recycling had: measured outside this suite on 200 genes with one
##   ineligible among them, a gene's sca.postvar was 1.122234 with the ineligible gene in
##   row 1 and 1.126142 with it in row 100:

set.seed(202)
perm <- sample(nrow(state_z$expression))
state_p <- state_z
state_p$expression <- state_p$expression[perm, , drop=FALSE]
state_p$features <- state_p$features[perm, , drop=FALSE]

res_p <- suppressMessages(test_deqms(state_p, cfg_z))
i_zp <- match(res_z$hits$gene, res_p$hits$gene)

report(!any(is.na(i_zp)) &&
  close_enough(res_z$hits$sca.P.Value, res_p$hits$sca.P.Value[i_zp]),
  "and permuting the rows of state$expression leaves every p-value where it was")

## a complete matrix reaches none of this, and says nothing about exclusions:

n0 <- log_len()
res_ok <- suppressMessages(test_deqms(state0, cfg_z))

report(!any(is.na(res_ok$hits$sca.t)) && !log_has(n0, "left out of DEqMS"),
  "a complete matrix leaves every gene in the prior, and the log quiet about it")

## what is refused. A prior needs at least two genes to be fitted from at all, which is
##   checked before limma::eBayes() gets the chance to refuse the same fit with a message
##   that says nothing about imputing or filtering:

report(errs_with(test_deqms(cripple(state0, unique(state0$features$gene)), cfg_z),
  "have a", "residual degree of freedom", "config$df_resid_min"),
  "test_deqms() refuses a fit in which no gene has residual degrees of freedom")

## and the spread of feature counts has to survive the exclusion, the up-front guard in
##   test_deqms() having counted every gene rather than the eligible ones:

n_by_gene <- table(state0$features$gene)
six <- names(n_by_gene)[n_by_gene %in% 6]

report(length(six) > 1 && errs_with(
  test_deqms(cripple(state0, setdiff(names(n_by_gene), six)), cfg_z),
  "all have the same number of features"),
  "and refuses a fit whose eligible genes all have the same feature count")

## DEqMS needs those counts spread rather than merely varied, and answers a shortage two
##   ways. Its loess comes back NaN for the genes at a crowded count, which are then
##   reported without a moderated statistic and counted in the log, and when it comes back
##   NaN for every gene the function fails outright with an error about a missing value in
##   an if() condition, which is refused. Two of the sixty genes here have one peptide and
##   the rest are given two, which gives the first:

set.seed(303)
npep_lop <- c(1, 1, rep(2, ngene - 2))
gene_lop <- rep(paste0("h", sprintf("%02d", seq_len(ngene))), times=npep_lop)
exprs_lop <- matrix(stats::rnorm(length(gene_lop) * nobs, 20, 1.5),
  nrow=length(gene_lop))
rownames(exprs_lop) <- paste0(gene_lop, ".p", unlist(lapply(npep_lop, seq_len)))
colnames(exprs_lop) <- state0$samples[[config0$obs_col]]

state_l <- list(
  expression=exprs_lop,
  features=data.frame(pep=rownames(exprs_lop), gene=gene_lop, stringsAsFactors=FALSE),
  samples=state0$samples
)

n0 <- log_len()
res_l <- suppressMessages(test_deqms(state_l, cfg_z))
na_l <- is.na(res_l$hits$sca.t)

report(sum(na_l) > ngene / 2 && all(is.finite(res_l$hits$P.Value)),
  "lopsided counts cost most genes their moderated statistic, limma's surviving")
report(log_has(n0, "non-finite posterior variance", "features per gene (count x genes):",
  "1x2 2x58"),
  "and the log counts them and names the feature counts, DEqMS's answer naming neither")

###############################################################################
cat("\n#############################################\n")
cat("\n## passes:", n_pass, "; failures:", n_fail, "; elapsed:",
  round(as.numeric(difftime(Sys.time(), t0, units="secs"))), "s\n")
cat("## log file:", log_file, "\n")

quit(save="no", status=if(n_fail > 0) 1 else 0)
