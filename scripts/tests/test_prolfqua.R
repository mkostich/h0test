## Tests for covariate handling in h0testr::test_prolfqua(): which covariates
##   are accepted, and whether the level ordering set by init_state() (declared
##   reference level first) survives into the fitted model.

usage <- function(msg=NULL) {
  if(!is.null(msg)) cat("ERROR:", msg, "\n\n", file=stderr())
  cat(
    "Test covariate handling in h0testr::test_prolfqua(): that factor and",
    "continuous covariates are both accepted and fit as their own kind, that the",
    "reference level declared in config$reference_levels is preserved, that a",
    "factor by continuous interaction can be tested while a main effect an",
    "interaction contains cannot, that the reported F has its error variance",
    "moderated across features by the same estimator limma uses, toward a flat",
    "prior or one fitted against mean intensity, and how missingness that leaves",
    "a covariate aliased for one feature is handled.",
    "",
    "Usage: Rscript test_prolfqua.R <r_dir>",
    "",
    "Required positional arguments:",
    "  <r_dir>  Path to the h0testr package R/ source directory; all .R files",
    "             there are sourced (the installed package is not used).",
    "",
    "Requires the prolfqua and limma packages to be installed.",
    "",
    "Output: one PASS/FAIL line per assertion to stdout, then a count of",
    "  passes and failures. Messages from expected errors are written to a",
    "  temporary log file, whose path is reported at the end.",
    "",
    "Exit codes: 0 all assertions passed; 1 one or more failed; 2 usage error;",
    "  3 prolfqua not installed; 4 limma not installed.",
    "",
    "Examples:",
    "  Rscript test_prolfqua.R C:/path/to/h0test/h0testr/R",
    "  Rscript test_prolfqua.R ../../h0test/h0testr/R",
    "  Rscript test_prolfqua.R C:/path/to/h0testr/R > test_prolfqua.out 2>&1",
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

if(!requireNamespace("limma", quietly=TRUE)) {
  cat("ERROR: limma not installed\n", file=stderr())
  quit(status=4)
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

## returns TRUE iff expr threw an error:
threw <- function(expr) inherits(try(expr, silent=TRUE), "try-error")

## TRUE iff pat appears in the log file at or after line `since`; f.err() stops
##   with "Stopping" and writes the detail to the log:
mark <- function() length(readLines(log_file))
logged <- function(pat, since=0) {
  txt <- readLines(log_file)
  if(since >= length(txt)) return(FALSE)
  any(grepl(pat, txt[(since + 1):length(txt)], fixed=TRUE))
}

## first per-feature lm fit from a test_prolfqua() result:
lm1 <- function(res) res$fit$modelDF$linear_model[[1]]

## the covariates reach the fit as numeric design matrix columns, so the fit has no
##   xlevels; the level ordering shows up instead as which columns model.matrix()
##   emitted, the reference level of each factor being the one with no column:
cols <- function(res) colnames(res$design$X)

###############################################################################
## shared data; prolfqua needs enough features to fit, so use sim_design():

set.seed(101)
nsamps <- 6
sim <- sim_design(
  sim_samples(factors=list(grp=c("ctl", "trt")), n_per_cell=nsamps),
  frm=~grp, test_term="grp", n_genes=20, n_genes_signif=5,
  effects=1, peps_per_gene=3, reps_per_sample=1,
  p_drop=0.1, mnar_c0=-Inf, mnar_c1=0, mcar_p=0
)
exprs <- sim$state$expression
gene <- sim$state$features$gene_id
feats <- data.frame(pep=rownames(exprs), gene=gene)
samps0 <- data.frame(
  obs=colnames(exprs),
  grp=c(rep("ctl", nsamps), rep("trt", nsamps))
)

mk_state <- function(samps) list(expression=exprs, features=feats, samples=samps)

## declare the non-alphabetical level as the reference, so that any silent
##   re-sorting downstream is visible:
cfg0 <- list(
  obs_id_col="obs", sample_id_col="obs", feat_id_col="pep", gene_id_col="gene",
  frm=~grp, test_term="grp", reference_levels=c(grp="trt"), log_file=log_file
)

###############################################################################
section("declared reference level survives into the fit")

out <- init_state(mk_state(samps0), cfg0, minimal=TRUE)
report(identical(levels(out$state$samples$grp), c("trt", "ctl")),
  "init_state() puts declared reference level first")

res <- test_prolfqua(out$state, out$config, is_log_transformed=FALSE)
report(identical(cols(res), c("(Intercept)", "grpctl")),
  "design has a column for the non-reference level only")
report("grpctl" %in% names(stats::coef(lm1(res))),
  "coefficients are relative to the declared reference level")
report(length(lm1(res)$xlevels) %in% 0,
  "the fit itself holds no factors, only numeric design columns")
report(nrow(res$hits) > 0 && all(res$hits$factor %in% "grp"),
  "anova hits returned for test_term")

###############################################################################
section("covariates that are categorical without being declared")

samps <- samps0
samps$flag <- rep(c(TRUE, FALSE), nsamps)
cfg <- cfg0
cfg$frm <- ~grp + flag

out <- init_state(mk_state(samps), cfg, minimal=TRUE)
report(identical(unname(out$config$covariate_types[["flag"]]), "factor"),
  "undeclared logical covariate classified as factor")

res <- try(test_prolfqua(out$state, out$config, is_log_transformed=FALSE),
  silent=TRUE)
report(!inherits(res, "try-error"),
  "test_prolfqua() accepts undeclared logical covariate")
if(!inherits(res, "try-error")) {
  report("flagTRUE" %in% cols(res) && !("flagFALSE" %in% cols(res)),
    "logical covariate coded with FALSE as the reference level")
}

samps <- samps0
samps$sex <- factor(rep(c("M", "F"), nsamps), levels=c("M", "F"))
cfg <- cfg0
cfg$frm <- ~grp + sex

out <- init_state(mk_state(samps), cfg, minimal=TRUE)
res <- try(test_prolfqua(out$state, out$config, is_log_transformed=FALSE),
  silent=TRUE)
report(!inherits(res, "try-error"),
  "test_prolfqua() accepts covariate already stored as a factor")
if(!inherits(res, "try-error")) {
  report("sexF" %in% cols(res) && !("sexM" %in% cols(res)),
    "pre-existing factor keeps its own level ordering, M as reference")
}

###############################################################################
section("continuous covariates")

samps <- samps0
samps$age <- seq(2, 24, length.out=nrow(samps0))
cfg <- cfg0
cfg$frm <- ~grp + age

out <- init_state(mk_state(samps), cfg, minimal=TRUE)
report(identical(unname(out$config$covariate_types[["age"]]), "numeric"),
  "undeclared numeric covariate classified as numeric")

## a continuous covariate is not registered in AnalysisTableAnnotation$factors,
##   which holds categorical annotations only, but that registration is inert on
##   this path: LFQData$new() takes setup=FALSE, so setup_analysis(), which is
##   what would coerce a $factors column to character, never runs, and the plain
##   stats::lm() underneath build_model() reads the column as continuous:

res <- try(test_prolfqua(out$state, out$config, is_log_transformed=FALSE),
  silent=TRUE)
report(!inherits(res, "try-error"),
  "test_prolfqua() accepts a continuous covariate as an adjustment")

if(!inherits(res, "try-error")) {
  report(!("age" %in% names(lm1(res)$xlevels)),
    "continuous covariate is not turned into a factor by the fit")
  report("age" %in% names(stats::coef(lm1(res))),
    "continuous covariate gets a single slope coefficient")
  report(is.numeric(lm1(res)$model$age),
    "continuous covariate reaches the fit as numeric")
  report(nrow(res$hits) > 0 && all(res$hits$factor %in% "grp"),
    "hits still come from the tested term, not the adjustment")
}

## and as the tested term: one term of one degree of freedom, which is what
##   test_prolfqua() was bounded to while it read a per-term anova table:

cfg_age <- cfg
cfg_age$test_term <- "age"
out <- init_state(mk_state(samps), cfg_age, minimal=TRUE)
res <- try(test_prolfqua(out$state, out$config, is_log_transformed=FALSE),
  silent=TRUE)
report(!inherits(res, "try-error"),
  "test_prolfqua() accepts a continuous covariate as the tested term")

if(!inherits(res, "try-error")) {
  report(nrow(res$hits) > 0 && all(res$hits$factor %in% "age"),
    "anova hits returned for the continuous test_term")
  report(all(res$hits$Df %in% 1),
    "continuous test_term is a single one degree of freedom term")
  report(all(is.finite(res$hits$FDR)),
    "continuous test_term yields usable p-values")
}

## a numeric covariate declared in config$reference_levels is categorical, so
##   it is accepted; few enough levels to leave the design estimable:
samps <- samps0
samps$dose <- rep(c(2, 4, 6), 4)
cfg <- cfg0
cfg$frm <- ~grp + dose
cfg$reference_levels <- c(grp="trt", dose="4")

out <- init_state(mk_state(samps), cfg, minimal=TRUE)
res <- try(test_prolfqua(out$state, out$config, is_log_transformed=FALSE),
  silent=TRUE)
report(!inherits(res, "try-error"),
  "numeric covariate declared in reference_levels is accepted")
if(!inherits(res, "try-error")) {
  report(all(c("dose2", "dose6") %in% cols(res)) && !("dose4" %in% cols(res)),
    "declared reference level is the absent column for numeric-as-factor")
}

###############################################################################
section("continuous covariate inside an interaction")

## an interaction between a factor and a continuous covariate is one term of one
##   degree of freedom, so naming the interaction itself is within the one term at
##   a time bound and runs today:

samps <- samps0
samps$age <- seq(2, 24, length.out=nrow(samps0))
cfg <- cfg0
cfg$frm <- ~grp * age
cfg$test_term <- "grp:age"

out <- init_state(mk_state(samps), cfg, minimal=TRUE)
res <- try(test_prolfqua(out$state, out$config, is_log_transformed=FALSE),
  silent=TRUE)
report(!inherits(res, "try-error"),
  "test_prolfqua() accepts a factor by continuous interaction as test_term")

if(!inherits(res, "try-error")) {
  report(nrow(res$hits) > 0 && all(res$hits$factor %in% "grp:age"),
    "hits come from the interaction term")
  report(all(res$hits$Df %in% 1),
    "factor by continuous interaction is one degree of freedom")
  ## make.names() over the design column names turns 'grpctl:age' into
  ##   'grpctl.age', since the fit is given syntactic column names:

  report(all(c("grpctl", "age", "grpctl.age") %in% names(stats::coef(lm1(res)))),
    "both main effects and the interaction are in the fit")
  report("grpctl:age" %in% cols(res),
    "design$X keeps the unmangled interaction column name")
}

## naming the main effect instead is a joint test over grp and grp:age by
##   marginality. The per-term anova table could not express that and refused it;
##   the full versus reduced comparison runs it, over both columns at once:

cfg_main <- cfg
cfg_main$test_term <- "grp"
out <- init_state(mk_state(samps), cfg_main, minimal=TRUE)
res <- try(test_prolfqua(out$state, out$config, is_log_transformed=FALSE),
  silent=TRUE)
report(!inherits(res, "try-error"),
  "test_prolfqua() runs a main effect that an interaction contains")

if(!inherits(res, "try-error")) {
  report(all(res$hits$Df %in% 2),
    "the joint test carries 2 degrees of freedom")
  report(identical(sort(cols(res)[res$design$cols_test]),
    sort(c("grpctl", "grpctl:age"))),
    "the columns tested are the main effect and the interaction")
  report(all(res$hits$factor %in% "grp"),
    "rows are labelled with config$test_term, not the columns")

  ## and it is the F-test stats::anova() gives for the same nested pair: dropping
  ##   grpctl and grpctl:age leaves the intercept and age. Against p.value.unmod,
  ##   because the reported p.value has the error variance moderated across features
  ##   while stats::anova() uses each feature's own:

  f1 <- rownames(exprs)[1]
  d <- data.frame(y=exprs[f1, ], grp=out$state$samples$grp,
    age=out$state$samples$age)
  a <- stats::anova(stats::lm(y ~ age, data=d), stats::lm(y ~ grp * age, data=d))
  report(isTRUE(all.equal(res$hits$p.value.unmod[res$hits$pep %in% f1],
    a[2, "Pr(>F)"])),
    "the joint test matches stats::anova() on the same nested pair")
}

## the single-term case is unchanged by the switch: the Type I sum of squares for a
##   term entered last is the difference in residual sums of squares between these
##   same two fits, so the numbers are the ones the anova table used to report:

cfg_one <- cfg0
cfg_one$frm <- ~grp + age
cfg_one$test_term <- "age"
out <- init_state(mk_state(samps), cfg_one, minimal=TRUE)
res <- test_prolfqua(out$state, out$config, is_log_transformed=FALSE)

ref <- sapply(rownames(exprs), function(f) {
  d <- data.frame(y=exprs[f, ], grp=samps$grp, age=samps$age)
  stats::anova(stats::lm(y ~ grp + age, data=d))["age", "Pr(>F)"]
})
report(isTRUE(all.equal(res$hits$p.value.unmod, unname(ref[match(res$hits$pep,
  names(ref))]))),
  "a single-term test still matches the Type I row for that term entered last")

###############################################################################
section("the error variance is moderated across features")

## the reported test borrows variance across features the way test_trend() does:
##   prolfqua::squeezeVarRob(robust=FALSE) returns what limma::squeezeVar() returns,
##   so the reported statistic is limma's moderated F, at whatever numerator degrees
##   of freedom config$test_term carries. res is still the single-term ~grp + age
##   test of 'age' from just above, so one degree of freedom:

h <- res$hits
md <- res$fit$modelDF
dfr <- md$df.residual[match(h$pep, md$pep)]
sig2 <- md$sigma[match(h$pep, md$pep)]^2

report(all(c("moderated", "s2.denom", "df.denom", "df.prior", "F.value.unmod",
  "p.value.unmod") %in% names(h)),
  "hits carries the moderation columns")
report(all(h$moderated), "moderation is on by default")

sv <- prolfqua::squeezeVarRob(sig2, df=dfr, robust=FALSE)
report(isTRUE(all.equal(sv$var.post, limma::squeezeVar(sig2, df=dfr)$var.post)) &&
  isTRUE(all.equal(sv$df.prior, limma::squeezeVar(sig2, df=dfr)$df.prior)),
  "prolfqua's shrinkage is limma::squeezeVar()'s, the one test_trend() uses")
report(isTRUE(all.equal(h$s2.denom, sv$var.post)),
  "the denominator variance is the posterior variance, not the feature's own")
report(all(h$s2.denom >= pmin(sig2, sv$var.prior)) &&
  all(h$s2.denom <= pmax(sig2, sv$var.prior)),
  "the posterior variance lies between the feature's own and the prior")
report(all(h$df.prior > 0) && all(is.finite(h$df.prior)),
  "the prior contributes finite degrees of freedom")
report(isTRUE(all.equal(h$df.denom, dfr + h$df.prior)),
  "the denominator degrees of freedom are the residual df plus the prior df")
report(max(abs(h$p.value - h$p.value.unmod)) > 0.01,
  "moderating moves the p-values, so these checks are not vacuous")

## the external check at one degree of freedom: limma::eBayes() on the same design
##   and the same data must give the same p-value, and the moderated F must be the
##   square of the moderated t. The fixture has no NAs, so every feature has the
##   same design and limma can be handed the matrix whole:

X <- res$design$X
rownames(X) <- as.character(out$state$samples$obs)
Y <- exprs[h$pep, rownames(X), drop=FALSE]
lf <- limma::eBayes(limma::lmFit(Y, X))
tt <- limma::topTable(lf, coef=res$design$cols_test, number=Inf, sort.by="none")

report(isTRUE(all.equal(h$p.value, tt$P.Value)),
  "at one df the moderated F-test is limma's moderated t-test")
report(isTRUE(all.equal(h$F.value, unname(lf$t[, res$design$cols_test]^2))),
  "and the moderated F is the square of the moderated t")

## above one degree of freedom the same construction applies, which is what makes a
##   switch on df unnecessary. Checked against a moderated Wald F built from the
##   coefficients and (X'X)^-1, which shares no arithmetic with the difference of
##   residual sums of squares that f.nested_f() uses:

wald_p <- function(Y, X, cols_test) {
  q <- length(cols_test)
  V <- solve(t(X) %*% X)[cols_test, cols_test, drop=FALSE]
  num <- rep(NA_real_, nrow(Y))
  s2 <- num
  dfr <- num
  for(i in 1:nrow(Y)) {
    fit <- stats::lm.fit(X, Y[i, ])
    dfr[i] <- length(fit$residuals) - fit$rank
    s2[i] <- sum(fit$residuals^2) / dfr[i]
    b <- fit$coefficients[cols_test]
    num[i] <- as.numeric(t(b) %*% solve(V) %*% b) / q
  }
  sv <- limma::squeezeVar(s2, df=dfr)
  stats::pf(num / sv$var.post, q, dfr + sv$df.prior, lower.tail=FALSE)
}

out2 <- init_state(mk_state(samps), cfg_main, minimal=TRUE)     ## ~grp*age, 'grp'
res2 <- test_prolfqua(out2$state, out2$config, is_log_transformed=FALSE)
h2 <- res2$hits
X2 <- res2$design$X
rownames(X2) <- as.character(out2$state$samples$obs)

report(all(h2$Df %in% 2) && all(h2$moderated),
  "the joint test is moderated too, with no switch on degrees of freedom")
report(isTRUE(all.equal(h2$p.value,
  wald_p(exprs[h2$pep, rownames(X2), drop=FALSE], X2, res2$design$cols_test))),
  "the joint moderated F matches an independently computed moderated Wald F")
report(isTRUE(all.equal(h2$p.value,
  stats::pf(h2$F.value, h2$Df, h2$df.denom, lower.tail=FALSE))),
  "the reported p-value is that F against the reported denominator df")

## and moderation can be turned off, which is the only way to get the p-value
##   stats::anova() would give:

cfg_nm <- cfg_one
cfg_nm$test_moderate <- FALSE
out_nm <- init_state(mk_state(samps), cfg_nm, minimal=TRUE)
res_nm <- test_prolfqua(out_nm$state, out_nm$config, is_log_transformed=FALSE)
hn <- res_nm$hits

report(!any(hn$moderated), "config$test_moderate=FALSE turns moderation off")
report(isTRUE(all.equal(hn$p.value, hn$p.value.unmod)),
  "and the reported test is then the unmoderated one")
report(all(is.na(hn$df.prior)) && isTRUE(all.equal(hn$s2.denom,
  md$sigma[match(hn$pep, md$pep)]^2)),
  "with the feature's own variance and no prior degrees of freedom")
report(isTRUE(all.equal(hn$p.value.unmod,
  h$p.value.unmod[match(hn$pep, h$pep)])),
  "moderating does not change the unmoderated column")

## config plumbing for the new key:

report(isTRUE(new_config()$test_moderate), "new_config() moderates by default")
report(isTRUE(check_config(list(test_moderate=FALSE, log_file=log_file))),
  "check_config() accepts a logical test_moderate")
report(threw(check_config(list(test_moderate="no", log_file=log_file))),
  "check_config() rejects a non-logical test_moderate")

###############################################################################
section("the prior variance can be fitted against mean intensity")

## config$test_trend fits the prior against each feature's mean intensity instead of
##   shrinking every feature toward one number. cfg_one is still the one degree of
##   freedom ~grp + age test of 'age', so limma can be asked for the same test:

cfg_tr <- cfg_one
cfg_tr$test_trend <- TRUE
out_tr <- init_state(mk_state(samps), cfg_tr, minimal=TRUE)
res_tr <- test_prolfqua(out_tr$state, out_tr$config, is_log_transformed=FALSE)
ht <- res_tr$hits

report(all(c("trend", "s2.prior") %in% names(ht)),
  "hits carries the columns describing the prior")
report(all(ht$moderated) && all(ht$trend),
  "config$test_trend fits the prior against the covariate")
report(!any(h$trend) && all(h$moderated),
  "and defaults to off, leaving the flat prior of the moderation above")
report(length(unique(ht$s2.prior)) > 1,
  "the prior variance varies by feature under the trend")
report(length(unique(h$s2.prior)) %in% 1,
  "and is one number shared by every feature without it")

## the covariate is limma's fit$Amean, the row mean of the response over the
##   observations where the feature was seen, and the shrinkage is again limma's:

amean <- rowMeans(exprs[ht$pep, , drop=FALSE], na.rm=TRUE)
mdt <- res_tr$fit$modelDF
dfr_t <- mdt$df.residual[match(ht$pep, mdt$pep)]
sig2_t <- mdt$sigma[match(ht$pep, mdt$pep)]^2
svt <- prolfqua::squeezeVarRob(sig2_t, df=dfr_t, covariate=amean, robust=FALSE)
svl <- limma::squeezeVar(sig2_t, df=dfr_t, covariate=amean)

report(isTRUE(all.equal(svt$var.post, svl$var.post)) &&
  isTRUE(all.equal(svt$df.prior, svl$df.prior)),
  "prolfqua's trended shrinkage is limma::squeezeVar(covariate=)'s")
report(isTRUE(all.equal(ht$s2.prior, as.numeric(svl$var.prior))),
  "the reported prior variance is that fitted trend, evaluated per feature")
report(isTRUE(all.equal(ht$s2.denom, svl$var.post)),
  "the denominator variance is the posterior variance under the trended prior")
report(isTRUE(all.equal(ht$df.denom, dfr_t + ht$df.prior)),
  "and the denominator df is the residual df plus the prior df of the trend")

## the external check, and the point of the option: limma::eBayes(trend=TRUE) is
##   what test_trend() runs, so with the trend on the two moderate identically. X and
##   Y are the design and data of the same test, from the section above:

lft <- limma::eBayes(limma::lmFit(Y, X), trend=TRUE)
ttt <- limma::topTable(lft, coef=res_tr$design$cols_test, number=Inf, sort.by="none")
i <- match(ht$pep, rownames(Y))

report(isTRUE(all.equal(unname(amean), unname(lft$Amean[i]))),
  "the covariate used is the one limma::lmFit() records as Amean")
report(isTRUE(all.equal(ht$p.value, ttt$P.Value[i])),
  "the trended moderated F-test is limma::eBayes(trend=TRUE)'s test")
report(isTRUE(all.equal(ht$F.value, unname(lft$t[i, res_tr$design$cols_test]^2))),
  "and its F is the square of that trended moderated t")
report(max(abs(ht$p.value - h$p.value[match(ht$pep, h$pep)])) > 1e-6,
  "the trend moves the p-values, so it is not a no-op here")
report(isTRUE(all.equal(ht$p.value.unmod, h$p.value.unmod[match(ht$pep, h$pep)])),
  "the trend does not touch the unmoderated column")

## a spline of up to four degrees of freedom needs more than four features to fit,
##   and prolfqua answers too few with an all-NA prior rather than by falling back,
##   which would lose every p-value. Caught before that, leaving the flat prior:

cfg_few <- cfg_tr
few <- list(expression=exprs[1:4, , drop=FALSE], features=feats[1:4, , drop=FALSE],
  samples=samps)
out_few <- init_state(few, cfg_few, minimal=TRUE)
res_few <- test_prolfqua(out_few$state, out_few$config, is_log_transformed=FALSE)

cfg_few2 <- cfg_few
cfg_few2$test_trend <- FALSE
out_few2 <- init_state(few, cfg_few2, minimal=TRUE)
res_few2 <- test_prolfqua(out_few2$state, out_few2$config, is_log_transformed=FALSE)

report(nrow(res_few$hits) %in% 4 && !any(is.na(res_few$hits$p.value)),
  "too few features to fit the trend still yields a p-value for each")
report(!any(res_few$hits$trend) && all(res_few$hits$moderated),
  "with the trend abandoned but the moderation kept")
report(isTRUE(all.equal(res_few$hits$p.value, res_few2$hits$p.value)),
  "and gives exactly what asking for the flat prior gives")

## and the trend sets the prior of a shrinkage that config$test_moderate can turn
##   off, in which case there is no prior to set:

cfg_both <- cfg_tr
cfg_both$test_moderate <- FALSE
out_both <- init_state(mk_state(samps), cfg_both, minimal=TRUE)
res_both <- test_prolfqua(out_both$state, out_both$config, is_log_transformed=FALSE)
hb <- res_both$hits

report(!any(hb$moderated) && !any(hb$trend),
  "config$test_moderate=FALSE overrides config$test_trend")
report(isTRUE(all.equal(hb$p.value, hb$p.value.unmod)) && all(is.na(hb$s2.prior)),
  "and reports the unmoderated test, with no prior")

## config plumbing for the new key:

report(identical(new_config()$test_trend, FALSE),
  "new_config() leaves the trend off")
report(isTRUE(check_config(list(test_trend=TRUE, log_file=log_file))),
  "check_config() accepts a logical test_trend")
report(threw(check_config(list(test_trend="yes", log_file=log_file))),
  "check_config() rejects a non-logical test_trend")

###############################################################################
section("continuous covariate under missingness")

## sim_design() with mnar_c0=-Inf drops nothing, so the fixture above has no NAs at all
##   and none of the assertions so far exercise missingness. Introduced here by
##   hand, so that what is missing and where is known:

report(sum(is.na(exprs)) %in% 0, "the shared fixture has no missing values")

## age tied at one value in three samples of each group, so that a feature
##   observed only there still supports grp but leaves age constant, and so
##   aliased, among the samples where it is observed:

age_tied <- c(rep(5, 3), 7, 9, 11, rep(5, 3), 7, 9, 11)
samps_m <- samps0
samps_m$age <- age_tied
cfg_m <- cfg0
cfg_m$frm <- ~grp + age
cfg_m$test_term <- "age"

e <- exprs
f_alias <- rownames(e)[1]
f_short <- rownames(e)[2]
f_full <- rownames(e)[3]
e[f_alias, !(age_tied %in% 5)] <- NA          ## observed only where age == 5
e[f_short, c(3, 9)] <- NA                     ## merely two observations short

report(sum(is.na(e)) > 0, "the missingness fixture does have missing values")

out <- init_state(list(expression=e, features=feats, samples=samps_m), cfg_m,
  minimal=TRUE)
m_miss <- mark()
res <- try(test_prolfqua(out$state, out$config, is_log_transformed=FALSE),
  silent=TRUE)
report(!inherits(res, "try-error"),
  "test_prolfqua() runs with a continuous covariate and missing values")

if(!inherits(res, "try-error")) {

  h <- res$hits
  md <- res$fit$modelDF
  dfr <- function(f) md$df.residual[md$pep %in% f]

  report(f_short %in% h$pep,
    "a feature merely short some observations is still tested")
  report(dfr(f_short) < dfr(f_full),
    "its residual degrees of freedom drop with the observations it lost")
  report(is.finite(h$p.value[h$pep %in% f_short]),
    "and it still yields a usable p-value")

  ## the aliased feature: stats::lm() returns NA for the slope it cannot estimate,
  ##   so removing the tested column does not reduce the rank of the fit and the
  ##   hypothesis is not estimable. Dropped explicitly and reported, rather than
  ##   vanishing because stats::anova() omitted the term (which is what happened
  ##   when the per-term table was read) or surfacing with an NaN F:

  report(md$isSingular[md$pep %in% f_alias],
    "a feature whose continuous covariate is constant is flagged singular")
  report(is.na(stats::coef(md$linear_model[[which(md$pep %in% f_alias)]])[["age"]]),
    "its slope is NA rather than a number")
  report(!(f_alias %in% h$pep), "and it is absent from the hits")
  report(logged("do not support the test of config$test_term", since=m_miss),
    "the drop is reported rather than silent")
  report(logged("WARNING", since=m_miss),
    "the drop is reported as a warning")
  report(logged(f_alias, since=m_miss),
    "the report names the dropped feature")
  report(length(unique(h$pep)) %in% (nrow(e) - 1),
    "every other feature is still tested")
}

###############################################################################

cat("\n## passes:", n_pass, "; failures:", n_fail, "; elapsed:",
  round(as.numeric(difftime(Sys.time(), t0, units="secs"))), "s\n")
cat("## log file:", log_file, "\n")

quit(save="no", status=if(n_fail > 0) 1 else 0)
