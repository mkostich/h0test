## Null calibration of the candidate msqrob2::msqrobAggregate() path, run before the
##   path is written into the package. The question is whether msqrob2's denominator
##   degrees of freedom, which are an effective residual df for the whole feature level
##   fit rather than a Satterthwaite ddf for the contrast, leave the test
##   anti-conservative in the way the first draft of test_method="prolfqua_lmer" was.
##
##   Three simulated truths, all with no group effect:
##     none     residual noise only.
##     shared   one per-observation shift vector, the same for every gene, which is
##              what a run or batch effect is. Note that this does not leave a true
##              null: the realized difference between the group means of that vector
##              is a real offset carried by every gene alike, so a rejection rate
##              above the level is not by itself a calibration failure. Reported
##              anyway, with that offset, because it is the case a reader will assume.
##     by_gene  an independent per-observation shift vector per gene, which is the
##              variance component (1|obs_col) actually models when the fit is per
##              gene. This is the honest null for the question being asked.
##
##   Arms: the existing test_method="msqrob" path (aggregate, then fit) as a
##   reference; msqrobAggregate with ridge off and on, each with and without the
##   random observation effect; and lmerTest, which is the engine
##   test_method="prolfqua_lmer" uses, so that the step 2 choice of
##   config$test_random_obs=TRUE is re-checked against the same fixture.
##
## Usage: Rscript calib_msqrob_agg.R [n_genes] [n_peps] [n_obs] > calib.out 2>&1

suppressMessages(library(msqrob2))
suppressMessages(library(QFeatures))
suppressMessages(library(lmerTest))
BiocParallel::register(BiocParallel::SerialParam())   ## bplapply(); see the segfault
                                                     ##   note for the lme4 stack
args <- commandArgs(trailingOnly=TRUE)
n_gene <- if(length(args) >= 1) as.integer(args[1]) else 300
n_pep <- if(length(args) >= 2) as.integer(args[2]) else 5
n_obs <- if(length(args) >= 3) as.integer(args[3]) else 12
alpha <- 0.05
set.seed(101)

cat("## genes:", n_gene, " peptides/gene:", n_pep, " observations:", n_obs,
  " alpha:", alpha, "\n\n")

feats <- data.frame(
  pep=sprintf("p%05d", 1:(n_gene * n_pep)),
  gene=rep(sprintf("g%04d", 1:n_gene), each=n_pep),
  stringsAsFactors=FALSE
)
samps <- data.frame(
  obs=sprintf("s%02d", 1:n_obs),
  grp=rep(c("ctl", "trt"), each=n_obs / 2),
  sex=rep(c("F", "M"), times=n_obs / 2),
  stringsAsFactors=FALSE
)
trt <- samps$grp %in% "trt"

sim <- function(run_sd, mode) {
  base <- rnorm(nrow(feats), 20, 1)
  e <- matrix(rnorm(nrow(feats) * n_obs, 0, 1), nrow=nrow(feats))
  y <- e + base
  off <- NA_real_
  if(mode %in% "shared") {
    s <- rnorm(n_obs, 0, run_sd)
    y <- y + rep(s, each=nrow(feats))
    off <- mean(s[trt]) - mean(s[!trt])
  } else if(mode %in% "by_gene") {
    s <- matrix(rnorm(n_gene * n_obs, 0, run_sd), nrow=n_gene)
    y <- y + s[rep(1:n_gene, each=n_pep), ]
    off <- stats::sd(rowMeans(s[, trt, drop=F]) - rowMeans(s[, !trt, drop=F]))
  }
  dimnames(y) <- list(feats$pep, samps$obs)
  return(list(y=y, off=off))
}

mk <- function(y) {
  x <- as.data.frame(y)
  x <- cbind(fnames=rownames(x), x)
  rownames(x) <- NULL
  o <- QFeatures::readQFeatures(assayData=x, quantCols=2:ncol(x), fnames="fnames",
    name="features", verbose=FALSE)
  for(nm in names(feats)) {
    SummarizedExperiment::rowData(o[["features"]])[[nm]] <- feats[[nm]]
  }
  SummarizedExperiment::colData(o)$grp <- factor(samps$grp, levels=c("ctl", "trt"))
  SummarizedExperiment::colData(o)$sex <- factor(samps$sex, levels=c("F", "M"))
  SummarizedExperiment::colData(o)$obs <- factor(samps$obs)
  return(o)
}

## the aggregator only fills the values of the aggregated assay, which take no part in
##   the fit: msqrobAggregate() fits on the un-aggregated assay. colMeans avoids
##   MsCoreUtils::robustSummary(); see f.combine_features_robust_summary():

agg_fun <- function(x, ...) base::colMeans(x, na.rm=TRUE)

say <- function(lbl, p, df, secs=NA) {
  cat(sprintf("%-44s reject %.3f  (%d/%d)  median df %s%s\n", lbl,
    mean(p < alpha, na.rm=TRUE), sum(p < alpha, na.rm=TRUE), sum(!is.na(p)),
    format(stats::median(df, na.rm=TRUE), digits=4),
    if(is.na(secs)) "" else sprintf("  %.0fs", secs)))
}

## the moderated t msqrob2::hypothesisTest() reports, computed from the same fitted
##   StatModel so that the ridge renamed parameters can be reached by name:

msq_p <- function(models, nom) {
  p <- rep(NA_real_, length(models))
  d <- rep(NA_real_, length(models))
  ## msqrobLmer() reports "lmer", msqrobLm() "lm" or "rlm", and a fit that failed
  ##   "fitError", so select by what did not fail rather than by what did:
  ok <- !vapply(models, msqrob2::getFitMethod, character(1)) %in% "fitError"
  for(i in which(ok)) {
    b <- msqrob2::getCoef(models[[i]])
    if(!(nom %in% names(b))) next
    L <- matrix(0, nrow=length(b), ncol=1, dimnames=list(names(b), nom))
    L[nom, 1] <- 1
    v <- drop(t(L) %*% msqrob2::getVcovUnscaled(models[[i]]) %*% L)
    se <- sqrt(msqrob2::getSigmaPosterior(models[[i]])^2 * v)
    d[i] <- msqrob2::getDfPosterior(models[[i]])
    p[i] <- 2 * stats::pt(-abs(b[nom] / se), df=d[i])
  }
  return(list(p=p, df=d))
}

cell <- function(lbl, y, frm, ridge, nom) {
  t0 <- Sys.time()
  o <- try(suppressMessages(suppressWarnings(
    msqrob2::msqrobAggregate(mk(y), i="features", fcol="gene", name="genes",
      formula=frm, ridge=ridge, aggregateFun=agg_fun))), silent=TRUE)
  if(inherits(o, "try-error")) {
    cat(sprintf("%-44s ERROR: %s\n", lbl,
      trimws(conditionMessage(attr(o, "condition")))))
    return(invisible(NULL))
  }
  r <- msq_p(SummarizedExperiment::rowData(o[["genes"]])$msqrobModels, nom)
  say(lbl, r$p, r$df, as.numeric(difftime(Sys.time(), t0, units="secs")))
}

## the reference: aggregate first, then fit at gene level, which is what
##   test_method="msqrob" does today:

ref <- function(lbl, y) {
  t0 <- Sys.time()
  o <- try(suppressMessages(suppressWarnings({
    o <- QFeatures::aggregateFeatures(mk(y), i="features", fcol="gene", na.rm=TRUE,
      name="genes", fun=base::colMeans)
    msqrob2::msqrob(object=o, i="genes", formula=~grp + sex, maxitRob=100)
  })), silent=TRUE)
  if(inherits(o, "try-error")) {
    cat(sprintf("%-44s ERROR: %s\n", lbl,
      trimws(conditionMessage(attr(o, "condition")))))
    return(invisible(NULL))
  }
  r <- msq_p(SummarizedExperiment::rowData(o[["genes"]])$msqrobModels, "grptrt")
  say(lbl, r$p, r$df, as.numeric(difftime(Sys.time(), t0, units="secs")))
}

## the engine test_method="prolfqua_lmer" uses: one lmerTest fit per gene, Satterthwaite
##   denominator degrees of freedom for the tested coefficient:

lmt <- function(lbl, y, obs_term) {
  t0 <- Sys.time()
  frm <- if(obs_term) y ~ grp + sex + (1 | pep) + (1 | obs) else
    y ~ grp + sex + (1 | pep)
  p <- rep(NA_real_, n_gene)
  d <- rep(NA_real_, n_gene)
  for(g in 1:n_gene) {
    i <- which(feats$gene %in% sprintf("g%04d", g))
    dat <- data.frame(
      y=as.vector(y[i, ]),
      pep=rep(feats$pep[i], times=n_obs),
      obs=rep(samps$obs, each=length(i)),
      grp=rep(samps$grp, each=length(i)),
      sex=rep(samps$sex, each=length(i)),
      stringsAsFactors=TRUE
    )
    fit <- try(suppressMessages(suppressWarnings(lmerTest::lmer(frm, data=dat))),
      silent=TRUE)
    if(inherits(fit, "try-error")) next
    cf <- try(suppressWarnings(summary(fit)$coefficients), silent=TRUE)
    if(inherits(cf, "try-error") || !("grptrt" %in% rownames(cf))) next
    p[g] <- cf["grptrt", "Pr(>|t|)"]
    d[g] <- cf["grptrt", "df"]
  }
  say(lbl, p, d, as.numeric(difftime(Sys.time(), t0, units="secs")))
}

run <- function(mode, run_sd) {

  s <- sim(run_sd, mode)
  cat("## truth:", mode, " per-observation shift sd:", run_sd,
    "(residual sd is 1)\n")
  if(mode %in% "shared") {
    cat("##   realized offset between the group means of the shared shift:",
      round(s$off, 3), "-- carried by every gene, so it is a real effect\n")
  } else if(mode %in% "by_gene") {
    cat("##   sd across genes of the offset between group means:",
      round(s$off, 3), "-- independent per gene, so the null holds\n")
  }
  y <- s$y

  ref("msqrob: aggregate, then gene level fit", y)
  cell("msqrobAggregate ridge=F (1|pep)", y, ~grp + sex + (1 | pep), FALSE, "grptrt")
  cell("msqrobAggregate ridge=F (1|pep)+(1|obs)", y,
    ~grp + sex + (1 | pep) + (1 | obs), FALSE, "grptrt")
  cell("msqrobAggregate ridge=T (1|pep)", y, ~grp + sex + (1 | pep), TRUE,
    "ridgegrptrt")
  cell("msqrobAggregate ridge=T (1|pep)+(1|obs)", y,
    ~grp + sex + (1 | pep) + (1 | obs), TRUE, "ridgegrptrt")
  lmt("lmerTest (1|pep)", y, FALSE)
  lmt("lmerTest (1|pep)+(1|obs)", y, TRUE)

  cat("\n")
}

run("none", 0)
run("by_gene", 1)
run("shared", 1)

cat("## done\n")
