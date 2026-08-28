#!/usr/bin/env Rscript

## msqrob_agg_reparam.R: localize the msqrob_agg shift between two parameterizations of the
##   same hypothesis, flagged in 1/contrast_matrix.R and summary entries 7 and 8.
##
## The hypothesis: ~sex with config$test_term "sex" and ~0 + sex with config$contrast
##   "sexM - sexF" state the same thing, the two designs spanning the same column space, so
##   every method's fixed effects should be invariant. Eight of nine methods are, to 1e-12
##   or better; msqrob_agg moves 0.55% in the p-value and 0.14% in the effect size.
##
## Three things are separated here, in order:
##   1. which reported quantity moves: the weighted sum of coefficients, the unscaled
##      standard error, the raw sigma and residual df, or the moderated sigma and df that
##      limma::squeezeVar() produces across genes inside msqrob2;
##   2. whether the genes that move are the ones h0testr fits at the gene level, which have
##      one observed feature and go to MASS::rlm() through msqrob2::msqrob(), or the ones
##      fitted as mixed models by lme4 through msqrob2::msqrobAggregate();
##   3. whether the movement is the stopping rule of an iterative fit rather than a
##      difference of estimate: each fit is repeated with the optimizer and the robust
##      iteration tolerances tightened, and the movement re-measured.
##
## Nothing in the package is modified. Reads h0testr, msqrob2, lme4 and MASS only.

suppressMessages(library(msqrob2))
suppressMessages(library(lme4))

args <- commandArgs(trailingOnly=TRUE)

if(length(args) != 2) {
  msg <- c(
    "",
    "msqrob_agg_reparam.R: localize the msqrob_agg deviation between ~sex with",
    "  config$test_term 'sex' and ~0 + sex with config$contrast 'sexM - sexF', two",
    "  parameterizations of one hypothesis whose fixed effects should be invariant.",
    "  Separates which reported quantity moves, which genes move, and whether the",
    "  movement is a stopping rule of an iterative fit or a difference of estimate.",
    "",
    "usage: Rscript msqrob_agg_reparam.R <r_dir> <outprefix>",
    "",
    "required positional arguments:",
    "  <r_dir>       path to the h0testr package R/ source directory; all .R files there",
    "                  are sourced, as in 1/contrast_matrix.R; the installed package is",
    "                  not used",
    "  <outprefix>   path prefix for the outputs; <outprefix>.txt is the report, the",
    "                  same text written to stdout, and <outprefix>.log the h0testr log",
    "",
    "outputs:",
    "  <outprefix>.txt   report: per-quantity deviations, per-gene deviations with the",
    "                      fit route each gene took, and the same deviations recomputed",
    "                      with the fit tolerances tightened",
    "  <outprefix>.log   h0testr's own messages from the two runs",
    "",
    "examples:",
    "  Rscript msqrob_agg_reparam.R C:/path/to/h0test/h0testr/R out/msqrob_reparam",
    "  Rscript msqrob_agg_reparam.R ../../h0test/h0testr/R ./scratch/agg1",
    "  Rscript msqrob_agg_reparam.R C:/path/to/h0testr/R ./reparam > reparam.log 2>&1",
    ""
  )
  cat(paste(msg, collapse="\n"), "\n", file=stderr(), sep="")
  quit(save="no", status=2)
}

r_dir <- args[1]

if(!dir.exists(r_dir)) {
  cat("ERROR: msqrob_agg_reparam.R: r_dir is not a directory: ", r_dir, "\n", sep="",
    file=stderr())
  quit(save="no", status=2)
}

for(f in list.files(r_dir, pattern="[.]R$", full.names=TRUE)) source(f)

outprefix <- args[2]
out_file <- paste0(outprefix, ".txt")
log_file <- paste0(outprefix, ".log")

dir_out <- dirname(out_file)

if(!dir.exists(dir_out)) {
  cat("ERROR: msqrob_agg_reparam.R: output directory does not exist: ", dir_out, "\n",
    sep="", file=stderr())
  quit(save="no", status=1)
}

unlink(c(out_file, log_file))

rpt <- character(0)

emit <- function(...) {
  line <- paste0(...)
  rpt <<- c(rpt, line)
  cat(line, "\n", sep="")
  return(invisible(NULL))
}

t0 <- Sys.time()

progress <- function(...) {
  cat("[", format(round(as.numeric(difftime(Sys.time(), t0, units="secs")))), "s] ",
    paste0(...), "\n", sep="", file=stderr())
  return(invisible(NULL))
}

## relative deviation between two vectors, on the scale of the larger, NA-safe. Used for
##   quantities with a meaningful zero-free scale (p-values, sigma, df):

rel_dev <- function(a, b) {
  ok <- is.finite(a) & is.finite(b)
  if(!any(ok)) return(NA_real_)
  den <- pmax(abs(a[ok]), abs(b[ok]))
  den[den < .Machine$double.eps] <- 1
  return(max(abs(a[ok] - b[ok]) / den))
}

## absolute deviation scaled by the spread of the quantity, for effect sizes, which are
##   centered on zero on a null fixture so a relative deviation is meaningless there:

abs_dev <- function(a, b) {
  ok <- is.finite(a) & is.finite(b)
  if(!any(ok)) return(NA_real_)
  return(max(abs(a[ok] - b[ok])) / max(1, max(abs(c(a[ok], b[ok])))))
}

###############################################################################
## fixture: identical to 1/capability_matrix.R and 1/contrast_matrix.R, so that anything
##   found here is comparable cell by cell with what those tables report.

set.seed(101)

nobs <- 24
ngene <- 40
npep_max <- 6

samps <- expand.grid(rep=1:2, batch=c("b1", "b2"), sex=c("F", "M"),
  grp=c("a", "b", "c"), stringsAsFactors=FALSE)
samps <- data.frame(
  obs=paste0("o", sprintf("%02d", 1:nobs)),
  grp=samps$grp,
  sex=samps$sex,
  batch=samps$batch,
  age=round(rnorm(nobs, 50, 10), 1),
  stringsAsFactors=FALSE
)

npep <- ((seq_len(ngene) - 1) %% npep_max) + 1          ## 1..6 peptides per gene
gene <- rep(paste0("g", sprintf("%02d", seq_len(ngene))), times=npep)
exprs <- matrix(rnorm(length(gene) * nobs, 20, 1.5), nrow=length(gene))
rownames(exprs) <- paste0(gene, ".p", unlist(lapply(npep, seq_len)))
colnames(exprs) <- samps$obs
feats <- data.frame(pep=rownames(exprs), gene=gene)

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
cfg0$is_log_transformed <- TRUE
cfg0$impute_method <- "none"
cfg0$n_features_min <- 5
cfg0$log_file <- log_file
cfg0$test_prior_df <- 5

state_raw <- list(expression=exprs, features=feats, samples=samps)

## the two runs. A is the config$test_term route, B the config$contrast route; L is the
##   weighted sum each states, over the fixed-effect names of its own design:

runs <- list(
  A=list(frm=~sex, term="sex", con="", L=c(sexM=1), lab="~sex, test_term sex"),
  B=list(frm=~0 + sex, term="", con="sexM - sexF", L=c(sexM=1, sexF=-1),
    lab="~0 + sex, contrast 'sexM - sexF'")
)

###############################################################################
## step 1: the two h0testr runs, and what each reports per gene.

fits <- list()
tabs <- list()

for(nm in names(runs)) {

  progress("h0testr run ", nm, ": ", runs[[nm]]$lab)

  ## f.subset_covariates() requires config$reference_levels to name only variables the
  ##   formula uses, so it is narrowed to the one variable in play here:

  cf <- cfg0
  cf$reference_levels <- c(sex="F")
  cf$frm <- runs[[nm]]$frm
  cf$test_term <- runs[[nm]]$term
  cf$contrast <- runs[[nm]]$con

  ini <- try(suppressMessages(initialize(state_raw, cf, minimal=TRUE)), silent=TRUE)

  if(inherits(ini, "try-error")) {
    cat("ERROR: msqrob_agg_reparam.R: initialize() failed for run ", nm, ": ",
      trimws(conditionMessage(attr(ini, "condition"))), "\n", sep="", file=stderr())
    quit(save="no", status=1)
  }

  res <- try(suppressMessages(test(ini$state, ini$config, method="msqrob_agg",
    is_log_transformed=TRUE, prior_df=5)), silent=TRUE)

  if(inherits(res, "try-error")) {
    cat("ERROR: msqrob_agg_reparam.R: test() failed for run ", nm, ": ",
      trimws(conditionMessage(attr(res, "condition"))), "\n", sep="", file=stderr())
    quit(save="no", status=1)
  }

  tabs[[nm]] <- res$standard[order(res$standard$feature), , drop=FALSE]
  fits[[nm]] <- res$fit
}

genes <- tabs$A$feature
stopifnot(identical(genes, tabs$B$feature))

n_pep <- as.vector(table(feats$gene)[genes])

emit("")
emit("msqrob_agg reparameterization shift: where it comes from")
emit(strrep("=", 78))
emit("")
emit("Fixture: ", nobs, " observations, ", ngene, " genes, 1-6 peptides per gene, all",
  " values simulated")
emit("  under the null, seed 101; the fixture 1/capability_matrix.R and",
  " 1/contrast_matrix.R use.")
emit("Run A: ", runs$A$lab)
emit("Run B: ", runs$B$lab)
emit("")
emit("Step 1. What test() reports, over all ", length(genes), " genes:")
emit("")
emit(sprintf("  %-28s %12s", "quantity", "deviation"))
emit(sprintf("  %-28s %12s", strrep("-", 28), strrep("-", 12)))
emit(sprintf("  %-28s %12.3g", "logfc (absolute, scaled)", abs_dev(tabs$A$logfc,
  tabs$B$logfc)))
emit(sprintf("  %-28s %12.3g", "stat (relative)", rel_dev(tabs$A$stat, tabs$B$stat)))
emit(sprintf("  %-28s %12.3g", "pval (relative)", rel_dev(tabs$A$pval, tabs$B$pval)))

###############################################################################
## step 2: the same comparison inside the fit, quantity by quantity. msqrob2 stores a
##   StatModel per gene: getCoef() the fixed effects followed by the BLUPs, getVcovUnscaled()
##   their joint unscaled covariance, getSigma() and getDF() the raw scale and residual df,
##   and getSigmaPosterior() and getDfPosterior() the same two after limma::squeezeVar()
##   has moderated them across the genes of the assay. hypothesisTest() forms
##   L'beta / (sigma_post * sqrt(L'VL)) against dfPosterior, so those five quantities are
##   everything the reported t and p are made of.

grab <- function(fit, L, genes) {

  se <- fit[["genes"]]
  dat <- SummarizedExperiment::rowData(se)
  models <- dat$msqrobModels
  ids <- as.character(dat[["gene"]])

  if(is.null(models)) return(NULL)

  idx <- match(genes, ids)
  nm_L <- names(L)

  one <- function(i) {
    m <- models[[i]]
    empty <- list(eff=NA_real_, sd_unscaled=NA_real_, sigma=NA_real_, df=NA_real_,
      sigma_post=NA_real_, df_post=NA_real_, type=NA_character_)
    if(is.na(i)) return(empty)
    ty <- msqrob2::getFitMethod(m)
    b <- try(msqrob2::getCoef(m), silent=TRUE)
    if(inherits(b, "try-error") || !all(nm_L %in% names(b))) {
      empty$type <- ty
      return(empty)
    }
    v <- msqrob2::getVcovUnscaled(m)[nm_L, nm_L, drop=FALSE]
    return(list(
      eff=sum(L * b[nm_L]),
      sd_unscaled=sqrt(drop(t(L) %*% v %*% L)),
      sigma=msqrob2::getSigma(m),
      df=msqrob2::getDF(m),
      sigma_post=msqrob2::getSigmaPosterior(m),
      df_post=msqrob2::getDfPosterior(m),
      type=ty
    ))
  }

  parts <- lapply(seq_along(idx), function(k) one(idx[k]))

  return(data.frame(
    gene=genes,
    eff=vapply(parts, function(x) x$eff, numeric(1)),
    sd_unscaled=vapply(parts, function(x) x$sd_unscaled, numeric(1)),
    sigma=vapply(parts, function(x) x$sigma, numeric(1)),
    df=vapply(parts, function(x) x$df, numeric(1)),
    sigma_post=vapply(parts, function(x) x$sigma_post, numeric(1)),
    df_post=vapply(parts, function(x) x$df_post, numeric(1)),
    type=vapply(parts, function(x) x$type, character(1)),
    stringsAsFactors=FALSE
  ))
}

progress("extracting fitted quantities")

pa <- grab(fits$A, runs$A$L, genes)
pb <- grab(fits$B, runs$B$L, genes)

emit("")

if(is.null(pa) || is.null(pb)) {

  emit("Step 2. The fitted models could not be read back out of test()$fit, so the")
  emit("  localization below is not available. Nothing is concluded from that.")

} else {

  quants <- list(
    list(k="eff", lab="L'beta, the weighted sum", fn=abs_dev),
    list(k="sd_unscaled", lab="sqrt(L'VL), unscaled se", fn=rel_dev),
    list(k="sigma", lab="sigma, raw", fn=rel_dev),
    list(k="df", lab="residual df, raw", fn=rel_dev),
    list(k="sigma_post", lab="sigma, moderated", fn=rel_dev),
    list(k="df_post", lab="df, moderated", fn=rel_dev)
  )

  emit("Step 2. The same comparison inside the fit, quantity by quantity. A relative")
  emit("  deviation except for L'beta, which is centered on zero here:")
  emit("")
  emit(sprintf("  %-28s %12s", "quantity", "deviation"))
  emit(sprintf("  %-28s %12s", strrep("-", 28), strrep("-", 12)))

  for(q in quants) {
    emit(sprintf("  %-28s %12.3g", q$lab, q$fn(pa[[q$k]], pb[[q$k]])))
  }

  ## step 3: which genes move, and by which route they were fitted. h0testr fits a gene
  ##   with fewer than two observed features at the gene level instead, MASS::rlm() through
  ##   msqrob2::msqrob(), because its random feature effect is not identified; the rest are
  ##   lme4 mixed models. The two routes have different iterative machinery, so which one
  ##   the moving genes took is the next thing to settle:

  d_eff <- abs(pa$eff - pb$eff)
  d_p <- abs(tabs$A$pval - tabs$B$pval) / pmax(tabs$A$pval, tabs$B$pval)

  ord <- order(d_p, decreasing=TRUE)

  emit("")
  emit("Step 3. Which genes move, worst ten by p-value deviation, with the route each")
  emit("  took: 'lmer' is the mixed model msqrobAggregate() fits, 'rlm' the gene level")
  emit("  fit h0testr substitutes when a gene has fewer than two observed features:")
  emit("")
  emit(sprintf("  %-6s %5s %-6s %11s %11s %11s %11s", "gene", "npep", "route",
    "d(pval)", "d(L'beta)", "d(df_raw)", "d(sigma)"))
  emit(sprintf("  %-6s %5s %-6s %11s %11s %11s %11s", strrep("-", 6), strrep("-", 5),
    strrep("-", 6), strrep("-", 11), strrep("-", 11), strrep("-", 11), strrep("-", 11)))

  for(i in utils::head(ord, 10)) {
    emit(sprintf("  %-6s %5d %-6s %11.3g %11.3g %11.3g %11.3g", genes[i], n_pep[i],
      paste0(substr(pa$type[i], 1, 4), if(!identical(pa$type[i], pb$type[i])) "!" else ""),
      d_p[i], d_eff[i], abs(pa$df[i] - pb$df[i]), abs(pa$sigma[i] - pb$sigma[i])))
  }

  by_route <- split(seq_along(genes), pa$type)

  emit("")
  emit("  By route, over all genes:")
  emit("")
  emit(sprintf("  %-8s %6s %13s %13s", "route", "genes", "worst d(pval)",
    "worst d(L'b)"))
  emit(sprintf("  %-8s %6s %13s %13s", strrep("-", 8), strrep("-", 6), strrep("-", 13),
    strrep("-", 13)))

  for(ty in names(by_route)) {
    ii <- by_route[[ty]]
    emit(sprintf("  %-8s %6d %13.3g %13.3g", ty, length(ii),
      max(d_p[ii], na.rm=TRUE), max(d_eff[ii], na.rm=TRUE)))
  }

  emit("")
  emit("  Route agrees gene for gene between the two runs: ",
    if(identical(pa$type, pb$type)) "yes" else "NO")
}

###############################################################################
## step 4: is it a stopping rule or an estimate? The fits are repeated here directly, in
##   lme4 and MASS, so the tolerances msqrob2 leaves at their defaults can be tightened.
##
##   Neither iteration is non-invariant in exact arithmetic. The robust loop reweights on
##   the residuals, which depend on the column space of the design and not on its
##   parameterization, and each reweighted fit projects onto that same space, so the whole
##   trajectory is the same trajectory. What is not invariant is the objective's value: the
##   REML criterion carries a log|X'V^-1 X| term, which under X -> XR gains 2 log|det R|, a
##   constant in the variance parameters. The optimum is therefore unmoved but the value at
##   it is shifted, and lme4's default optimizer stops on an absolute tolerance,
##   nloptwrap's ftol_abs, so the two parameterizations stop at different points on the
##   same surface. If that is what this is, tightening the tolerance shrinks the deviation
##   toward floating point; if the deviation is indifferent to the tolerance, it is not.

long_gene <- function(g) {
  ii <- which(feats$gene %in% g)
  return(data.frame(
    y=as.vector(exprs[ii, , drop=FALSE]),
    pep=rep(feats$pep[ii], times=nobs),
    obs=rep(samps$obs, each=length(ii)),
    sex=rep(samps$sex, each=length(ii)),
    stringsAsFactors=FALSE
  ))
}

## msqrob2's robust loop, verbatim in structure from msqrob2:::.noridge_msqrobLmer(), with
##   the lmer control and the loop tolerance exposed:

fit_one <- function(dat, frm, robust, ctrl, tol, maxit=100) {
  m <- try(lme4::lmer(frm, dat, control=ctrl), silent=TRUE)
  if(inherits(m, "try-error")) return(NULL)
  m@frame$`(weights)` <- rep(1, dim(m@frame)[1])
  if(robust) {
    sse_old <- m@devcomp$cmp["pwrss"]
    it <- 0
    while(it < maxit) {
      it <- it + 1
      res <- stats::resid(m)
      m@frame$`(weights)` <- MASS::psi.huber(res / stats::mad(res, 0))
      m <- lme4::refit(m)
      sse <- m@devcomp$cmp["pwrss"]
      if(abs(sse_old - sse) / sse_old <= tol) break
      sse_old <- sse
    }
  }
  return(m)
}

eff_of <- function(m, L) {
  b <- lme4::fixef(m)
  if(!all(names(L) %in% names(b))) return(NA_real_)
  return(sum(L * b[names(L)]))
}

frms <- list(
  A=stats::as.formula("y ~ sex + (1|pep) + (1|obs)"),
  B=stats::as.formula("y ~ 0 + sex + (1|pep) + (1|obs)")
)

ctrl_loose <- lme4::lmerControl(calc.derivs=FALSE)
ctrl_tight <- lme4::lmerControl(calc.derivs=FALSE,
  optCtrl=list(ftol_abs=1e-14, ftol_rel=1e-14, xtol_abs=1e-14, xtol_rel=1e-14,
    maxeval=100000))

settings <- list(
  list(key="lmer_default", robust=FALSE, ctrl=ctrl_loose, tol=1e-6,
    lab="plain lmer, lme4 defaults"),
  list(key="lmer_tight", robust=FALSE, ctrl=ctrl_tight, tol=1e-6,
    lab="plain lmer, optimizer tightened"),
  list(key="rob_default", robust=TRUE, ctrl=ctrl_loose, tol=1e-6,
    lab="robust, both as msqrob2 has them"),
  list(key="rob_tight_opt", robust=TRUE, ctrl=ctrl_tight, tol=1e-6,
    lab="robust, optimizer tightened only"),
  list(key="rob_tight_both", robust=TRUE, ctrl=ctrl_tight, tol=1e-14,
    lab="robust, optimizer and loop tightened")
)

genes_mixed <- genes[n_pep >= 2]
genes_thin <- genes[n_pep < 2]

progress("refitting ", length(genes_mixed), " genes under ", length(settings),
  " tolerance settings")

devs <- lapply(settings, function(s) numeric(0))
names(devs) <- vapply(settings, function(s) s$key, character(1))

for(k in seq_along(genes_mixed)) {

  dat <- long_gene(genes_mixed[k])

  for(s in settings) {
    ma <- fit_one(dat, frms$A, s$robust, s$ctrl, s$tol)
    mb <- fit_one(dat, frms$B, s$robust, s$ctrl, s$tol)
    if(is.null(ma) || is.null(mb)) next
    devs[[s$key]] <- c(devs[[s$key]],
      abs(eff_of(ma, runs$A$L) - eff_of(mb, runs$B$L)))
  }

  if(k %% 5 == 0) progress("  ", k, "/", length(genes_mixed), " genes refitted")
}

emit("")
emit("Step 4. The mixed model refitted directly, ", length(genes_mixed), " genes with two")
emit("  or more peptides, comparing L'beta between the two parameterizations. Worst")
emit("  absolute deviation over those genes:")
emit("")
emit(sprintf("  %-38s %12s", "fit", "worst d(L'b)"))
emit(sprintf("  %-38s %12s", strrep("-", 38), strrep("-", 12)))

for(s in settings) {
  v <- devs[[s$key]]
  emit(sprintf("  %-38s %12.3g", s$lab, if(length(v)) max(v) else NA_real_))
}

## and the gene level route, MASS::rlm() on the aggregated values, whose default accuracy
##   is 1e-4, looser than anything above:

emit("")

if(length(genes_thin)) {

  progress("refitting ", length(genes_thin), " single-peptide genes at the gene level")

  agg <- t(vapply(genes_thin, function(g) {
    ii <- which(feats$gene %in% g)
    return(colMeans(exprs[ii, , drop=FALSE], na.rm=TRUE))
  }, numeric(nobs)))

  X <- list(
    A=stats::model.matrix(~sex, samps),
    B=stats::model.matrix(~0 + sex, samps)
  )

  rlm_eff <- function(y, acc) {
    out <- rep(NA_real_, 2)
    for(j in 1:2) {
      nm <- names(runs)[j]
      m <- try(MASS::rlm(X[[nm]], y, method="M", maxit=100, acc=acc), silent=TRUE)
      if(inherits(m, "try-error")) return(NA_real_)
      b <- m$coefficients
      L <- runs[[nm]]$L
      if(!all(names(L) %in% names(b))) return(NA_real_)
      out[j] <- sum(L * b[names(L)])
    }
    return(abs(out[1] - out[2]))
  }

  d_def <- apply(agg, 1, rlm_eff, acc=1e-4)
  d_tig <- apply(agg, 1, rlm_eff, acc=1e-14)

  emit("  The gene level route refitted directly, ", length(genes_thin),
    " single-peptide genes, MASS::rlm()")
  emit("  as msqrob2 calls it. Worst absolute deviation in L'beta:")
  emit("")
  emit(sprintf("  %-38s %12s", "fit", "worst d(L'b)"))
  emit(sprintf("  %-38s %12s", strrep("-", 38), strrep("-", 12)))
  emit(sprintf("  %-38s %12.3g", "rlm, acc as msqrob2 has it (1e-4)",
    max(d_def, na.rm=TRUE)))
  emit(sprintf("  %-38s %12.3g", "rlm, acc tightened to 1e-14", max(d_tig, na.rm=TRUE)))

} else {
  emit("  No single-peptide genes in the fixture, so the gene level route is not")
  emit("  exercised here.")
}

###############################################################################
## step 5: the loop traced, which is what settles it. Two things are recorded per robust
##   iteration for the genes that move most and one that does not: the variance parameters
##   lme4 arrived at under each parameterization, and a cross fit, the second design fitted
##   with the first design's weights. If the reweighted fit itself were parameterization
##   dependent the cross fit would separate from the first design's own; if instead the two
##   parameterizations are being handed different variance parameters, the cross fit stays
##   with it and theta is where they part.

trace_gene <- function(g, maxit=6) {

  dat <- long_gene(g)

  ma <- lme4::lmer(frms$A, dat, control=ctrl_loose)
  mb <- lme4::lmer(frms$B, dat, control=ctrl_loose)
  ma@frame$`(weights)` <- rep(1, dim(ma@frame)[1])
  mb@frame$`(weights)` <- rep(1, dim(mb@frame)[1])

  rec <- NULL

  for(it in seq_len(maxit)) {

    ra <- stats::resid(ma)
    rb <- stats::resid(mb)
    wa <- MASS::psi.huber(ra / stats::mad(ra, 0))
    wb <- MASS::psi.huber(rb / stats::mad(rb, 0))

    ma@frame$`(weights)` <- wa
    ma <- lme4::refit(ma)
    mb@frame$`(weights)` <- wb
    mb <- lme4::refit(mb)

    ## design B, weights from A: one reweighted fit, both parameterizations, same weights:

    cross <- lme4::lmer(frms$B, dat, control=ctrl_loose)
    cross@frame$`(weights)` <- wa
    cross <- lme4::refit(cross)

    th_a <- lme4::getME(ma, "theta")
    th_b <- lme4::getME(mb, "theta")

    rec <- rbind(rec, data.frame(it=it, d_eff=eff_of(ma, runs$A$L) - eff_of(mb, runs$B$L),
      d_cross=eff_of(ma, runs$A$L) - eff_of(cross, runs$B$L),
      th_a=max(th_a), th_b=max(th_b), th_a_min=min(th_a), th_b_min=min(th_b)))
  }

  return(rec)
}

worst <- genes_mixed[order(vapply(genes_mixed, function(g) {
  i <- match(g, genes)
  return(abs(pa$eff[i] - pb$eff[i]))
}, numeric(1)), decreasing=TRUE)]

show <- c(utils::head(worst, 2), utils::tail(worst, 1))

emit("")
emit("Step 5. The robust loop traced for the two genes that move most and the one that")
emit("  moves least. d(eff) is between the two parameterizations; d(cross) is design B")
emit("  refitted with design A's own weights, against design A; theta_max and theta_min")
emit("  are the largest and smallest variance parameter lme4 arrived at:")

for(g in show) {

  progress("tracing gene ", g)
  tr <- trace_gene(g)

  emit("")
  emit("  gene ", g, ", ", n_pep[match(g, genes)], " peptides:")
  emit("")
  emit(sprintf("  %3s %11s %11s %11s %11s %11s %11s", "it", "d(eff)", "d(cross)",
    "th_max A", "th_max B", "th_min A", "th_min B"))
  emit(sprintf("  %3s %11s %11s %11s %11s %11s %11s", strrep("-", 3), strrep("-", 11),
    strrep("-", 11), strrep("-", 11), strrep("-", 11), strrep("-", 11), strrep("-", 11)))

  for(i in seq_len(nrow(tr))) {
    emit(sprintf("  %3d %11.2e %11.2e %11.4f %11.4f %11.2e %11.2e", tr$it[i],
      tr$d_eff[i], tr$d_cross[i], tr$th_a[i], tr$th_b[i], tr$th_a_min[i],
      tr$th_b_min[i]))
  }
}

emit("")
emit(strrep("=", 78))

writeLines(rpt, con=out_file)

progress("done; report written to ", out_file)

quit(save="no", status=0)
