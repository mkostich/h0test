## Tests what an exported h0testr::test_*() engine does when it is called directly
##   rather than through h0testr::test_h0(), on a state whose factor covariates arrive as
##   character columns. That is the ordinary shape of a state read back from the .tsv
##   files h0testr writes, since utils::read.table() returns character.
##   h0testr::f.design_X() hands state$samples to stats::model.matrix() as they are, and
##   a character column there is levelled by sorting, so config$reference_levels reached
##   the design only through h0testr::init_state() or h0testr::test_h0(): a direct caller
##   got coefficients named for the alphabetically first level instead of the declared
##   reference, with the sign of every effect flipped and nothing said about it.
##   The central assertion is therefore an invariance: for each engine, the result of a
##   call on character covariates equals the result of the same call on the factor
##   columns init_state() resolves, and the reported effect has the sign the declared
##   reference level implies rather than the one sorting implies.
##   Also checks that every exported engine validates its config and state, and that a
##   config naming one column as both the feature and the gene id is refused by the
##   engine's own message, which names the method to use instead, rather than by the
##   state check that also has grounds to refuse it.

usage <- function(msg=NULL) {
  if(!is.null(msg)) cat("ERROR:", msg, "\n\n", file=stderr())
  cat(
    "Test the exported h0testr::test_*() engines under direct calls: that",
    "config$reference_levels reaches the design when the covariates arrive as",
    "character columns rather than as the factors h0testr::init_state() resolves, so",
    "that a direct call agrees with h0testr::test_h0() about which level is the reference",
    "and about the sign of the effect; that each engine validates its config and state",
    "instead of failing further in; and that an already aggregated config is refused by",
    "the engine's own message rather than by the state check.",
    "",
    "Usage: Rscript test_direct_call.R <r_dir> [--engines=<list>]",
    "",
    "Required positional arguments:",
    "  <r_dir>  Path to the h0testr package R/ source directory; all .R files",
    "             there are sourced (the installed package is not used).",
    "",
    "Optional named arguments:",
    "  --engines=<list>  Comma separated subset of lm,trend,voom,deqms,msqrob,proda,",
    "                      prolfqua to exercise; default is all of them.",
    "",
    "Output: one 'PASS: <what>' or 'FAIL: <what>' line per assertion on stdout,",
    "  section headers with elapsed seconds, and a final count of passes and",
    "  failures. Exit code 0 if every assertion passed, 1 if any failed, 2 on a",
    "  usage error.",
    "",
    "Examples:",
    "  Rscript test_direct_call.R C:/path/to/h0testr/R",
    "  Rscript test_direct_call.R C:/path/to/h0testr/R --engines=lm,trend",
    "  Rscript test_direct_call.R ../../h0test/h0testr/R --engines=msqrob",
    "",
    sep="\n", file=stderr()
  )
  quit(save="no", status=2)
}

args <- commandArgs(trailingOnly=TRUE)
if(length(args) < 1 || length(args) > 2) usage("wrong number of arguments")

r_dir <- args[1]
if(!dir.exists(r_dir)) usage(paste("r_dir not a directory:", r_dir))

engines_all <- c("lm", "trend", "voom", "deqms", "msqrob", "proda", "prolfqua")
engines <- engines_all

if(length(args) %in% 2) {
  if(!grepl("^--engines=", args[2])) usage(paste("unrecognized argument:", args[2]))
  engines <- trimws(strsplit(sub("^--engines=", "", args[2]), ",")[[1]])
  bad <- setdiff(engines, engines_all)
  if(length(bad)) usage(paste("unknown engine(s):", paste(bad, collapse=", ")))
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
mark <- function() length(readLines(log_file))
logged <- function(pat, since=0) {
  txt <- readLines(log_file)
  if(since >= length(txt)) return(FALSE)
  any(grepl(pat, txt[(since + 1):length(txt)], fixed=TRUE))
}

###############################################################################
## shared data: 20 genes of 1 to 4 peptides, 12 observations in two groups of six,
##   with a crossed second covariate. The declared reference level of grp is "trt",
##   which is not the alphabetically first level, so sorting and the declaration
##   disagree and the coefficient they name is not the same one. The peptide count
##   varies by gene because DEqMS fits its variance prior against that count and
##   test_deqms() refuses a state where every gene has the same one:

set.seed(101)

nobs <- 12
ngene <- 20
npep_max <- 4

samps <- data.frame(
  obs=paste0("o", sprintf("%02d", 1:nobs)),
  grp=rep(c("ctl", "trt"), each=nobs / 2),
  sex=rep(c("F", "M"), nobs / 2),
  stringsAsFactors=FALSE
)

## the first ten genes are up in the treated group, so an effect with a known sign
##   exists: with reference "trt" the tested coefficient is grpctl and the effect is
##   negative; with the reference sorting would pick it is grptrt and positive:

gene_eff <- matrix(0, nrow=ngene, ncol=nobs)
gene_eff[1:10, samps$grp %in% "trt"] <- 3

npep <- ((seq_len(ngene) - 1) %% npep_max) + 1         ## 1..4 peptides per gene
gene <- rep(paste0("g", sprintf("%02d", seq_len(ngene))), times=npep)

exprs <- matrix(stats::rnorm(length(gene) * nobs, 20, 1.5), nrow=length(gene))
exprs <- exprs + gene_eff[match(gene, unique(gene)), , drop=FALSE]
rownames(exprs) <- paste0(gene, ".p", unlist(lapply(npep, seq_len)))
colnames(exprs) <- samps$obs
feats <- data.frame(pep=rownames(exprs), gene=gene, stringsAsFactors=FALSE)

cfg0 <- new_config()
cfg0$obs_id_col <- cfg0$obs_col <- "obs"
cfg0$sample_id_col <- "obs"
cfg0$feat_id_col <- cfg0$feat_col <- "pep"
cfg0$gene_id_col <- "gene"
cfg0$reference_levels <- c(grp="trt", sex="F")
cfg0$estimability <- "test"
cfg0$df_resid_min <- 2
cfg0$save_state <- FALSE
cfg0$permute_var <- ""
cfg0$is_log_transformed <- TRUE     ## simulated on a log scale; normalize() not called
cfg0$impute_method <- "none"
cfg0$n_features_min <- 5           ## the default of 1000 would filter out everything
cfg0$log_file <- log_file
cfg0$frm <- ~grp + sex
cfg0$test_term <- "grp"

out <- suppressMessages(init_state(list(expression=exprs, features=feats,
  samples=samps), cfg0, minimal=TRUE))
state_fac <- out$state
config0 <- out$config

## init_state() resolved the declared reference level into the factor columns; the
##   direct-call state is the same data with those columns back as character, which is
##   what read.table() hands back and what an assembled state usually carries:

state_chr <- state_fac
for(nom in c("grp", "sex")) {
  state_chr$samples[[nom]] <- as.character(state_fac$samples[[nom]])
}

report(identical(levels(state_fac$samples$grp), c("trt", "ctl")),
  "init_state() puts the declared reference level first")
report(is.character(state_chr$samples$grp),
  "the direct-call state carries grp as character")
report(identical(sort(unique(state_chr$samples$grp)), c("ctl", "trt")),
  "sorting that column would pick the other level as the reference")

## a gene level state for the engines that report one row per gene:

out_g <- suppressMessages(combine_features(state_fac, config0, method="medianPolish",
  rescale=FALSE))
state_gene_fac <- out_g$state
config_gene <- out_g$config

state_gene_chr <- state_gene_fac
for(nom in c("grp", "sex")) {
  state_gene_chr$samples[[nom]] <- as.character(state_gene_fac$samples[[nom]])
}

###############################################################################
## the invariance each engine is checked against: the same call on character and on
##   factor covariates, compared on the columns that would move if the two designs
##   disagreed about the reference level. Sign is checked as well as equality, since
##   two runs that were both wrong in the same way would be equal:

chk_engine <- function(nom, res_chr, res_fac, pcol, fcol, ncol_expect=NULL,
    sign_check=TRUE) {

  h_chr <- res_chr$hits
  h_fac <- res_fac$hits

  report(isTRUE(all.equal(h_chr[[pcol]], h_fac[[pcol]])),
    paste(nom, ": a direct call on character covariates gives the p-values of the",
      "factor call"))

  if(!is.null(fcol) && fcol %in% names(h_chr)) {

    report(isTRUE(all.equal(h_chr[[fcol]], h_fac[[fcol]])),
      paste(nom, ": and the same effect estimates"))

    ## the ten genes up in the treated group, under reference level "trt": the
    ##   coefficient is ctl - trt, so the effect is negative. Sorting would name
    ##   grptrt instead and report it positive. Skipped for test_voom(), whose
    ##   response is the log counts per million limma::voom() derives rather than the
    ##   values handed to it: the fixture here is on a log scale, so voom's library
    ##   size scaling divides the simulated effect out and the sign of what it reports
    ##   says nothing about the reference level. The invariance above still does:

    if(sign_check) {
      ## whichever column this engine keys its table by: the gene when it reports one
      ##   row per gene, the feature otherwise, and 'name' from proDA::test_diff():

      ids <- NULL
      for(idc in c("gene", "feature", "name", "pep")) {
        if(is.null(ids) && idc %in% names(h_chr)) ids <- h_chr[[idc]]
      }
      if(is.null(ids)) ids <- rownames(h_chr)
      up <- grepl("^g0[1-9]$|^g10$", sub("[.]p[0-9]+$", "", as.character(ids)))
      report(sum(up) > 0 && all(h_chr[[fcol]][up] < 0),
        paste(nom, ": the effect has the sign the declared reference level implies"))
    }
  }

  if(!is.null(ncol_expect)) {
    report(ncol_expect %in% names(h_chr),
      paste(nom, ": the reported coefficient is named for the declared reference"))
    report(!any(grepl("^grptrt$", names(h_chr))),
      paste(nom, ": and not for the level sorting would have picked"))
  }
}

###############################################################################
if("lm" %in% engines) {
  section("test_lm(): direct call")

  cfg <- config0
  r_chr <- try(suppressMessages(test_lm(state_chr, cfg)), silent=TRUE)
  r_fac <- try(suppressMessages(test_lm(state_fac, cfg)), silent=TRUE)

  report(!inherits(r_chr, "try-error"), "test_lm: character covariates accepted")

  if(!inherits(r_chr, "try-error") && !inherits(r_fac, "try-error")) {
    chk_engine("test_lm", r_chr, r_fac, "pval", "grpctl", ncol_expect="grpctl")
  }
}

###############################################################################
if("trend" %in% engines) {
  section("test_trend(): direct call")

  r_chr <- try(suppressMessages(test_trend(state_chr, config0)), silent=TRUE)
  r_fac <- try(suppressMessages(test_trend(state_fac, config0)), silent=TRUE)

  report(!inherits(r_chr, "try-error"), "test_trend: character covariates accepted")

  if(!inherits(r_chr, "try-error") && !inherits(r_fac, "try-error")) {
    chk_engine("test_trend", r_chr, r_fac, "P.Value", "logFC")
    report(identical(colnames(r_chr$fit$coefficients),
      colnames(r_fac$fit$coefficients)),
      "test_trend: the two fits name the same coefficients")
  }
}

###############################################################################
if("voom" %in% engines) {
  section("test_voom(): direct call")

  r_chr <- try(suppressMessages(test_voom(state_chr, config0)), silent=TRUE)
  r_fac <- try(suppressMessages(test_voom(state_fac, config0)), silent=TRUE)

  report(!inherits(r_chr, "try-error"), "test_voom: character covariates accepted")

  if(!inherits(r_chr, "try-error") && !inherits(r_fac, "try-error")) {
    chk_engine("test_voom", r_chr, r_fac, "P.Value", "logFC", sign_check=FALSE)
  }
}

###############################################################################
if("deqms" %in% engines) {
  section("test_deqms(): direct call")

  r_chr <- try(suppressMessages(test_deqms(state_chr, config0)), silent=TRUE)
  r_fac <- try(suppressMessages(test_deqms(state_fac, config0)), silent=TRUE)

  report(!inherits(r_chr, "try-error"), "test_deqms: character covariates accepted")

  if(!inherits(r_chr, "try-error") && !inherits(r_fac, "try-error")) {
    chk_engine("test_deqms", r_chr, r_fac, "sca.P.Value", "logFC")
  }
}

###############################################################################
if("msqrob" %in% engines) {
  section("test_msqrob(): direct call")

  r_chr <- try(suppressMessages(test_msqrob(state_chr, config0)), silent=TRUE)
  r_fac <- try(suppressMessages(test_msqrob(state_fac, config0)), silent=TRUE)

  report(!inherits(r_chr, "try-error"),
    "test_msqrob: character covariates accepted, rather than erroring on a design mismatch")

  if(!inherits(r_chr, "try-error") && !inherits(r_fac, "try-error")) {
    chk_engine("test_msqrob", r_chr, r_fac, "pval", "logFC")
  }
}

###############################################################################
if("proda" %in% engines) {
  section("test_proda(): direct call")

  r_chr <- try(suppressMessages(test_proda(state_chr, config0)), silent=TRUE)
  r_fac <- try(suppressMessages(test_proda(state_fac, config0)), silent=TRUE)

  report(!inherits(r_chr, "try-error"), "test_proda: character covariates accepted")

  if(!inherits(r_chr, "try-error") && !inherits(r_fac, "try-error")) {
    chk_engine("test_proda", r_chr, r_fac, "pval", "diff")
  }
}

###############################################################################
if("prolfqua" %in% engines) {
  section("test_prolfqua(): direct call")

  r_chr <- try(suppressMessages(test_prolfqua(state_chr, config0)), silent=TRUE)
  r_fac <- try(suppressMessages(test_prolfqua(state_fac, config0)), silent=TRUE)

  report(!inherits(r_chr, "try-error"), "test_prolfqua: character covariates accepted")

  if(!inherits(r_chr, "try-error") && !inherits(r_fac, "try-error")) {
    chk_engine("test_prolfqua", r_chr, r_fac, "pval", NULL)
  }
}

###############################################################################
section("every exported engine validates its config and state")

## a config the workflow's own check refuses, and a state whose expression columns do
##   not match its samples. An engine that does not check either one fails somewhere
##   further in, on a message about whatever it happened to reach first:

cfg_bad <- config0
cfg_bad$normalization_method <- "no_such_method"

state_bad <- state_fac
colnames(state_bad$expression) <- paste0("x", seq_len(ncol(state_bad$expression)))

engine_fns <- list(
  test_lm=function(st, cf) test_lm(st, cf),
  test_trend=function(st, cf) test_trend(st, cf),
  test_voom=function(st, cf) test_voom(st, cf),
  test_deqms=function(st, cf) test_deqms(st, cf),
  test_msqrob=function(st, cf) test_msqrob(st, cf),
  test_proda=function(st, cf) test_proda(st, cf),
  test_prolfqua=function(st, cf) test_prolfqua(st, cf)
)

for(nom in names(engine_fns)) {
  fn <- engine_fns[[nom]]
  report(threw(suppressMessages(fn(state_fac, cfg_bad))),
    paste(nom, ": a config check_config() refuses is refused here too"))
  report(threw(suppressMessages(fn(state_bad, config0))),
    paste(nom, ": a state whose columns do not match its samples is refused"))
}

###############################################################################
section("an already aggregated config is refused by the engine's own message")

## config$feat_id_col and config$gene_id_col naming one column is what
##   combine_features() leaves behind, and the mixed and aggregate paths cannot use
##   it. f.check_state() also has grounds to refuse this pairing against a feature
##   level state, so the check that reads config alone runs first: of the two
##   refusals it is the one that names the method to use instead:

cfg_agg <- config0
cfg_agg$feat_id_col <- cfg_agg$gene_id_col <- "gene"

m0 <- mark()
report(threw(suppressMessages(test_msqrob(state_fac, cfg_agg, aggregate=TRUE))),
  "test_msqrob: one column named as both the feature and the gene id is an error")
report(logged("the aggregate path needs feature level input", m0),
  "test_msqrob: and the message says what is wrong")
report(logged("use that instead", m0),
  "test_msqrob: and names the method to use instead")

m0 <- mark()
report(threw(suppressMessages(test_prolfqua(state_fac, cfg_agg, mixed=TRUE))),
  "test_prolfqua: one column named as both the feature and the gene id is an error")
report(logged("the mixed path needs feature level input", m0),
  "test_prolfqua: and the message says what is wrong")
report(logged("use that instead", m0),
  "test_prolfqua: and names the method to use instead")

###############################################################################
cat("\n## passes:", n_pass, "; failures:", n_fail, "; elapsed:",
  round(as.numeric(difftime(Sys.time(), t0, units="secs"))), "s\n")

quit(save="no", status=if(n_fail > 0) 1 else 0)
