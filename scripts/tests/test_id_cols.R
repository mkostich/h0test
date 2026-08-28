## Tests which column of state$features identifies the rows a test method returns.
##   config$feat_col always matches rownames(state$expression), so a row-wise engine is
##   keyed by it whether the state holds precursors or has been through
##   combine_features(); only deqms and msqrob, which aggregate for themselves, are
##   keyed by config$gene_id_col. f.format_lm() used to read config$gene_id_col for an
##   engine that is row-wise, so test(method="lm") failed outright on precursor-level
##   data. Runs all seven engines at both levels, and checks not only that they run but
##   that each result row is joined to the metadata of the feature it actually describes.

usage <- function(msg=NULL) {
  if(!is.null(msg)) cat("ERROR:", msg, "\n\n", file=stderr())
  cat(
    "Test feature/gene id handling in h0testr: f.gene_ids(), f.gene_features(),",
    "f.gene_counts(), f.test_id_col(), the gene-level form of a feature metadata",
    "table, and all seven test_method engines run through test() twice, once on",
    "precursor-level data with config$feat_id_col != config$gene_id_col and once on",
    "the same data after combine_features(). Also checks that test_deqms() refuses a",
    "run in which every gene has the same number of features, and that f.tune2()",
    "records such a combination as untested instead of stopping the sweep.",
    "",
    "Usage: Rscript test_id_cols.R <r_dir>",
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
    "  3 a required test package is not installed.",
    "",
    "Examples:",
    "  Rscript test_id_cols.R C:/path/to/h0test/h0testr/R",
    "  Rscript test_id_cols.R ../../h0test/h0testr/R",
    "  Rscript test_id_cols.R C:/path/to/h0testr/R > test_id_cols.out 2>&1",
    sep="\n", file=stderr()
  )
  quit(status=2)
}

args <- commandArgs(trailingOnly=TRUE)
if(length(args) != 1) usage("wrong number of arguments")

r_dir <- args[1]
if(!dir.exists(r_dir)) usage(paste("r_dir not a directory:", r_dir))

for(pkg in c("limma", "DEqMS", "msqrob2", "proDA", "prolfqua", "MsCoreUtils")) {
  if(!requireNamespace(pkg, quietly=TRUE)) {
    cat("ERROR:", pkg, "not installed\n", file=stderr())
    quit(status=3)
  }
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

## an expected error: TRUE iff the call fails and the log picks up every pattern:

errs_with <- function(expr, ...) {
  n0 <- log_len()
  out <- try(suppressMessages(expr), silent=TRUE)
  if(!inherits(out, "try-error")) return(FALSE)
  if(length(list(...)) %in% 0) return(TRUE)
  return(log_has(n0, ...))
}

runs_with <- function(expr, ...) {
  n0 <- log_len()
  out <- try(suppressMessages(expr), silent=TRUE)
  if(inherits(out, "try-error")) return(FALSE)
  if(length(list(...)) %in% 0) return(TRUE)
  return(log_has(n0, ...))
}

close_enough <- function(a, b, tol=1e-8) {
  if(length(a) != length(b)) return(FALSE)
  i <- !is.na(a) & !is.na(b)
  if(!any(i)) return(FALSE)
  if(any(is.na(a) != is.na(b))) return(FALSE)
  return(max(abs(a[i] - b[i])) < tol)
}

###############################################################################
## Twenty observations, fully crossed in sex and batch; sex is the test term and
##   resolves to one coefficient, so every engine including deqms can run it. Thirty
##   genes with one to four peptides each, so the counts DEqMS fits its variance prior
##   against actually vary, and one further gene whose gene id is blank, so that the
##   unknown_* naming is exercised end to end. Per-gene variances vary, giving the
##   moderating engines a prior to estimate. Two marker columns: one that varies within
##   a gene and so is a property of the peptide, one that does not and so is a property
##   of the gene; the join checks below read them to confirm a result row was matched to
##   the right feature rather than merely to some feature.

set.seed(202)

nobs <- 20
ngene <- 30

samps <- expand.grid(rep=1:5, batch=c("b1", "b2"), sex=c("F", "M"),
  stringsAsFactors=FALSE)
samps$obs <- paste0("o", sprintf("%02d", 1:nobs))
samps <- samps[, c("obs", "sex", "batch")]

npep <- rep(c(1, 2, 3, 4), length.out=ngene)
gene <- rep(paste0("g", sprintf("%02d", 1:ngene)), times=npep)
pep <- paste0(gene, ".p", unlist(lapply(npep, seq_len)))

## two peptides with no gene assignment, which f.gene_ids() turns into one unknown_*
##   gene each rather than pooling them:

gene <- c(gene, "", "")
pep <- c(pep, "orphan.p1", "orphan.p2")
nfeat <- length(pep)

gene_eff <- matrix(0, nrow=nfeat, ncol=nobs)
gene_eff[gene %in% paste0("g", sprintf("%02d", 1:8)), samps$sex %in% "M"] <- 3

sds <- rep(sqrt(1 / rgamma(nfeat, shape=3, rate=3)), times=1)
exprs <- matrix(rnorm(nfeat * nobs, 20, rep(sds, times=nobs)), nrow=nfeat)
exprs <- exprs + gene_eff
rownames(exprs) <- pep
colnames(exprs) <- samps$obs

## the gene each peptide ends up in, which for the two orphans is the unknown_* gene
##   f.gene_ids() gives it; the gene-level marker is built from this so that a join
##   check covers the orphans as well as the ordinary genes:

gene_res <- ifelse(gene %in% "", paste0("unknown_", pep), gene)

feats <- data.frame(pep=pep, gene=gene,
  pep_tag=paste0("T_", pep),                       ## varies within a gene
  gene_tag=paste0("G_", gene_res),                 ## constant within a gene
  stringsAsFactors=FALSE)

cfg0 <- new_config()
cfg0$obs_id_col <- cfg0$obs_col <- "obs"
cfg0$sample_id_col <- "obs"
cfg0$feat_id_col <- cfg0$feat_col <- "pep"
cfg0$gene_id_col <- "gene"
cfg0$frm <- ~sex + batch
cfg0$test_term <- "sex"
cfg0$reference_levels <- c(sex="F", batch="b1")
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
state_pep <- out$state
config_pep <- out$config

methods <- c("lm", "trend", "deqms", "msqrob", "proda", "prolfqua", "voom")

## the gene ids f.gene_ids() assigns, which the two orphan peptides turn into their own
##   genes; every gene-level expectation below is built from these:

genes_all <- f.gene_ids(state_pep$features, config_pep, "suite")

###############################################################################
section("f.gene_ids()")

report(identical(genes_all[1:3], c("g01", "g02", "g02")),
  "f.gene_ids() returns the gene column where it is populated")

report(identical(genes_all[nfeat - 1:0], c("unknown_orphan.p1", "unknown_orphan.p2")),
  "f.gene_ids() gives a feature with a blank gene id a gene of its own")

fx <- state_pep$features
fx$gene[1] <- NA
report(identical(f.gene_ids(fx, config_pep, "suite")[1], "unknown_g01.p1"),
  "f.gene_ids() treats a missing gene id the same as a blank one")

report(errs_with(f.gene_ids(fx[, c("pep", "pep_tag")], config_pep, "suite"),
  "!(config$gene_id_col %in% names(feats))"),
  "f.gene_ids() errors when config$gene_id_col is not a column of the table")

###############################################################################
section("f.gene_features()")

gf <- f.gene_features(state_pep$features, config_pep, "suite")

report(nrow(gf$features) %in% length(unique(genes_all)),
  "f.gene_features() returns one row per gene")

report(identical(gf$features$gene, unique(genes_all)),
  "f.gene_features() keeps the genes in order of first appearance")

report(identical(gf$genes, genes_all),
  "f.gene_features() reports the per-input-row gene ids it grouped by")

report("gene_tag" %in% names(gf$features),
  "f.gene_features() keeps a column that is constant within every gene")

report(!("pep_tag" %in% names(gf$features)),
  "f.gene_features() drops a column that varies within a gene")

report(!("pep" %in% names(gf$features)),
  "f.gene_features() drops config$feat_id_col, which always varies within a gene")

n0 <- log_len()
invisible(suppressMessages(f.gene_features(state_pep$features, config_pep, "suite")))
report(log_has(n0, "feature metadata", "pep_tag"),
  "f.gene_features() names the columns it dropped")

report(identical(gf$features$gene_tag, paste0("G_", gf$features$gene)),
  "the kept column holds the value shared by the features of its gene")

report(identical(as.integer(gf$features[[config_pep$n_feats_col]]),
  as.integer(table(genes_all)[gf$features$gene])),
  "f.gene_features() adds the number of features aggregated into each gene")

## an all-NA column agrees with itself within every gene and so is a gene-level column:

fx <- state_pep$features
fx$whatever <- NA
report("whatever" %in% names(f.gene_features(fx, config_pep, "suite")$features),
  "f.gene_features() counts NA as agreeing with NA within a gene")

## a table that has been through here already is one row per gene, so recounting would
##   report 1 for every gene and lose the counts of the features that went into it:

fx <- gf$features
report(identical(f.gene_features(fx, config_pep, "suite")$features[[
  config_pep$n_feats_col]], fx[[config_pep$n_feats_col]]),
  "f.gene_features() preserves an existing feature count instead of recounting")

###############################################################################
section("f.gene_level_method() and f.test_id_col()")

report(identical(f.gene_level_method(methods),
  methods %in% c("deqms", "msqrob")),
  "f.gene_level_method() is TRUE for deqms and msqrob and FALSE for the rest")

report(all(sapply(setdiff(methods, c("deqms", "msqrob")),
  function(m) f.test_id_col(m, config_pep) %in% "pep")),
  "f.test_id_col() keys a row-wise method by config$feat_col")

report(all(sapply(c("deqms", "msqrob"),
  function(m) f.test_id_col(m, config_pep) %in% "gene")),
  "f.test_id_col() keys deqms and msqrob by config$gene_id_col")

## after combine_features() the two are the same column, so every method agrees:

cfg_agg <- config_pep
cfg_agg$feat_col <- cfg_agg$feat_id_col <- cfg_agg$gene_id_col
report(all(sapply(methods, function(m) f.test_id_col(m, cfg_agg) %in% "gene")),
  "f.test_id_col() gives the gene id for every method once the state is aggregated")

###############################################################################
section("f.gene_counts()")

cnt <- f.gene_counts(state_pep, config_pep, "suite")

report(identical(unname(cnt[gf$features$gene]),
  as.integer(table(genes_all)[gf$features$gene])),
  "f.gene_counts() counts the features of each gene when no count column is present")

report(all(c("unknown_orphan.p1", "unknown_orphan.p2") %in% names(cnt)) &&
  all(cnt[c("unknown_orphan.p1", "unknown_orphan.p2")] %in% 1),
  "f.gene_counts() covers the unknown_* genes rather than leaving them unnamed")

state_gf <- list(expression=state_pep$expression, features=gf$features,
  samples=state_pep$samples)
report(identical(unname(f.gene_counts(state_gf, cfg_agg, "suite")[gf$features$gene]),
  as.integer(gf$features[[config_pep$n_feats_col]])),
  "f.gene_counts() reads config$n_feats_col when the table carries it")

report(length(unique(cnt)) > 1,
  "the fixture's features per gene vary, as DEqMS requires")

###############################################################################
section("combine_features()")

for(agg in c("medianPolish", "robustSummary")) {

  out <- suppressMessages(combine_features(state_pep, config_pep, method=agg))
  st <- out$state
  cf <- out$config

  report(identical(as.integer(st$features[[cf$n_feats_col]]),
    as.integer(table(genes_all)[st$features$gene])),
    paste0("combine_features(", agg, ") records the features per gene"))

  report(all(rownames(st$expression) == st$features[[cf$feat_col]]) &&
    identical(cf$feat_col, cf$gene_id_col),
    paste0("combine_features(", agg, ") leaves feat_col matching the row names"))

  report(("gene_tag" %in% names(st$features)) &&
    !("pep_tag" %in% names(st$features)),
    paste0("combine_features(", agg, ") keeps gene-level and drops feature-level cols"))
}

## the aggregated state every engine is also run against below:

out <- suppressMessages(combine_features(state_pep, config_pep,
  method="medianPolish"))
state_agg <- out$state
config_agg <- out$config

report(nrow(state_agg$expression) %in% length(unique(genes_all)),
  "the aggregated state has one row per gene")

###############################################################################
## The bug itself: every engine run through test() on precursor-level data, where
##   config$feat_id_col and config$gene_id_col are different columns. Before the fix
##   test(method="lm") failed here, because f.format_lm() labelled its rows with gene
##   ids while test() checked them against config$feat_col.

section("test() on precursor-level data: all seven engines")

pep_means <- rowMeans(state_pep$expression, na.rm=TRUE)
gene_means <- tapply(pep_means, genes_all, mean, na.rm=TRUE)

res_pep <- list()

for(m in methods) {

  ok <- try(suppressMessages(test(state_pep, config_pep, method=m)), silent=TRUE)
  report(!inherits(ok, "try-error"),
    paste("test() runs method", m, "on precursor-level data"))
  if(inherits(ok, "try-error")) next
  res_pep[[m]] <- ok

  gene_level <- f.gene_level_method(m)
  want <- if(gene_level) unique(genes_all) else state_pep$features$pep

  report(setequal(ok$standard$feature, want) &&
    nrow(ok$standard) %in% length(want),
    paste("method", m, "reports one row per",
      if(gene_level) "gene" else "peptide"))

  ## the join: each result row carries the metadata of the feature it describes, read
  ##   from the marker column that identifies that level uniquely:

  if(gene_level) {
    join_ok <- identical(ok$original$gene_tag, paste0("G_", ok$standard$feature)) &&
      !("pep_tag" %in% names(ok$original))
  } else {
    join_ok <- identical(ok$original$pep_tag, paste0("T_", ok$standard$feature))
  }
  report(join_ok,
    paste("method", m, "joins each result row to the right feature metadata"))

  report(all(!is.na(ok$standard$expr)),
    paste("method", m, "reports an average expression for every row"))
}

## where f.fill_standard() computes the average expression rather than the engine
##   reporting its own, it is the mean over the right level:

if(!is.null(res_pep$lm)) {
  report(close_enough(res_pep$lm$standard$expr,
    unname(pep_means[res_pep$lm$standard$feature])),
    "method lm reports the mean over the peptide's observations")
}

if(!is.null(res_pep$msqrob)) {
  report(close_enough(res_pep$msqrob$standard$expr,
    unname(gene_means[res_pep$msqrob$standard$feature])),
    "method msqrob reports the mean over the peptides of the gene")
}

## and the test actually found the planted effect, so these are results and not just
##   correctly labelled noise. Every gene carrying the effect has to be represented among
##   the hits and the hits have to be dominated by them, rather than consisting of them
##   alone: adj_pval < 0.05 controls a false discovery rate, it does not promise no false
##   discovery, and over the 75 peptides here one null peptide lands just inside the
##   threshold. This assertion read "only" while test_lm() defaulted fdr.method to "BY",
##   whose extra conservatism happened to hold that peptide out; the default is "BH" now,
##   as it is for every other engine, so the assertion is written to what an FDR
##   threshold actually promises. An NA adjusted p-value means the feature was not tested
##   at all, which is not a hit:

if(!is.null(res_pep$lm)) {
  sdf <- res_pep$lm$standard
  hit <- sdf$feature[!is.na(sdf$adj_pval) & sdf$adj_pval < 0.05]
  hit_gene <- genes_all[match(hit, state_pep$features$pep)]
  eff <- paste0("g", sprintf("%02d", 1:8))
  report(length(hit) > 0 && all(eff %in% hit_gene) && mean(hit_gene %in% eff) >= 0.75,
    "method lm on peptides finds the peptides of the genes carrying the effect")
}

###############################################################################
section("test() on aggregated data: all seven engines")

for(m in methods) {

  ok <- try(suppressMessages(test(state_agg, config_agg, method=m)), silent=TRUE)
  report(!inherits(ok, "try-error"),
    paste("test() runs method", m, "on aggregated data"))
  if(inherits(ok, "try-error")) next

  report(setequal(ok$standard$feature, unique(genes_all)) &&
    nrow(ok$standard) %in% length(unique(genes_all)),
    paste("method", m, "reports one row per gene on aggregated data"))

  report(identical(ok$original$gene_tag, paste0("G_", ok$standard$feature)),
    paste("method", m, "joins each aggregated row to the right gene metadata"))
}

## deqms on aggregated data works only because combine_features() recorded the feature
##   counts; recounting the aggregated table would have given every gene 1 and left the
##   engine with no spread to fit its prior against:

report(config_agg$n_feats_col %in% names(state_agg$features) &&
  length(unique(state_agg$features[[config_agg$n_feats_col]])) > 1,
  "the aggregated state still carries the varying feature counts deqms needs")

###############################################################################
## Uniform feature counts: nothing for DEqMS's variance prior to be fitted against.
##   Refused by the engine, and recorded as an untested combination by f.tune2() so
##   that a sweep finishes its matrix.

section("uniform feature counts")

## thirty genes, two peptides each, and no orphans, so every count is 2:

npep2 <- rep(2, ngene)
gene2 <- rep(paste0("h", sprintf("%02d", 1:ngene)), times=npep2)
pep2 <- paste0(gene2, ".p", unlist(lapply(npep2, seq_len)))
exprs2 <- matrix(rnorm(length(pep2) * nobs, 20, 1), nrow=length(pep2))
rownames(exprs2) <- pep2
colnames(exprs2) <- samps$obs
feats2 <- data.frame(pep=pep2, gene=gene2, stringsAsFactors=FALSE)

out <- initialize(list(expression=exprs2, features=feats2, samples=samps), cfg0,
  minimal=TRUE)
state_u <- out$state
config_u <- out$config

report(errs_with(test_deqms(state_u, config_u),
  "every gene has the same number of features"),
  "test_deqms() refuses a run in which every gene has the same feature count")

report(errs_with(test(state_u, config_u, method="deqms"),
  "every gene has the same number of features"),
  "the refusal reaches test() rather than being swallowed")

## f.tune2() checks the same condition and returns the row shape a real result has, so
##   the sweep carries on and tune_check() reads the combination as never tested:

config_u2 <- config_u
config_u2$test_method <- "deqms"
config_u2$normalization_method <- "none"

row_u <- try(suppressMessages(f.tune2(state_u, config_u2)), silent=TRUE)

report(!inherits(row_u, "try-error") && nrow(row_u) %in% 1 &&
  is.na(row_u$nhits) && is.na(row_u$ntests),
  "f.tune2() records uniform feature counts as an untested combination")

report(!inherits(row_u, "try-error") &&
  identical(names(row_u), names(f.tune2_na_row(config_u2))),
  "the row f.tune2() returns has the columns tune_check() expects")

## and any other failure of the test step is caught the same way, so that one engine
##   that cannot fit a combination costs that cell and not the rest of the sweep. The
##   failure used here is a missing config$test_prior_df, which test() needs for
##   test_method "proda" and refuses without: an unpredicted failure of the test step
##   that no guard above sees coming. It used to be a test_method that does not exist,
##   which check_config() now refuses at the first step of the sweep instead; see below:

config_u3 <- config_pep
config_u3$test_method <- "proda"
config_u3$test_prior_df <- NULL
config_u3$normalization_method <- "none"

n0 <- log_len()
row_b <- try(suppressMessages(f.tune2(state_pep, config_u3)), silent=TRUE)

report(!inherits(row_b, "try-error") && nrow(row_b) %in% 1 && is.na(row_b$ntests),
  "f.tune2() records any other failure of the test step as untested too")

report(log_has(n0, "failed on this", "proda"),
  "f.tune2() says in the log which combination it gave up on")

## a config$test_method that is not one of h0testr::test_methods() is a different kind
##   of failure: no other parameter of the sweep can rescue it, and every step calls
##   check_config() first, so it stops the run at the first step rather than costing one
##   cell. Refused at configuration time on purpose; see check_config():

config_u4 <- config_pep
config_u4$test_method <- "no_such_method"
config_u4$normalization_method <- "none"

n0 <- log_len()
row_c <- try(suppressMessages(f.tune2(state_pep, config_u4)), silent=TRUE)

report(inherits(row_c, "try-error"),
  "a test_method that does not exist stops the sweep rather than costing one cell")

report(log_has(n0, "unexpected test_method", "no_such_method"),
  "and check_config() says so before any step of it runs")

###############################################################################

cat("\n## passes:", n_pass, "; failures:", n_fail, "; elapsed:",
  round(as.numeric(difftime(Sys.time(), t0, units="secs"))), "s\n")
cat("## log file:", log_file, "\n")

quit(save="no", status=if(n_fail > 0) 1 else 0)
