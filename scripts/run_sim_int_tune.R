usage <- function(msg=NULL) {
  if(!is.null(msg)) cat("ERROR:", msg, "\n\n", file=stderr())
  cat(
    "Run h0testr::tune() over a simulated sex * age fixture, once unpermuted and once per",
    "permutation, then h0testr::tune_check() over the result files. The configuration",
    "comes from the sim.rds that make_sim_int.R wrote, so it matches the fixture.",
    "",
    "Usage: Rscript run_sim_int_tune.R <r_dir> <dir_in> <out_dir> [--n_perm=<int>]",
    "         [--permute_var=<name>] [--seed=<int>] [--norms=<a,b,c>] [--imputes=<a,b,c>]",
    "         [--tests=<a,b,c>] [--iquants=<a,b>] [--iscales=<a,b>] [--prefix=<str>]",
    "         [--suffix=<str>] [--fdr_cutoff=<num>]",
    "",
    "Required positional arguments:",
    "  <r_dir>    Path to the h0testr package R/ source directory; every .R file there is",
    "               sourced, so the installed package is not used.",
    "  <dir_in>   Fixture directory from make_sim_int.R; must hold expression.tsv,",
    "               features.tsv, samples.tsv and sim.rds.",
    "  <out_dir>  Directory the results are written to; created if absent.",
    "",
    "Optional named arguments:",
    "  --n_perm=<int>        Permuted runs, alongside the one unpermuted run; default 3.",
    "  --permute_var=<name>  Sample variable shuffled in the permuted runs; default sex,",
    "                          the only single variable whose shuffle breaks sex:age.",
    "  --seed=<int>          Seed; run i uses seed + i. Default 1.",
    "  --norms=<a,b,c>       Comma separated normalization_methods for tune(); default",
    "                          log2,RLE,none.",
    "  --imputes=<a,b,c>     Comma separated impute_methods; default",
    "                          none,unif_sample_lod,sample_lod.",
    "  --tests=<a,b,c>       Comma separated test_methods; default lm,trend,voom.",
    "  --iquants=<a,b>       Comma separated impute_quantiles; default 0. Only the",
    "                          quantile and scale grids are settable here; tune()'s",
    "                          defaults are used for spans, npcs and ks.",
    "  --iscales=<a,b>       Comma separated impute_scales; default 1.",
    "  --prefix=<str>        Result filename prefix; default sim_int.",
    "  --suffix=<str>        Result filename suffix; default .sexage.tune.tsv.",
    "  --fdr_cutoff=<num>    Cutoff tune_check() sorts on; default 0.05.",
    "",
    "Output: in <out_dir>, one result table per run plus the check and the timings.",
    "  <prefix>0<suffix>       Unpermuted tune() result, one row per combination.",
    "  <prefix>{1..n}<suffix>  Permuted tune() results, same shape.",
    "  tune_check.tsv          tune_check() output, sorted, with fdr and the permuted",
    "                            hit summaries.",
    "  timing.tsv              run, permute_var, combinations, seconds, sec_per_comb.",
    "  h0testr.<run>.log       The package's own log for that run, one file per run.",
    "  run.log                 This script's stdout, which is also shown on stdout.",
    "",
    "Progress: a line per run start and finish; the per-combination detail is in the",
    "  h0testr.<run>.log for the run in progress, which grows as the sweep proceeds.",
    "",
    "Exit codes: 0 on success, 2 on a usage error, 3 if <dir_in> lacks a required file,",
    "  4 if <out_dir> cannot be created, 5 if sim.rds cannot be read, 6 if a tune() run",
    "  fails, 7 if tune_check() fails, 8 if an output file cannot be written.",
    "",
    "Examples:",
    "  Rscript run_sim_int_tune.R C:/path/to/h0testr/R ./sim_int ./sim_int_out",
    "  Rscript run_sim_int_tune.R ../../h0test/h0testr/R ./sim_int ./sim_int_out --n_perm=1",
    "  Rscript run_sim_int_tune.R C:/path/to/h0testr/R ./sim_int ./sim_full \\",
    "    --norms=log2,RLE,TMM,vsn,none --imputes=none,unif_sample_lod,qrilc,knn \\",
    "    --tests=lm,trend,deqms --iquants=0,0.05",
    "",
    sep="\n", file=stderr()
  )
  quit(save="no", status=2)
}

args <- commandArgs(trailingOnly=TRUE)
if(length(args) < 3) usage("wrong number of arguments")

r_dir <- args[1]
dir_in <- args[2]
out_dir <- args[3]
if(!dir.exists(r_dir)) usage(paste("r_dir not a directory:", r_dir))
if(!dir.exists(dir_in)) usage(paste("dir_in not a directory:", dir_in))

opts <- list(n_perm="3", permute_var="sex", seed="1", norms="log2,RLE,none",
  imputes="none,unif_sample_lod,sample_lod", tests="lm,trend,voom", iquants="0",
  iscales="1", prefix="sim_int.", suffix=".sexage.tune.tsv", fdr_cutoff="0.05")

for(arg in args[-(1:3)]) {
  if(!grepl("^--[a-z_]+=", arg)) usage(paste("unrecognized argument:", arg))
  nom <- sub("^--([a-z_]+)=.*$", "\\1", arg)
  if(!(nom %in% names(opts))) usage(paste("unrecognized argument:", arg))
  opts[[nom]] <- sub("^--[a-z_]+=", "", arg)
}

as_int <- function(nom) {
  v <- suppressWarnings(as.numeric(opts[[nom]]))
  if(is.na(v) || v != round(v)) usage(paste0("--", nom, " not a whole number: ", opts[[nom]]))
  return(as.integer(v))
}

as_nums <- function(nom) {
  v <- suppressWarnings(as.numeric(strsplit(opts[[nom]], ",", fixed=TRUE)[[1]]))
  if(!length(v) || any(is.na(v))) usage(paste0("--", nom, " not numeric: ", opts[[nom]]))
  return(v)
}

n_perm <- as_int("n_perm")
seed <- as_int("seed")
if(n_perm < 1) usage(paste("--n_perm must be at least 1; got:", opts$n_perm))
fdr_cutoff <- as_nums("fdr_cutoff")
if(length(fdr_cutoff) != 1 || fdr_cutoff <= 0 || fdr_cutoff > 1) {
  usage(paste("--fdr_cutoff must be one value in (0, 1]; got:", opts$fdr_cutoff))
}

iquants <- as_nums("iquants")
iscales <- as_nums("iscales")
norms <- strsplit(opts$norms, ",", fixed=TRUE)[[1]]
imputes <- strsplit(opts$imputes, ",", fixed=TRUE)[[1]]
tests <- strsplit(opts$tests, ",", fixed=TRUE)[[1]]
if(!nzchar(opts$suffix)) usage("--suffix must not be empty; tune_check() requires one")

need <- c("expression.tsv", "features.tsv", "samples.tsv", "sim.rds")
miss <- need[!file.exists(file.path(dir_in, need))]
if(length(miss)) {
  cat("ERROR: dir_in", dir_in, "lacks:", paste(miss, collapse=", "), "\n", file=stderr())
  quit(save="no", status=3)
}

for(f in list.files(r_dir, pattern="[.]R$", full.names=TRUE)) source(f)

if(!dir.exists(out_dir) && !dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)) {
  cat("ERROR: cannot create out_dir:", out_dir, "\n", file=stderr())
  quit(save="no", status=4)
}

log_out <- file.path(out_dir, "run.log")
sink(log_out, split=TRUE)

bail <- function(status, ...) {
  cat("ERROR:", paste0(...), "\n")
  while(sink.number() > 0) sink()
  cat("ERROR:", paste0(...), "\n", file=stderr())
  quit(save="no", status=status)
}

say <- function(...) {
  cat(format(Sys.time(), "%H:%M:%S"), " ", paste0(...), "\n", sep="")
  utils::flush.console()
}

###############################################################################
## configuration, from the fixture rather than rebuilt, so the two cannot drift:

sim <- try(readRDS(file.path(dir_in, "sim.rds")), silent=TRUE)
if(inherits(sim, "try-error") || is.null(sim$config)) {
  bail(5, "cannot read sim.rds in ", dir_in, ": ", as.character(sim))
}

cfg0 <- sim$config
cfg0$dir_in <- dir_in
cfg0$dir_out <- out_dir
cfg0$data_file_in <- "expression.tsv"
cfg0$feature_file_in <- "features.tsv"
cfg0$sample_file_in <- "samples.tsv"
cfg0$save_state <- FALSE

if(!(opts$permute_var %in% names(sim$state$samples))) {
  bail(3, "--permute_var names no column of samples.tsv: ", opts$permute_var,
    "; columns: ", paste(names(sim$state$samples), collapse=", "))
}

## the combinations a sweep of these lists costs, by the same branches tune() takes;
##   the grids tune() is not given here are read from its own defaults:

fm <- formals(tune)
n_grid <- c(loess_logit=length(eval(fm$impute_spans)),
  bpca=length(eval(fm$impute_npcs)), ppca=length(eval(fm$impute_npcs)),
  svdImpute=length(eval(fm$impute_npcs)), knn=length(eval(fm$impute_ks)),
  lls=length(eval(fm$impute_ks)))

n_cells <- sum(sapply(imputes, function(im) {
  if(im %in% c("unif_global_lod", "unif_sample_lod", "min_det")) return(length(iquants))
  if(im %in% c("qrilc", "rnorm_feature")) return(length(iscales))
  if(im %in% "min_prob") return(length(iquants) * length(iscales))
  if(im %in% names(n_grid)) return(n_grid[[im]])
  return(1)
}))
n_comb <- length(norms) * length(tests) * n_cells

say("fixture: ", dir_in)
say("normalizations: ", paste(norms, collapse=","))
say("imputations: ", paste(imputes, collapse=","), "  (quantiles ",
  paste(iquants, collapse=","), "; scales ", paste(iscales, collapse=","), ")")
say("tests: ", paste(tests, collapse=","))
say("combinations per run: ", n_comb, "; runs: ", n_perm + 1, "; total: ",
  n_comb * (n_perm + 1))

## one R session attaches every engine; running several of these in it has segfaulted:

risky <- intersect(tests, c("deqms", "msqrob", "msqrob_agg", "prolfqua", "prolfqua_lmer"))
if(length(risky) > 1) {
  say("WARNING: tests names ", length(risky), " of the heavy engines (",
    paste(risky, collapse=", "), "); those have segfaulted when run in one R session, ",
    "so consider one engine per invocation")
}

###############################################################################
## the runs:

timing <- NULL

for(i in 0:n_perm) {

  cfg <- cfg0
  cfg$permute_var <- if(i %in% 0) "" else opts$permute_var
  cfg$log_file <- file.path(out_dir, paste0("h0testr.", i, ".log"))
  f_out <- file.path(out_dir, paste0(opts$prefix, i, opts$suffix))

  say("run ", i, " of ", n_perm, ": permute_var '", cfg$permute_var, "'; log ",
    basename(cfg$log_file))

  set.seed(seed + i)
  t0 <- Sys.time()
  rslt <- try(tune(cfg, normalization_methods=norms, impute_methods=imputes,
    impute_quantiles=iquants, impute_scales=iscales, test_methods=tests), silent=TRUE)
  secs <- as.numeric(difftime(Sys.time(), t0, units="secs"))

  if(inherits(rslt, "try-error")) bail(6, "run ", i, " failed: ", as.character(rslt))

  ok <- try(utils::write.table(rslt, f_out, sep="\t", quote=FALSE, row.names=FALSE),
    silent=TRUE)
  if(inherits(ok, "try-error")) bail(8, "cannot write ", f_out, ": ", as.character(ok))

  say("run ", i, " done: ", nrow(rslt), " rows in ", round(secs, 1), "s (",
    round(secs / max(nrow(rslt), 1), 2), "s per combination); tested ",
    sum(!is.na(rslt$ntests)), " of ", nrow(rslt))

  timing <- rbind(timing, data.frame(run=i, permute_var=cfg$permute_var,
    combinations=nrow(rslt), seconds=round(secs, 1),
    sec_per_comb=round(secs / max(nrow(rslt), 1), 3), stringsAsFactors=FALSE))
}

ok <- try(utils::write.table(timing, file.path(out_dir, "timing.tsv"), sep="\t",
  quote=FALSE, row.names=FALSE), silent=TRUE)
if(inherits(ok, "try-error")) bail(8, "cannot write timing.tsv: ", as.character(ok))

###############################################################################
## the check, its messages going to stdout and so to run.log:

say("tune_check over ", n_perm + 1, " files in ", out_dir)

chk <- try(tune_check(out_dir, opts$prefix, opts$suffix, list(log_file=""),
  fdr_cutoff=fdr_cutoff), silent=TRUE)
if(inherits(chk, "try-error")) bail(7, "tune_check failed: ", as.character(chk))

ok <- try(utils::write.table(chk, file.path(out_dir, "tune_check.tsv"), sep="\t",
  quote=FALSE, row.names=FALSE), silent=TRUE)
if(inherits(ok, "try-error")) bail(8, "cannot write tune_check.tsv: ", as.character(ok))

say("wrote tune_check.tsv: ", nrow(chk), " rows, ", sum(!is.na(chk$fdr)), " with an fdr")
cat("\ntop combinations:\n")
print(utils::head(chk, 10))
cat("\ntotal elapsed: ", round(sum(timing$seconds), 1), "s\n", sep="")

while(sink.number() > 0) sink()
quit(save="no", status=0)
