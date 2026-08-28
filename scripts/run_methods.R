usage <- function(msg=NULL) {
  if(!is.null(msg)) cat("ERROR:", msg, "\n\n", file=stderr())
  cat(
    "Run one h0testr::run() pipeline on a simulated fixture with save_state TRUE, so every",
    "intermediate state and both result tables are written as files for inspection. One test",
    "method per invocation, because DEqMS, msqrob2 and prolfqua have segfaulted when loaded",
    "into the same R session; call this once per method, serially.",
    "",
    "Usage: Rscript run_methods.R <r_dir> <dir_in> <out_dir> --method=<name> [--frm=<formula>]",
    "         [--test_term=<term>] [--norm=<name>] [--impute=<name>] [--iquant=<num>]",
    "         [--fdr=<num>]",
    "",
    "Required positional arguments:",
    "  <r_dir>    Path to the h0testr package R/ source directory; every .R file there is",
    "               sourced, so the installed package is not used.",
    "  <dir_in>   Fixture directory from make_sim_int.R; needs expression.tsv, features.tsv,",
    "               samples.tsv, sim.rds and truth.tsv.",
    "  <out_dir>  Parent directory for results; a subdirectory named for the method is",
    "               created in it and everything for this run goes there.",
    "",
    "Required named argument:",
    "  --method=<name>     One of lm, trend, voom, proda, prolfqua, deqms, msqrob,",
    "                        prolfqua_lmer, msqrob_agg. The four gene level methods (deqms,",
    "                        msqrob, prolfqua_lmer, msqrob_agg) take feature level input, so",
    "                        combine_features is dropped from run_order for them.",
    "",
    "Optional named arguments:",
    "  --frm=<formula>     Model formula, no spaces; default ~sex+age.",
    "  --test_term=<term>  Term tested; default sex.",
    "  --norm=<name>       Normalization method; default log2.",
    "  --impute=<name>     Imputation method; default unif_sample_lod.",
    "  --iquant=<num>      Imputation quantile; default 0.",
    "  --fdr=<num>         Cutoff for counting hits, in (0, 1); default 0.05.",
    "",
    "Output: <out_dir>/<method>/ holding, from h0testr itself, three TSVs per pipeline step",
    "  (1.initial, 2.prepped, then one per run_order step), the two result tables named",
    "  <n>.results.reformat.tsv and <n>.results.original.tsv where n is length(run_order)+3,",
    "  and h0testr.log. This script adds summary.tsv, one row of method, gene_level, rows,",
    "  hits at the cutoff, planted genes recovered, seconds, and the columns of both tables.",
    "",
    "Exit codes: 0 on success, 2 on a usage error, 3 if <dir_in> lacks a needed file, 4 if",
    "  the output subdirectory cannot be created, 5 if sim.rds cannot be read, 6 if run()",
    "  fails, 7 if the summary cannot be written.",
    "",
    "Examples:",
    "  Rscript run_methods.R C:/path/to/h0testr/R ./sim_int ./sim_int_runs --method=lm",
    "  Rscript run_methods.R ../../h0test/h0testr/R ./sim_int ./runs --method=voom --fdr=0.01",
    "  Rscript run_methods.R C:/path/to/h0testr/R ./sim_int ./runs_int --method=trend",
    "    --frm=~sex*age --test_term=sex:age",
    "",
    sep="\n", file=stderr()
  )
  quit(save="no", status=2)
}

args <- commandArgs(trailingOnly=TRUE)
if(length(args) < 4) usage("wrong number of arguments")

r_dir <- args[1]
dir_in <- args[2]
out_dir <- args[3]
if(!dir.exists(r_dir)) usage(paste("r_dir not a directory:", r_dir))
if(!dir.exists(dir_in)) usage(paste("dir_in not a directory:", dir_in))

opts <- list(method="", frm="~sex+age", test_term="sex", norm="log2",
  impute="unif_sample_lod", iquant=0, fdr=0.05)
nums <- c("iquant", "fdr")

for(arg in args[-(1:3)]) {
  if(!grepl("^--[a-z_]+=", arg)) usage(paste("unrecognized argument:", arg))
  nom <- sub("^--([a-z_]+)=.*$", "\\1", arg)
  val <- sub("^--[a-z_]+=", "", arg)
  if(!(nom %in% names(opts))) usage(paste("unrecognized argument:", arg))
  if(nom %in% nums) {
    num <- suppressWarnings(as.numeric(val))
    if(is.na(num)) usage(paste("value not numeric:", arg))
    opts[[nom]] <- num
  } else {
    if(!nzchar(val)) usage(paste("value empty:", arg))
    opts[[nom]] <- val
  }
}

methods_ok <- c("lm", "trend", "voom", "proda", "prolfqua", "deqms", "msqrob",
  "prolfqua_lmer", "msqrob_agg")
if(!(opts$method %in% methods_ok))
  usage(paste0("--method= must name one of ", paste(methods_ok, collapse=", "),
    "; got: '", opts$method, "'"))
if(opts$fdr <= 0 || opts$fdr >= 1)
  usage(paste("--fdr= must be between 0 and 1, exclusive:", opts$fdr))

frm <- try(stats::as.formula(opts$frm), silent=TRUE)
if(inherits(frm, "try-error"))
  usage(paste("--frm= is not a formula:", opts$frm))

need <- c("expression.tsv", "features.tsv", "samples.tsv", "sim.rds", "truth.tsv")
miss <- need[!file.exists(file.path(dir_in, need))]
if(length(miss)) {
  cat("ERROR: dir_in", dir_in, "lacks:", paste(miss, collapse=", "), "\n", file=stderr())
  quit(save="no", status=3)
}

for(f in list.files(r_dir, pattern="[.]R$", full.names=TRUE)) source(f)

dir_run <- file.path(out_dir, opts$method)
if(!dir.exists(dir_run) && !dir.create(dir_run, recursive=TRUE, showWarnings=FALSE)) {
  cat("ERROR: cannot create output directory:", dir_run, "\n", file=stderr())
  quit(save="no", status=4)
}

say <- function(...) {
  cat(format(Sys.time(), "%H:%M:%S"), " ", paste0(...), "\n", sep="")
  utils::flush.console()
}

###############################################################################
## config: the fixture's own, with the run for this method laid over it

sim <- try(readRDS(file.path(dir_in, "sim.rds")), silent=TRUE)
if(inherits(sim, "try-error")) {
  cat("ERROR: cannot read sim.rds:", as.character(sim), "\n", file=stderr())
  quit(save="no", status=5)
}

gene_level <- f.gene_level_method(opts$method)

config <- sim$config
config$dir_in <- dir_in
config$dir_out <- dir_run
config$data_file_in <- "expression.tsv"
config$feature_file_in <- "features.tsv"
config$sample_file_in <- "samples.tsv"
config$log_file <- file.path(dir_run, "h0testr.log")
config$save_state <- TRUE
config$frm <- frm
config$test_term <- opts$test_term
config$normalization_method <- opts$norm
config$impute_method <- opts$impute
config$impute_quantile <- opts$iquant
config$test_method <- opts$method

## gene level methods report per gene from feature level input, so they must not
##   see the features already aggregated:
config$run_order <- c("normalize", "combine_replicates",
  if(!gene_level) "combine_features", "filter_state", "impute")

say("method ", opts$method, " (", if(gene_level) "gene" else "observation",
  " level), frm ", deparse(frm), ", test_term ", opts$test_term)
say("run_order: ", paste(config$run_order, collapse=", "))
say("output to ", dir_run)

t0 <- Sys.time()
rslt <- try(run(config), silent=TRUE)
secs <- round(as.numeric(difftime(Sys.time(), t0, units="secs")), 1)

if(inherits(rslt, "try-error")) {
  cat("ERROR: run() failed for method", opts$method, "after", secs, "s:\n",
    as.character(rslt), "\n", file=stderr())
  quit(save="no", status=6)
}

say("run() finished in ", secs, " s")

###############################################################################
## what the result table says, against what was planted

std <- rslt$standard
n_rows <- if(is.null(std)) 0L else nrow(std)
n_hits <- NA_integer_
n_recov <- NA_integer_
n_planted <- NA_integer_

if(!is.null(std) && "adj_pval" %in% names(std)) {
  hits <- std$feature[!is.na(std$adj_pval) & std$adj_pval < opts$fdr]
  n_hits <- length(hits)

  truth <- utils::read.table(file.path(dir_in, "truth.tsv"), header=TRUE, sep="\t",
    quote="", as.is=TRUE, check.names=FALSE)

  ## a truth column belongs to the term iff it involves the same variables, so
  ##   "sex" claims sexM but not sexM:agemid, which "sex:age" claims instead:
  vars <- all.vars(sim$config$frm)
  col_vars <- function(nm) {
    parts <- unlist(strsplit(nm, ":", fixed=TRUE))
    sort(unique(unlist(lapply(parts, function(p) vars[startsWith(p, vars)]))))
  }
  want <- sort(unlist(strsplit(opts$test_term, ":", fixed=TRUE)))
  cols <- setdiff(names(truth), "gene_id")
  cols <- cols[sapply(cols, function(nm) identical(col_vars(nm), want))]

  if(length(cols)) {
    planted <- truth$gene_id[rowSums(abs(truth[, cols, drop=FALSE]) > 0) > 0]
    n_planted <- length(planted)
    n_recov <- length(intersect(hits, planted))
  }
}

say("rows ", n_rows, "; hits at fdr<", opts$fdr, " ", n_hits, "; planted recovered ",
  n_recov, " of ", n_planted)

ok <- try(silent=TRUE, {
  smry <- data.frame(method=opts$method, gene_level=gene_level,
    frm=paste(deparse(frm), collapse=""), test_term=opts$test_term,
    norm=opts$norm, impute=opts$impute, rows=n_rows, hits=n_hits,
    recovered=n_recov, planted=n_planted, seconds=secs,
    std_cols=paste(names(std), collapse=","),
    orig_cols=paste(names(rslt$original), collapse=","),
    stringsAsFactors=FALSE)
  utils::write.table(smry, file.path(dir_run, "summary.tsv"), sep="\t", quote=FALSE,
    row.names=FALSE)
})

if(inherits(ok, "try-error")) {
  cat("ERROR: cannot write summary:", as.character(ok), "\n", file=stderr())
  quit(save="no", status=7)
}

say("wrote ", file.path(dir_run, "summary.tsv"))
quit(save="no", status=0)
