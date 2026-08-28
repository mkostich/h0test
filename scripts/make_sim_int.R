usage <- function(msg=NULL) {
  if(!is.null(msg)) cat("ERROR:", msg, "\n\n", file=stderr())
  cat(
    "Simulate a peptide level dataset with a sex * age design and write it as h0testr",
    "input files, together with the true coefficients and a summary of what was made.",
    "Sex has 2 levels and age 3, both unordered factors; the term under test is the",
    "sex:age interaction. Effects are planted on the two main effects and on the",
    "interaction, one model matrix column per gene, with random sign.",
    "",
    "Usage: Rscript make_sim_int.R <r_dir> <out_dir> [--seed=<int>] [--n_per_cell=<int>]",
    "         [--n_genes=<int>] [--peps_per_gene=<int>] [--p_drop=<num>]",
    "         [--reps_per_sample=<int>] [--n_signif=<int>] [--log_cv_mean=<num>]",
    "         [--log_cv_sd=<num>] [--mnar_c0=<num>] [--fc_sex=<num>] [--fc_age=<num>]",
    "         [--fc_int=<num>]",
    "",
    "Required positional arguments:",
    "  <r_dir>    Path to the h0testr package R/ source directory; every .R file there is",
    "               sourced, so the installed package is not used.",
    "  <out_dir>  Directory the fixture is written to; created if absent. Existing files",
    "               of the same names are overwritten.",
    "",
    "Optional named arguments:",
    "  --seed=<int>             Seed, so the fixture is reproducible; default 1.",
    "  --n_per_cell=<int>       Samples per sex-by-age cell; default 6, so 36 samples.",
    "  --n_genes=<int>          Genes; default 300.",
    "  --peps_per_gene=<int>    Peptides simulated per gene before dropout; default 4.",
    "  --p_drop=<num>           Probability a peptide is dropped, each gene keeping at",
    "                             least one; default 0.1667, so about 1000 features.",
    "  --reps_per_sample=<int>  Technical replicates per sample; default 2, so 72",
    "                             observations from 36 samples.",
    "  --n_signif=<int>         Genes given an effect, per term; default 30, so 10% of",
    "                             300 on each of sex, age and sex:age.",
    "  --log_cv_mean=<num>      Mean of log(feature CV); default -2, so a CV near 0.135.",
    "                             The package default, -0.75, leaves the interaction",
    "                             undetectable at this replication.",
    "  --log_cv_sd=<num>        SD of log(feature CV), >= 0; default 0.25.",
    "  --mnar_c0=<num>          Intercept of logit(p(missing)) ~ c0 + c1 * log(intensity),",
    "                             c1 being -0.5; higher means more missing, -Inf none.",
    "                             Default 1.5, measured to give about 4% missing and 300",
    "                             complete features, so impute_method none has data while",
    "                             the imputers still have work to do.",
    "  --fc_sex=<num>           Fold change planted on sex, > 0; default 4.",
    "  --fc_age=<num>           Fold change planted on age, > 0; default 3.",
    "  --fc_int=<num>           Fold change planted on sex:age, > 0; default 3, chosen so",
    "                             the interaction is detectable, not to be plausible.",
    "",
    "Output: seven files in <out_dir>, and the same summary on stdout as in the last of",
    "  them. Planted coefficients are log2 of the --fc_* fold changes.",
    "  expression.tsv        Features by observations, integers with NA, row names in",
    "                          column 1 and no name in the header, as load_data() reads.",
    "  features.tsv          feature_id, gene_id.",
    "  samples.tsv           observation_id, sample_id, sex, age.",
    "  truth.tsv             gene_id and one column per model matrix column, of true",
    "                          log2 coefficients; 0 where no effect was planted.",
    "  truth_features.tsv    feature_id, gene_id, feat_mean, feat_cv.",
    "  sim.rds               The whole sim_design() return value, including the config",
    "                          that matches the fixture, for reuse without regenerating.",
    "  fixture_summary.txt   Dimensions, missingness, features per gene, complete",
    "                          feature count, and effects planted per column.",
    "",
    "Exit codes: 0 on success, 2 on a usage error, 3 if <out_dir> cannot be created,",
    "  4 if the simulation fails, 5 if an output file cannot be written.",
    "",
    "Examples:",
    "  Rscript make_sim_int.R C:/path/to/h0testr/R ./sim_int",
    "  Rscript make_sim_int.R ../../h0test/h0testr/R ./sim_int --seed=7 --fc_int=2",
    "  Rscript make_sim_int.R C:/path/to/h0testr/R ./sim_hard --log_cv_mean=-0.75",
    "    --mnar_c0=4.65 --n_per_cell=10",
    "",
    sep="\n", file=stderr()
  )
  quit(save="no", status=2)
}

args <- commandArgs(trailingOnly=TRUE)
if(length(args) < 2) usage("wrong number of arguments")

r_dir <- args[1]
out_dir <- args[2]
if(!dir.exists(r_dir)) usage(paste("r_dir not a directory:", r_dir))

opts <- list(seed=1, n_per_cell=6, n_genes=300, peps_per_gene=4, p_drop=1/6,
  reps_per_sample=2, n_signif=30, log_cv_mean=-2, log_cv_sd=0.25, mnar_c0=1.5,
  fc_sex=4, fc_age=3, fc_int=3)
ints <- c("seed", "n_per_cell", "n_genes", "peps_per_gene", "reps_per_sample", "n_signif")
fcs <- c("fc_sex", "fc_age", "fc_int")

for(arg in args[-(1:2)]) {
  if(!grepl("^--[a-z_]+=", arg)) usage(paste("unrecognized argument:", arg))
  nom <- sub("^--([a-z_]+)=.*$", "\\1", arg)
  val <- sub("^--[a-z_]+=", "", arg)
  if(!(nom %in% names(opts))) usage(paste("unrecognized argument:", arg))
  num <- suppressWarnings(as.numeric(val))
  if(is.na(num)) usage(paste("value not numeric:", arg))
  if(nom %in% ints && (!is.finite(num) || num != round(num)))
    usage(paste("value not a whole number:", arg))
  if(nom %in% fcs && (!is.finite(num) || num <= 0))
    usage(paste("fold change must be a finite number above 0:", arg))
  if(nom %in% "log_cv_sd" && (!is.finite(num) || num < 0))
    usage(paste("log_cv_sd must be a finite number at or above 0:", arg))
  if(nom %in% "log_cv_mean" && !is.finite(num))
    usage(paste("log_cv_mean must be finite:", arg))
  opts[[nom]] <- if(nom %in% ints) as.integer(num) else num
}

for(f in list.files(r_dir, pattern="[.]R$", full.names=TRUE)) source(f)

if(!dir.exists(out_dir) && !dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)) {
  cat("ERROR: cannot create out_dir:", out_dir, "\n", file=stderr())
  quit(save="no", status=3)
}

say <- function(...) {
  cat(format(Sys.time(), "%H:%M:%S"), " ", paste0(...), "\n", sep="")
  utils::flush.console()
}

###############################################################################
## the fixture:

say("simulating: ", opts$n_genes, " genes, ", opts$peps_per_gene, " peptides each, ",
  opts$n_per_cell, " samples per cell, ", opts$reps_per_sample, " reps per sample")

set.seed(opts$seed)

sim <- try(silent=TRUE, {

  samps <- sim_samples(factors=list(sex=c("F", "M"), age=c("young", "mid", "old")),
    n_per_cell=opts$n_per_cell)

  sim_design(samps, frm=~sex * age, test_term="sex:age",
    n_genes=opts$n_genes, peps_per_gene=opts$peps_per_gene, p_drop=opts$p_drop,
    reps_per_sample=opts$reps_per_sample, cv_reps=0.1,
    log_cv_mean=opts$log_cv_mean, log_cv_sd=opts$log_cv_sd, mnar_c0=opts$mnar_c0,
    n_genes_signif=c(sex=opts$n_signif, age=opts$n_signif, "sex:age"=opts$n_signif),
    effects=c(sex=log2(opts$fc_sex), age=log2(opts$fc_age),
      "sex:age"=log2(opts$fc_int)))
})

if(inherits(sim, "try-error")) {
  cat("ERROR: simulation failed:", as.character(sim), "\n", file=stderr())
  quit(save="no", status=4)
}

say("simulated ", nrow(sim$state$expression), " features by ",
  ncol(sim$state$expression), " observations")

###############################################################################
## the summary, which is both printed and saved:

e <- sim$state$expression
n_per_gene <- table(sim$feat_gene)
n_complete <- sum(apply(e, 1, function(v) !any(is.na(v))))
na_per_obs <- colMeans(is.na(e))
planted <- colSums(sim$truth != 0)

lines <- c(
  paste("fixture written by make_sim_int.R at", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  paste0("  seed=", opts$seed, " n_per_cell=", opts$n_per_cell, " n_genes=", opts$n_genes,
    " peps_per_gene=", opts$peps_per_gene, " p_drop=", signif(opts$p_drop, 4),
    " reps_per_sample=", opts$reps_per_sample, " n_signif=", opts$n_signif),
  paste0("  log_cv_mean=", opts$log_cv_mean, " log_cv_sd=", opts$log_cv_sd,
    " mnar_c0=", opts$mnar_c0, " fc_sex=", opts$fc_sex, " fc_age=", opts$fc_age,
    " fc_int=", opts$fc_int),
  "",
  paste0("features: ", nrow(e), " over ", length(n_per_gene), " genes"),
  paste0("observations: ", ncol(e), " from ",
    length(unique(sim$state$samples$sample_id)), " samples in ",
    nrow(unique(sim$state$samples[, c("sex", "age")])), " sex-by-age cells"),
  "",
  "features per gene:",
  paste0("  ", paste0(names(table(n_per_gene)), ":", as.integer(table(n_per_gene)),
    collapse="  ")),
  "",
  paste0("missing: ", round(100 * mean(is.na(e)), 2), "% of values; per observation ",
    round(100 * min(na_per_obs), 2), "% to ", round(100 * max(na_per_obs), 2), "%"),
  paste0("complete features: ", n_complete, " of ", nrow(e), " (",
    round(100 * n_complete / nrow(e), 1), "%), which is what voom keeps when nothing",
    " is imputed"),
  "",
  "effects planted, by model matrix column (genes with a non-zero true log2 coefficient):",
  paste0("  ", names(planted), ": ", as.integer(planted), collapse="\n"),
  paste0("effect sizes (log2): sex ", round(log2(opts$fc_sex), 4), ", age ",
    round(log2(opts$fc_age), 4), ", sex:age ", round(log2(opts$fc_int), 4)),
  paste0("feature CV: median ", round(exp(opts$log_cv_mean), 4), ", so about ",
    round(sqrt(log(1 + exp(opts$log_cv_mean)^2)) / log(2), 3), " sd on log2")
)

###############################################################################
## the files:

ok <- try(silent=TRUE, {

  utils::write.table(e, file.path(out_dir, "expression.tsv"), sep="\t", quote=FALSE,
    row.names=TRUE, col.names=TRUE)
  utils::write.table(sim$state$features, file.path(out_dir, "features.tsv"), sep="\t",
    quote=FALSE, row.names=FALSE)
  utils::write.table(sim$state$samples, file.path(out_dir, "samples.tsv"), sep="\t",
    quote=FALSE, row.names=FALSE)

  truth <- data.frame(gene_id=rownames(sim$truth), sim$truth, check.names=FALSE,
    stringsAsFactors=FALSE)
  utils::write.table(truth, file.path(out_dir, "truth.tsv"), sep="\t", quote=FALSE,
    row.names=FALSE)

  truth_feats <- data.frame(feature_id=names(sim$feat_gene),
    gene_id=as.character(sim$feat_gene), feat_mean=sim$feat_mean[names(sim$feat_gene)],
    feat_cv=sim$feat_cv[names(sim$feat_gene)], stringsAsFactors=FALSE)
  utils::write.table(truth_feats, file.path(out_dir, "truth_features.tsv"), sep="\t",
    quote=FALSE, row.names=FALSE)

  saveRDS(sim, file.path(out_dir, "sim.rds"))
  writeLines(lines, file.path(out_dir, "fixture_summary.txt"))
})

if(inherits(ok, "try-error")) {
  cat("ERROR: writing outputs failed:", as.character(ok), "\n", file=stderr())
  quit(save="no", status=5)
}

say("wrote 7 files to ", out_dir)
cat("\n", paste(lines, collapse="\n"), "\n", sep="")

quit(save="no", status=0)
