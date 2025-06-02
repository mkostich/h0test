## this file h0tune.1.R is called (sourced) in place of h0test.1.R for the 
##   purposes of tuning parameters; this script is meant to be invoked before 
##     invoking h0test.2.R, then h0test.3.R
##   this file preps data structures for h0test.2.R and h0test.3.R;
##   all configuration is done at the top of this file and in the calling
##     environment; 
##   tunable parameters should be set in calling environment and commmeted 
##     out here;
##   like h0test.1.R except f.save_main() only prints message (does not save
##     any files), and tunable paramter settings are commented out so they 
##     can be set in calling environment; also, configuration reporting is 
##     turned off;
##
## overall flow:
##   1) prefilter for all NA or all 0 rows (features) or columns (observations)
##   2) normalize, and if not part of normalization, transform (e.g. log2)
##   3) combine technical replicates (e.g. take median)
##   4) filter features and observations based on number of values > 0
##   5) impute missing values
##   6) test hypotheses, etc.
##
## this script:
##   a) CUSTOMIZE: configure variables
##   b) define functions
##   c) read data (exprs), sample metadata (meta), and feature metadata (feats);
##        exprs is numeric matrix, meta and feats are data frames; 
##        rownames(exprs) == feats[, FEAT_ID_COL]
##        colnames(exprs) == meta[, OBS_COL]
##   d) configure types and levels of covariates in meta;
##   e) add feature and sample metadata used for filtering:
##         i) samples per feature added as feats[, N_SAMPLES_EXPR_COL];
##        ii) median feature expression in expressing samples added as feats[, MEDIAN_RAW_COL];
##       iii) features per sample added as meta[, N_FEATURES_EXPR_COL];
##   f) prefilter features and samples:
##         i) prefilter features: n_expressing_samples >= 1; 
##        ii) prefilter samples: n_expressed_features >= 1;
##   g) permute PERMUTE_COL, if PERMUTE_COL != ""
##   Note: does not save anything to disk, but does define:
##     f.save_main(exprs, feats, meta, output_dir, prefix, suffix)

###############################################################################
## CONFIGURE THIS: variable initialization:

## TUNE SCRIPT: f.save_main() just prints message;
##   parameters to be tuned are commented out below, and are meant to be 
##     set by an outer loop, which tries different options;
##   define the following commented out parameters before calling this script:
##     TEST_TERM, PERMUTE_VAR, NORM_METHOD, NORM_QUANTILE, 
##       IMPUTE_METHOD, IMPUTE_SCALE, TEST_METHOD;

DIR_IN="/flashscratch/kostim/ps3/carter_glyco/2" ## where DATA_FILE_IN, FETURE_FILE_IN, and SAMPLE_FILE_IN live
DATA_FILE_IN <- "exprs_glyco_20250508c.tsv"      ## lfq normalized quant matrix .tsv; row features, column samples
FEATURE_FILE_IN <- "feats_glyco_20250508c.tsv"   ## feature annotation .tsv; row features
SAMPLE_FILE_IN <- "meta_glyco_20250508c.tsv"       ## sample annotation .tsv; row samples
DIR_OUT="."                                      ## output directory

## formula for testing:
FRM <- ~ age + strain + gender + age:strain      ## formula with variable of interest and covariates
## TEST_TERM <- "age:strain"                        ## term in FRM on which test is to be performed
## PERMUTE_VAR <- ""                                ## variable to permute; "" for no permutation (normal execution)
SAMPLE_FACTORS <- list(
  age=c("4mo", "12mo", "24mo"),
  strain=c("B6", "129", "A_J", "Balbc_J", "CAST", "NOD", "NZO", "PWK", "WSB"),
  gender=c("Male", "Female")
)

## meta and feats column names:
N_SAMPLES_EXPR_COL <- "n_samps_expr"             ## new col; n samples expressing feature
MEDIAN_RAW_COL <- "median_raw"                   ## new col; median feature expression in expressing samples
N_FEATURES_EXPR_COL <- "n_feats_expr"            ## new col; n features express
FEAT_ID_COL <- "MasterProteinAccessions"         ## feats[, FEAT_ID_COL] == rownames(exprs)
OBS_COL <- "obs_f"                               ## unique id for observations in meta
SAMPLE_COL <- "sample"                           ## unique id for samples in meta

## tunable parameters:
## NORM_METHOD <- "quantile"       ## in c("vsn","cpm","quantile","qquantile","TMM","TMMwsp","RLE","upperquartile")
## NORM_QUANTILE <- 0.75           ## for quantile normalization; 0.5 is median; 0.75 is upper quartile;
TRANSFORM_METHOD <- "log2"      ## in c('log2', 'log10', 'none')
N_SAMPLES_MIN <- 2              ## min samples/feature w/ feature expression > 0 to keep feature
N_FEATURES_MIN <- 1000          ## min features/sample w/ expression > 0 to keep sample
## IMPUTE_METHOD <- "sample_lod"   ## in c("unif_global_lod", "unif_sample_lod", "sample_lod", "rnorm_feature")
## IMPUTE_SCALE <- 1               ## for rnorm_feature, adjustment on sd of distribution [1: no change];
## TEST_METHOD <- "voom"           ## in c("voom", "trend")

## output file naming:
LOG_FILE <- "log.txt"           ## log file path; or "" for log to console                 
DATA_MID_OUT <- ".expression"   ## midfix for output expression files
FEATURE_MID_OUT <- ".features"  ## midfix for output feature files
SAMPLE_MID_OUT <- ".samples"      ## midfix for output metadata file
RESULT_MID_OUT <- ".results"    ## prefix for output results file
SUFFIX_OUT <- ".tsv"            ## suffix for output files

## required helper:
SRC_TEST="/projects/compsci/jgeorge/kostim/resources/scripts/mk_protools/h0test/h0test.0.R" 

## misc:
PROBS <- c(0, 0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99, 1.0)
WIDTH <- 110


###############################################################################
## set up:

options(width=WIDTH)
setwd(DIR_OUT)
source(SRC_TEST)

f.save_main <- function(exprs, feats, meta, output_dir, prefix, suffix) {
  cat("file saving turned off\n")
  return(NULL)
  file_out <- paste0(output_dir, "/", prefix, DATA_MID_OUT, suffix)
  f.log("writing exprs data to", file_out)
  f.save_tsv(exprs, file_out)

  file_out <- paste0(output_dir, "/", prefix, FEATURE_MID_OUT, suffix)
  f.log("writing feature metadata to", file_out)
  f.save_tsv(feats, file_out)

  file_out <- paste0(output_dir, "/", prefix, SAMPLE_MID_OUT, suffix)
  f.log("writing sample metadata to", file_out)
  f.save_tsv(meta, file_out)
}

f.quantile <- function(v, probs=PROBS, digits=3, na.rm=T) {
  if(LOG_FILE %in% "") {
    print(round(quantile(v, probs=probs, na.rm=na.rm), digits=digits))
  } else {
    capture.output(round(quantile(v, probs=probs, na.rm=na.rm), digits=digits),
      file=LOG_FILE, append=T)
  }
}


###############################################################################
## check TEST_TERM compatible with FRM:

f.msg("checking TEST_TERM compatible with FRM")
trms <- as.character(FRM)[2]
trms <- gsub("[[:space:]]+", "", trms)
trms <- unlist(strsplit(trms, split="\\+"))
if(!(TEST_TERM %in% trms)) f.err("TEST_TERM", TEST_TERM, "not in FRM", as.character(FRM)[2])
rm(trms)


###############################################################################
## read data:

f.log("reading data")
feats1 <- read.table(paste(DIR_IN, FEATURE_FILE_IN, sep="/"), header=T, sep="\t", quote="", as.is=T)
meta1 <- read.table(paste(DIR_IN, SAMPLE_FILE_IN, sep="/"), header=T, sep="\t", quote="", as.is=T)
exprs1 <- read.table(paste(DIR_IN, DATA_FILE_IN, sep="/"), header=T, sep="\t", quote="", as.is=T)
exprs1 <- as.matrix(exprs1)

if(!(typeof(exprs1) %in% "double")) f.err("!(typeof(exprs1) %in% 'double')")
if(!all(rownames(exprs1) == feats1[, FEAT_ID_COL, drop=T])) f.err("!all(rownames(exprs1)==feats1[,FEAT_ID_COL])")
if(!all(colnames(exprs1) == meta1[, OBS_COL, drop=T])) f.err("!all(colnames(exprs1)==meta1[, OBS_COL])")


###############################################################################
## set up types and levels of covariates:

for(nom in c(N_SAMPLES_EXPR_COL, MEDIAN_RAW_COL, N_FEATURES_EXPR_COL)) {
  if(nom %in% names(meta1)) f.err("nom %in% names(meta1); nom:", nom)
}
for(nom in c(OBS_COL, SAMPLE_COL)) {
  if(!(nom %in% names(meta1))) f.err("!(nom %in% names(meta1)); nom:", nom)
}
if(!(FEAT_ID_COL %in% names(feats1))) f.err("!(FEAT_ID_COL %in% names(feats1))")

## variables referred to in formula:
vars <- as.character(FRM)[2]
vars <- gsub("[[:space:]]+", "", vars)
vars <- gsub("[\\:\\*\\-]", "+", vars)
vars <- unlist(strsplit(vars, split="\\+"))
vars <- sort(unique(vars))

## make sure all needed variables in meta:
if(!all(vars %in% names(meta1))) f.err("!all(vars %in% names(meta1))")
if(!all(names(SAMPLE_FACTORS) %in% vars)) f.err("!all(names(SAMPLE_FACTORS) %in% vars)")

f.msg("subsetting metadata")
meta1 <- meta1[, c(OBS_COL, SAMPLE_COL, vars)]
rm(vars)

f.msg("setting factor levels")

for(nom in names(SAMPLE_FACTORS)) {
  ## check for potential misconfiguration first:
  if(!(nom %in% names(meta1))) f.err("!(nom %in% names(meta1)); nom:", nom)
  lvls1 <- SAMPLE_FACTORS[[nom]]
  lvls2 <- sort(unique(as.character(meta1[[nom]])))
  if(!all(lvls1 %in% lvls2)) f.err("!all(lvls1 %in% lvls2) for nom:", nom)
  if(!all(lvls2 %in% lvls1)) f.err("!all(lvls2 %in% lvls1) for nom:", nom)
  meta1[[nom]] <- factor(as.character(meta1[[nom]]), levels=SAMPLE_FACTORS[[nom]])
}


###############################################################################
## add stats to sample and feature metadata: 

n <- f.samples_per_feature(exprs1)
if(!all(names(n) == feats1[, FEAT_ID_COL, drop=T])) f.err("!all(names(n) == feats1[, FEAT_ID_COL])")
feats1[, N_SAMPLES_EXPR_COL] <- n

m <- f.feature_median_expression(exprs1)
if(!all(names(m) == feats1[, FEAT_ID_COL, drop=T])) f.err("!all(names(m) == feats1[, FEAT_ID_COL])")
feats1[, MEDIAN_RAW_COL] <- m
rm(m)

n <- f.features_per_sample(exprs1)
if(!all(names(n) == meta1[, OBS_COL, drop=T])) f.err("!all(names(n) == meta1[, OBS_COL])")
meta1[, N_FEATURES_EXPR_COL] <- n

n <- apply(exprs1, 1, function(v) sum(v > 0, na.rm=T))
f.msg("samples per feature:")
f.quantile(n, probs=PROBS, digits=0)

n <- apply(exprs1, 2, function(v) sum(v > 0, na.rm=T))
f.msg("features per sample")
f.quantile(n, probs=PROBS, digits=0)
rm(n)

## SAVE 1: initial
f.msg("N features: ", nrow(exprs1), "; N observations: ", ncol(exprs1))
f.msg("signal distribution:")
f.quantile(c(exprs1), probs=PROBS, digits=0) 
f.msg("min(exprs1):", min(c(exprs1), na.rm=T), "; mean(exprs1):", mean(c(exprs1), na.rm=T))
f.msg("num NAs: ", sum(is.na(c(exprs1))))
f.save_main(exprs1, feats1, meta1, DIR_OUT, "1.initial", SUFFIX_OUT)


###############################################################################
## prefilter:

f.log("prefilter features and samples")

f.msg("before filtering features", nrow(exprs1), "features, and", ncol(exprs1), "observations")
filter_list <- f.filter_features(exprs1, meta1, feats1, n_samples_min=1)
exprs1 <- filter_list$exprs
meta1 <- filter_list$meta
feats1 <- filter_list$feats
f.msg("after filtering features", nrow(exprs1), "features, and", ncol(exprs1), "observations")

filter_list <- f.filter_samples(exprs1, meta1, feats1, n_features_min=1)
exprs1 <- filter_list$exprs
meta1 <- filter_list$meta
feats1 <- filter_list$feats
f.msg("after filtering samples", nrow(exprs1), "features, and", ncol(exprs1), "observations")

if(!all(rownames(exprs1) == feats1[, FEAT_ID_COL, drop=T])) f.err("!all(rownames(exprs1) == feats1[, FEAT_ID_COL])")
if(!all(colnames(exprs1) == meta1[, OBS_COL, drop=T])) f.err("!all(colnames(exprs1) == meta1[, OBS_COL])")


###############################################################################
## permute (or not):

if(!(PERMUTE_VAR %in% "")) {
  f.msg("permuting", PERMUTE_VAR)
  tmp <- meta1[!duplicated(meta1[, SAMPLE_COL]), c(SAMPLE_COL, PERMUTE_VAR), drop=T]
  rownames(tmp) <- tmp[, SAMPLE_COL]
  tmp[, PERMUTE_VAR] <- sample(tmp[, PERMUTE_VAR, drop=T], nrow(tmp), replace=F)
  meta1[, PERMUTE_VAR] <- tmp[meta1[, SAMPLE_COL], PERMUTE_VAR, drop=T]
} else f.msg("skipping permutation; PERMUTE_VAR %in% ''")

if(!all(rownames(exprs1) == feats1[, FEAT_ID_COL, drop=T])) f.err("!all(rownames(exprs1) == feats1[, FEAT_ID_COL])")
if(!all(colnames(exprs1) == meta1[, OBS_COL, drop=T])) f.err("!all(colnames(exprs1) == meta1[, OBS_COL])")

f.save_main(exprs1, feats1, meta1, DIR_OUT, "2.prepped", SUFFIX_OUT)

