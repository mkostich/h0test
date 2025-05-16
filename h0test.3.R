## this file ...3.R meant to be invoked after invoking ...1.R, and ...2.R
##
## flow:
##   1) prefilter for all NA or all 0 rows (features) or columns (observations)
##   2) normalize, and if not part of normalization, transform (e.g. log2)
##   3) combine technical replicates (e.g. take median)
##   4) filter features and observations based on number of values > 0
##   5) impute missing values
##   6) test hypotheses, etc.
##
## this script:
##   a) filter features and samples at biosample level;
##   b) like ...1.R, adds filter related stats to feats[, N_SAMPLES_EXPR_COL], 
##        feats[, MEDIAN_RAW_COL], and meta[, N_FEATURES_EXPR_COL];
##   c) impute;
##   d) test;

###############################################################################
## filter at (bio)sample level: 

f.msg("copying variables")
exprs <- exprs2
feats <- feats2
meta <- meta2

f.log("filtering at sample level")

f.msg("before filtering features", nrow(exprs), "features, and", ncol(exprs), "observations")
filter_list <- f.filter_features(exprs, meta, feats, n_samples_min=N_SAMPLES_MIN)
exprs <- filter_list$exprs
meta <- filter_list$meta
feats <- filter_list$feats
f.msg("after filtering features", nrow(exprs), "features, and", ncol(exprs), "observations")
if(!all(rownames(exprs) == feats[, FEAT_ID_COL, drop=T])) f.err("!all(rownames(exprs) == feats[, FEAT_ID_COL])")
if(!all(colnames(exprs) == meta[, SAMPLE_COL, drop=T])) f.err("!all(colnames(exprs) == meta[, SAMPLE_COL])")

f.msg("before filtering samples", nrow(exprs), "features, and", ncol(exprs), "observations")
filter_list <- f.filter_samples(exprs, meta, feats, n_features_min=N_FEATURES_MIN)
exprs <- filter_list$exprs
meta <- filter_list$meta
feats <- filter_list$feats
f.msg("after filtering samples", nrow(exprs), "features, and", ncol(exprs), "observations")
if(!all(rownames(exprs) == feats[, FEAT_ID_COL, drop=T])) f.err("!all(rownames(exprs) == feats[, FEAT_ID_COL])")
if(!all(colnames(exprs) == meta[, SAMPLE_COL, drop=T])) f.err("!all(colnames(exprs) == meta[, SAMPLE_COL])")


###############################################################################
## add filter stats: 

n <- f.samples_per_feature(exprs)
if(!all(names(n) == feats[, FEAT_ID_COL, drop=T])) f.err("!all(names(n) == feats[, FEAT_ID_COL, drop=T])")
feats[, N_SAMPLES_EXPR_COL] <- n
f.msg("samples per feature:", round(quantile(n, probs=PROBS, na.rm=T)))

m <- f.feature_median_expression(exprs)
if(!all(names(m) == feats[, FEAT_ID_COL, drop=T])) f.err("!all(names(m) == feats[, FEAT_ID_COL, drop=T])")
feats[, MEDIAN_RAW_COL] <- m
rm(m)

n <- f.features_per_sample(exprs)
if(!all(names(n) == meta[, SAMPLE_COL])) f.err("!all(names(n) == meta[, SAMPLE_COL])")
meta[, N_FEATURES_EXPR_COL] <- n
f.msg("features per sample:", round(quantile(n, probs=PROBS, na.rm=T)))
rm(n)

f.msg("N features: ", nrow(exprs), "; N observations: ", ncol(exprs))
f.msg("signal distribution: ", quantile(c(exprs), probs=PROBS, na.rm=T))
f.msg("min(exprs):", min(c(exprs), na.rm=T), "; mean(exprs):", mean(c(exprs), na.rm=T))
f.msg("num NAs: ", sum(is.na(c(exprs))))
f.save_main(exprs, feats, meta, DIR_OUT, "5.filtered", SUFFIX_OUT)


###############################################################################
## impute:

f.log("imputing")
f.msg("before imputing, N features: ", nrow(exprs), "; N observations: ", ncol(exprs))
f.msg("before imputing, signal distribution:")
f.quantile(c(exprs), probs=PROBS, digits=3)
f.msg("before imputing, min(exprs):", min(c(exprs), na.rm=T), "; mean(exprs):", mean(c(exprs), na.rm=T))
f.msg("before imputing, num NAs: ", sum(is.na(c(exprs))))

if(IMPUTE_METHOD %in% "unif_global_lod") {
  exprs <- f.impute_unif_global_lod(exprs, impute_quantile=IMPUTE_QUANTILE)
} else if(IMPUTE_METHOD %in% "sample_lod") {
  exprs <- f.impute_sample_lod(exprs)
} else if(IMPUTE_METHOD %in% "unif_sample_lod") {
  exprs <- f.impute_unif_sample_lod(exprs, impute_quantile=IMPUTE_QUANTILE)
} else if(IMPUTE_METHOD %in% "rnorm_feature") {
  exprs <- f.impute_rnorm_feature(exprs, scale.=IMPUTE_SCALE)
} else if(IMPUTE_METHOD %in% "none") {
  f.msg("skipping imputation: IMPUTE_METHOD %in% 'none'")
} else f.err("unexpected IMPUTE_METHOD:", IMPUTE_METHOD)

f.msg("after imputing, N features: ", nrow(exprs), "; N observations: ", ncol(exprs))
f.msg("after imputing, signal distribution:")
f.quantile(c(exprs), probs=PROBS, digits=3)
f.msg("after imputing, min(exprs):", min(c(exprs), na.rm=T), "; mean(exprs):", mean(c(exprs), na.rm=T))
f.msg("after imputing, num NAs: ", sum(is.na(c(exprs))))
if(!all(rownames(exprs) == feats[, FEAT_ID_COL, drop=T])) f.err("!all(rownames(exprs) == feats[, FEAT_ID_COL])")
if(!all(colnames(exprs) == meta[, SAMPLE_COL, drop=T])) f.err("!all(colnames(exprs) == meta[, SAMPLE_COL])")
f.save_main(exprs, feats, meta, DIR_OUT, "6.imputed", SUFFIX_OUT)


###############################################################################
## test:

f.log("testing")

if(any(duplicated(feats[, FEAT_ID_COL, drop=T]))) f.err("any(duplicated(feats[, FEAT_ID_COL, drop=T]))")
rownames(feats) <- feats[, FEAT_ID_COL, drop=T]

if(TEST_METHOD %in% "voom") {
  tbl <- f.test_voom(exprs, meta, FRM, TEST_TERM, normalize.method="none")
} else if(TEST_METHOD %in% "trend") {
  tbl <- f.test_trend(exprs, meta, FRM, TEST_TERM)
} else if(TEST_METHOD %in% "none") {
  f.msg("skipping testing: TEST_METHOD %in% 'none'")
} else f.err("unexpected TEST_METHOD:", TEST_METHOD)

if(!all(rownames(tbl) %in% rownames(feats))) f.err("!all(rownames(tbl) %in% rownames(feats))")
tbl <- cbind(feats[rownames(tbl), ], tbl)
rownames(tbl) <- NULL
file_out <- paste0(DIR_OUT, "/", "7", RESULT_MID_OUT, SUFFIX_OUT)
f.log("writing results to", file_out)
f.save_tsv(tbl, file_out)

