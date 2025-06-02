## this file ...2.R meant to be invoked after invoking ...1.R, and before ...3.R
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
##   a) normalize exprs
##   b) transform (e.g. log2) exprs (if not done as part of normalize)
##   c) combine technical replicate observations; reduces number of columns
##        in exprs and number of rows in samps; exprs, samps, and feats in sync;

###############################################################################
## normalize:

f.msg("copying variables")
exprs2 <- exprs1
samps2 <- samps1
feats2 <- feats1

f.log("normalizing data")
f.msg("before normalization, N features: ", nrow(exprs2), "; N observations: ", ncol(exprs2))
f.msg("before normalization, signal distribution: ")
f.quantile(c(exprs2), probs=PROBS, digits=0)
f.msg("before normalization, min(exprs1):", min(c(exprs2), na.rm=T), "; mean(exprs2):", mean(c(exprs2), na.rm=T))
f.msg("before normalization, num NAs: ", sum(is.na(c(exprs2))))

if(NORM_METHOD %in% c("TMM", "TMMwsp", "RLE", "upperquartile")) {
  exprs2 <- f.normalize_edger(exprs2, method=NORM_METHOD, p=NORM_QUANTILE)
} else if(NORM_METHOD %in% "quantile") {
  exprs2 <- f.normalize_quantile(exprs2, norm_quantile=NORM_QUANTILE)
} else if(NORM_METHOD %in% "cpm") {
  exprs2 <- f.normalize_cpm(exprs2, multiplier=1e6)
} else if(NORM_METHOD %in% "vsn") {
  exprs2 <- f.normalize_vsn(exprs2)
} else if(NORM_METHOD %in% "qquantile") {
  exprs2 <- f.normalize_qquantile(exprs2)
} else if(NORM_METHOD %in% "none") {
  f.msg("skipping normalization: NORM_METHOD %in% 'none'")
} else f.err("unexpected NORM_METHOD:", NORM_METHOD)

if(!all(rownames(exprs2) == feats2[, FEAT_ID_COL, drop=T])) f.err("!all(rownames(exprs2) == feats2[, FEAT_ID_COL])")
if(!all(colnames(exprs2) == samps2[, OBS_COL, drop=T])) f.err("!all(colnames(exprs2) == samps2[, OBS_COL])")

f.msg("after normalization, N features: ", nrow(exprs2), "; N observations: ", ncol(exprs2))
f.msg("after normalization, signal distribution:")
f.quantile(c(exprs2), probs=PROBS, digits=3)
f.msg("after normalization, min(exprs2):", min(c(exprs2), na.rm=T), "; mean(exprs2):", mean(c(exprs2), na.rm=T))
f.msg("after normalization, num NAs: ", sum(is.na(c(exprs2))))


###############################################################################
## transform:

if(!(NORM_METHOD %in% "vsn" || TRANSFORM_METHOD %in% "none")) {
  f.log("transforming data")
  if(TRANSFORM_METHOD %in% "log2") {
    exprs2 <- log2(exprs2 + 1)
  } else if (TRANSFORM_METHOD %in% "log10") {
    exprs2 <- log10(exprs2 + 1)
  } else f.err("unexpected TRANSFORM_METHOD:", TRANSFORM_METHOD)

  if(!all(rownames(exprs2) == feats2[, FEAT_ID_COL, drop=T])) {
    f.err("!all(rownames(exprs2) == feats2[, FEAT_ID_COL])")
  } 
  if(!all(colnames(exprs2) == samps2[, OBS_COL, drop=T])) {
    f.err("!all(colnames(exprs2) == samps2[, OBS_COL])")
  } 
} else f.msg("skipping transformation: NORM_METHOD %in% 'vsn' || TRANSFORM_METHOD %in% 'none'")

## SAVE 2: normalize + transform:

f.msg("after transform, min(exprs2):", min(c(exprs2), na.rm=T), "; mean(exprs2):", mean(c(exprs2), na.rm=T))
f.msg("N features: ", nrow(exprs2), "; N observations: ", ncol(exprs2))
f.msg("signal distribution: ")
f.quantile(c(exprs2), probs=PROBS, digits=3)
f.msg("min(exprs2):", min(c(exprs2), na.rm=T), "; mean(exprs2):", mean(c(exprs2), na.rm=T))
f.msg("num NAs: ", sum(is.na(c(exprs2))))
f.save_main(exprs2, feats2, samps2, DIR_OUT, "3.normalized", SUFFIX_OUT)


###############################################################################
## aggregate technical replicate observations:

f.log("combining technical replicates")

if(!all(colnames(exprs2) == samps2[, OBS_COL, drop=T])) f.err("!all(colnames(exprs2)==samps2[, OBS_COL])")
if(!all(samps2$obs_f == colnames(exprs2))) f.err("!all(samps2$obs_f == colnames(exprs2))")

f <- function(v, s) tapply(v, s, median, na.rm=T)
exprs2 <- t(apply(exprs2, 1, f, samps2[, SAMPLE_COL, drop=T]))
rm(f)
samps2 <- samps2[!duplicated(samps2[, SAMPLE_COL, drop=T]), ]
samps2[, OBS_COL] <- NULL
if(!all(samps2[, SAMPLE_COL, drop=T] %in% colnames(exprs2))) {
  f.err("!all(samps2[, SAMPLE_COL] %in% colnames(exprs2))")
}

exprs2 <- exprs2[, samps2[, SAMPLE_COL, drop=T]]
if(!all(rownames(exprs2) == feats2[, FEAT_ID_COL, drop=T])) f.err("!all(rownames(exprs2) == feats2[, FEAT_ID_COL])")
if(!all(colnames(exprs2) == samps2[, SAMPLE_COL, drop=T])) f.err("!all(colnames(exprs2) == samps2[, SAMPLE_COL])")

## SAVE 3: collapse
f.msg("N features: ", nrow(exprs2), "; N observations: ", ncol(exprs2))
f.msg("signal distribution: ")
f.quantile(c(exprs2), probs=PROBS, digits=3)
f.msg("min(exprs2):", min(c(exprs2), na.rm=T), "; mean(exprs2):", mean(c(exprs2), na.rm=T))
f.msg("num NAs: ", sum(is.na(c(exprs2))))
f.save_main(exprs2, feats2, samps2, DIR_OUT, "4.combined", SUFFIX_OUT)

