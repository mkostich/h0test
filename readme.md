# h0testR
## Hypothesis testing for ms-based proteomics

This package contains functions for analysis of non-negative continuous 
expression matrices, such as those produced by mass-spec-based proteomics.

The package was initially developed as a test platform for evaluating 
performance of different workflow configurations and parameter
settings. 

---

## Install

```
## in bash: download under ~/opt/h0test:
mkdir -p ~/opt
cd ~/opt
git clone https://github.com/mkostich/h0test

## in R: install dependencies: 
install.packages("randomForest")
install.packages("glmnet")
install.packages("BiocManager")
BiocManager::install("limma")
BiocManager::install("edgeR")

## install h0testr package and view manual:
install.packages("~/opt/h0test/h0testr", repo=NULL, type="source")
help(package="h0testr")
```

---

## Run

```
## get a default configuration object:
config <- h0testr::f.new_config()
str(config)                           ## view default settings

## customize input/output files:
config$out_dir <- "."                                ## directory for output
config$in_dir <- "/path/to/input/file/dir"           ## directory to input
config$feature_file_in <- "my_feature_metadata.tsv"  ## expected in config$in_dir
config$sample_file_in <- "my_sample_metadata.tsv"    ## expected in config$in_dir
config$data_file_in <- "my_signal_data.tsv"          ## expected in config$in_dir

## customize file cross-referencing:
config$feat_id_col <- "gene_id"       ## column in feature metadata matching rownames of sample_file_in
config$obs_id_col <- "replicate_id"   ## column in sample metadata matching column names of sample_file_in
config$sample_id_col <- "sample_id"   ## column in sample metadata; set obs_id_col == sample_id_col if no technical replicates

## customize testing configuration:
config$frm=~age+gender+age:gender     ## formula with variable of interest and covariates
config$test_term="age:gender"         ## term in config$frm on which test is to be performed
config$sample_factors=list(           ## set levels of factor variables in sample metadata
  age=c("young", "old"),    
  gender=c("Male", "Female")
)

## run workflow, generating hit table:
tbl <- h0testr::run(config)
head(tbl)
```
  
---

## Dependencies

- R, with utils, stats, methods, randomForest, glmnet, limma, edgeR
- Developed and tested with R 4.3.1, randomForest 4.7.1.2, glmnet 4.1.8, limma 3.56.2, and edgeR 3.42.4
- Dependencies (other than this package) available on JAX HPC in this apptainer container:<br/>
  `/projects/compsci/jgeorge/kostim/resources/containers/msproteomics_r.20250605a.sif`
- Should work on Linux, Mac, or Windows

---

## Inputs

Inputs are three tab-delimited text files in directory `config$dir_in`. The names of 
the files can be configured using parameters `config$feature_file_in`, 
`config$sample_file_in`, and `config$data_file_in`:

- `$feature_file_in`: feature metadata, with unique key `config$feat_id_col` 
     matching rownames of `$data_file_in`. Other columns contain supplementary 
     feature metadata that will be as feature annotation in final results table.
     Should not have rownames.
- `$sample_file_in`: observation metadata, without rownames, with:
  - `config$obs_id_col`: unique key; technical replicate identifier matching colnames of `data_file_in`.
  - `config$sample_id_col`: biosample identifier; not unique if technical replication
  - if you don't have technical replicates, set so `config$obs_id_col == config$sample_id_col`.
  - all variables referred to in the formula specified by configuration variable `config$frm`
- `$data_file_in`: raw (not log transformed) non-negative signal data table 
    with feature rownames and observation colnames.
  - rownames match the `config$feat_id_col` of `config$feature_file_in`.
  - colnames match the `config$obs_id_col` of `config$sample_file_in'.

---

## Configuration
  
Runs are configured by editing the default configuration returned by the 
`f.new_config()` function. The most frequently changed settings are towards 
the top of the returned list of defaults, and are briefly described below:

```
feature_file_in="features.tsv"  ## feature annotation .tsv; row features, metadata columns
sample_file_in="samples.tsv"    ## sample annotation .tsv; row observations, metadata columns
data_file_in="expression.tsv"   ## quantification matrix .tsv; rowname features, colname observations
dir_in="."                      ## feature_file_in, sample_file_in, and data_file_in found here (character)
dir_out="."                     ## output directory (character)
## formula for testing: actual formula can have '+' and ':'; not tested w/ e.g. '*' yet.
frm=~age+gender+age:gender      ## formula with variables of interest and covariates; all in sample_file_in
test_term="age:gender"          ## term (character) in $frm on which test is to be performed
permute_var=""                  ## name (character) of variable in $frm to permute; "" for no permutation (normal execution)
sample_factors=list(            ## set levels of factor variables in $frm
  age=c("young", "old"),        ## by default, numeric treated as numeric; if levels set here, treated as factor
  gender=c("Male", "Female")    ## by default, character treated as factor with alphabetically ordered levels
)
## samps and feats column names:
feat_id_col="gene_id"           ## column (character) in feature_file_in matching rownames of data_file_in
obs_id_col="observation_id"     ## column in sample_file_in matching colnames of data_file_in
sample_id_col="sample_id"       ## column in sample_file_in with sample ids; not unique if tech reps; same as obs_id_col if no tech reps
```

Here is the complete default configuration. The normalization, imputaton, 
filtering and testing options can be set further down in the list:

```
config <- list(

  ## input/output paths:
  feature_file_in="features.tsv",      ## feature annotation .tsv; row features
  sample_file_in="samples.tsv",        ## sample annotation .tsv; row observations
  data_file_in="expression.tsv",       ## quantification matrix .tsv; rowname features, colname observations
  dir_in=".",                          ## data_file_in, feature_file_in, and sample_file_in found here
  dir_out=".",                         ## output directory
  
  ## formula for testing: actual formula can have '+' and ':'; not tested w/ e.g. '*' yet.
  frm=~age+gender+age:gender,          ## formula with variable of interest and covariates
  test_term="age:gender",              ## term (character) in $frm on which test is to be performed
  permute_var="",                      ## name (character) of variable to permute; "" for no permutation (normal execution)
  sample_factors=list(                 ## set levels of factor variables in $frm
    age=c("young", "old"),             ## by default, numeric treated as numeric; if levels set here, treated as factor
    gender=c("Male", "Female")         ## by default, character treated as factor with alphabetically ordered levels
  ),
  
  ## samps and feats column names (new cols are introduced by the code):
  feat_id_col="gene_id",               ## column (character) in feature_file_in matching rownames of data_file_in
  obs_id_col="observation_id",         ## column in sample_file_in matching colnames of data_file_in
  sample_id_col="sample_id",           ## column in sample_file_in with sample ids; not unique if tech reps; same as obs_id_col if no tech reps
  obs_col="",                          ## for internal use; leave ""; samps[, obs_col] == colnames(exprs) throughout script
  n_samples_expr_col="n_samps_expr",   ## new col for feature metadata; n samples expressing feature
  median_raw_col="median_raw",         ## new col for feature metadata; median feature expression in expressing samples
  n_features_expr_col="n_feats_expr",  ## new col for sample metadata; n features expressed
  
  ## output file naming:
  log_file="",                         ## log file path (character); or "" for log to console                 
  feature_mid_out=".features",         ## midfix for output feature files
  sample_mid_out=".samples",           ## midfix for output samples file
  data_mid_out=".expression",          ## midfix for output expression files
  result_mid_out=".results",           ## midfix for output results file
  suffix_out=".tsv",                   ## suffix for output files  
  
  ## tunable options: defaults are usually ok, except:
  ##   for dia: usually works ok: RLE:unif_sample_lod:0.05 for norm_method:impute_method:impute_quantile
  ##   for dda: usually works ok: quantile:0.75:unif_sample_lod:0 for norm_method:norm_quantile:impute_method:impute_quantile
  norm_method="quantile",              ## in c("TMM", "TMMwsp", "RLE", "upperquartile", "quantile", "cpm", "vsn", "qquantile", "log2", "none")
  norm_quantile=0.75,                  ## for quantile normalization; 0.5 is median; 0.75 is upper quartile;
  n_samples_min=2,                     ## min samples/feature w/ feature expression > 0 to keep feature
  n_features_min=1000,                 ## min features/sample w/ expression > 0 to keep sample
  ## in: c("sample_lod", "unif_sample_lod", "unif_global_lod", "rnorm_feature", "glm_binom", "loess_logit", "glmnet", "rf", "none")
  impute_method="sample_lod",
  impute_quantile=0.01,                ## quantile for unif_ imputation methods
  impute_scale=1,                      ## for rnorm_feature, adjustment on sd of distribution [1: no change];
  impute_n_pts=1e7,                    ## granularity of imputed values for f.impute_glm_binom and f.impute_loess_logit
  impute_span=0.25,                    ## loess span for f.impute_loess_logit 
  test_method="trend",                 ## in c("voom", "trend")
  ## run_order character vector with elements from {"normalize", "combine_reps", "filter", "impute"}:
  run_order=c("normalize", "combine_reps", "filter", "impute"),   ## determines order of workflow operations
  
  ## misc; 
  probs=c(0, 0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99, 1.0),
  width=110,
  verbose=T
)
```

---

## Workflow overview

The workflow always begins with loading data and ends with testing of 
hypotheses. Intermediate steps can be configured using `config$run_order`. 
The default `config$run_order` of 
`c("normalize", "combine_reps", "filter", "impute")` yields the following
workflow:


1) Load data: read data from disk, format, [permute], check.
2) Inter-sample normalization: includes `log2(x+1)` transform
3) Combine technical replicates: using median
4) Filter features and samples: based on expression levels
5) Impute missing values: `0` and `NA` treated as missing
6) Test hypotheses

---

## Outputs

A series of intermediate files and a final results file are written to 
the output directory specified by `config$dir_out`. All these files are in 
tab-delimited text format. Each set of intermediate files consists of an 
expression matrix file, a feature metadata file, and a sample metadata file. 
The names of the files can be customized using `config$data_mid_out`, 
`config$feature_mid_out`, and `config$sample_mid_out` , respectively. The 
common file suffix can be customized by changing `config$suffix_out`.

Prefix numbering of output files reflects the order in which they were 
generated, which is determined by `config$run_order`. With the default settings 
for `config$data_mid_out`, `config$feature_mid_out`, `config$sample_mid_out`, 
`config$suffix_out`, and `config$run_order`, output files in `config$dir_out` 
will be:

```
## sample covariates reduced to set needed for downstream analysis;
##   features/observation and observations/feature recorded:
1.initial.expression.tsv
1.initial.features.tsv
1.initial.samples.tsv

## prefiltered uninformative features and observations; 
##   permute if PERMUTE_VAR not "":
2.prepped.expression.tsv
2.prepped.features.tsv
2.prepped.samples.tsv

## inter-sample normalized and transformed (e.g. log2):
3.normalized.expression.tsv
3.normalized.features.tsv
3.normalized.samples.tsv

## technical replicate data aggregated into sample-level data:
4.combined.expression.tsv
4.combined.features.tsv
4.combined.samples.tsv

## filtered for min samples/feature and features/sample:
5.filtered.expression.tsv
5.filtered.features.tsv
5.filtered.samples.tsv

## imputed missing values:
6.imputed.expression.tsv
6.imputed.features.tsv
6.imputed.samples.tsv

## test results:
7.results.tsv
```

The `results.tsv` with `config$test_term` set to `age` would look something like this:

```
results <- read.table("7.results.tsv", header=T, sep="\t", quote="", as.is=T)

> head(results)
   accession ... n_samps_expr median_raw     age12     age24  AveExpr        F      P.Value    adj.P.Val
1     K3BVX3 ...         202   7.449487 0.4403154 0.9871583 7.500060 67.11553 1.836730e-21 1.011487e-17
2     Q8JRU9 ...         202   7.106755 1.7900813 2.7292021 7.028900 58.77616 1.660120e-19 4.571139e-16
3     K9W4L5 ...         202   9.545797 0.6914303 1.2014321 9.552100 43.85464 1.118592e-15 2.053363e-12
4     E3TMK8 ...         202   9.141165 0.6385290 1.0643483 9.151548 42.81598 2.152679e-15 2.963701e-12
5     K6SUP6 ...         202   8.665571 0.6286665 0.8722942 8.642624 38.05398 4.674231e-14 5.148198e-11
6     P4TXL2 ...         202   9.486986 0.2715985 0.6570489 9.487649 37.33104 7.544069e-14 6.924198e-11

## where ... represents additional feature columns inserted from config$feature_file_in.
```

---

## NORM_METHOD

Values returned by each `config$norm_method` are `log2(x+1)` transformed 
unless `config$norm_method %in% "none"`:

**cpm**: counts per million; for each sample: `multiplier * (intensities / sum(intensities, na.rm=T))`

**quantile**: for each sample: `multiplier * intensities / quantile(intensities, probs=config$norm_quantile, na.rm=T))`

**RLE**: Relative log expression, using `edgeR::calcNormFactors()` then `edgeR::cpm()`. 
  See [https://doi.org/10.1038/npre.2010.4282.1](https://doi.org/10.1038/npre.2010.4282.1 "Anders and Huber, 2010")

**TMM**: Trimmed mean of medians normalization, using `edgeR::calcNormFactors()` then `edgeR::cpm()`. 
  See [https://doi.org/10.1186/gb-2010-11-3-r25](https://doi.org/10.1186/gb-2010-11-3-r25 "Robinson and Oshlack, 2010")
  
**TMMwsp**: TMM with singleton pairing, using `edgeR::calcNormFactors()` then `edgeR::cpm()`. 
  May be better than TMM when many zeros.

**qquantile**: classical quantile normalization, using `limma::normalizeQuantiles()`

**vsn**: variance stabilizing transformation, using `limma::normalizeVSN()`

**upperquartile**: upper quartile expression, using `edgeR::calcNormFactors()` then `edgeR::cpm()`.
  See [https://doi.org/10.1186/1471-2105-11-94](https://doi.org/10.1186/1471-2105-11-94 "Bullard et al., 2010")

**log2**: values are only `log(x+1)` transformed.

**none**: values are unchanged.

---

## IMPUTE_METHOD

Zero and `NA` values are both treated as missing. Methods `glmnet` and 
  `rf` are slow, and may take days to complete on large datasets.

**sample_lod**: for each sample, replace all missing values with the minimum 
  observed value in that sample.
  
**unif_sample_lod**: for each sample, calculate the `config$impute_quantile` 
  quantile `q` of observed values. Replace each `NA` value by randomly 
  sampling between `0` and `q`: `runif(1, min=0, max=q)`.

**unif_global_lod**: for each feature remaining after filtering, calculate the 
  lowest observed value. Calculate the `config$impute_quantile` quantile `q` 
  of the lowest observed values across all features. Replace each `NA` value 
  by randomly sampling between `0` and `q`: `runif(1, min=0, max=q)`.

**rnorm_feature**: for each feature, calculate the mean `m` and standard 
  deviation `s`. Replace each `NA` value with a random draw from a normal 
  distribution centered at `m` with standard deviation  
  `config$impute_scale * s`: `rnorm(1, mean=m, sd=s*config$impute_scale)
  
**glm_binom**: fit generalized linear model with logit link and binomial 
  errors to `p(missing) ~ log2(intensity)`, where `p(missing)` estimated 
  from number of missing observations for a feature, and `intensity` is
  median intensity in expressing (feature not missing) samples. Impute by
  drawing randomly from the resulting probability density.

**loess_logit**: fit locally linear model to 
  `logit(p(missing)) ~ log2(intensity)`, where `p(missing)` and `intensity` 
  estimated as for imputation by `glm_binom`. Impute by randomly drawing 
  from the resulting probability density.

**glmnet**: use the `glmnet::cv.glmnet` function to model missing feature 
  intensity as a function of the intensity of other features.

**rf**: use the `randomForest::randomForest` function to model missing 
  feature intensity as a function of the intensity of other features. 

**none**: no imputation. Depending on test method, features and/or samples
  with missing values may be dropped from the analysis.

---

## TEST_METHOD

Method used for hypothesis testing. 

**voom**: `limma::voom() -> limma::lmFit() -> limma::eBayes(trend=F) -> limma::topTable()`

**trend**: `limma::lmFit() -> limma::eBayes(trend=T) -> limma::topTable()`

---

## Tuning (alpha)

In practice, we run one run with no permuted variable, then around 20 runs 
with permuted variables. The naive code for doing so is shown below, which may 
be suitable for small to medium datasets. For larger datasets, we would 
normally launch a separate instance in its own directory, with one instance
for the unpermuted data, and 20 instances with permuted data:

```
## R:

rm(list=ls())

config <- h0testr::f.new_config()
config$feature_file_in <- "feats.tsv"
config$sample_file_in <- "samps.tsv"
config$data_file_in <- "exprs.tsv", 
config$feat_id_col <- "gene"
config$obs_id_col <- "replicate"
config$sample_id_col <- "sample"
config$frm <- ~age+strain+age:strain
config$test_term <- "age:strain"

## do this run with unpermuted variables once:
config$permute_var <- ""
tbl <- h0testr::f.tune(config)
write.table(tbl, file="0.age_strain.tune.tsv", sep="\t", quote=F, row.names=F)

## do this run, with permuted variable in config$test_term N times:
config$permute_var <- "age"
N <- 20
for(idx in 1:N) {
  tbl <- h0testr::f.tune(config)
  write.table(tbl, file=paste0(idx, ".age_sex.tune.tsv"), 
    sep="\t", quote=F, row.names=F)
}
```

The results of each run will look something like this:

```
> tbl
             norm norm_quant          impute imp_quant scale  test perm nhits ntests     time
1             TMM       0.75      sample_lod      0.01     1  voom       1993   7897 08:25:16
2             TMM       0.75 unif_sample_lod      0.00     1  voom       2063   7897 08:25:21
3             TMM       0.75 unif_sample_lod      0.01     1  voom       1513   7897 08:25:27
4             TMM       0.75 unif_sample_lod      0.05     1  voom        941   7897 08:25:32
5             TMM       0.75 unif_sample_lod      0.10     1  voom        633   7897 08:25:37
6             TMM       0.75 unif_sample_lod      0.25     1  voom        406   7897 08:25:42
7             TMM       0.75       glm_binom      0.25     1  voom        179   7897 08:25:48
8             TMM       0.75     loess_logit      0.25     1  voom        175   7897 08:25:53
```

Collect results. From parent directory, in R (only base packages required for this part):

```
rm(list=ls())

tbl <- f.tune_check(dir_in="/path/to/tuning/results", sfx=".age_strain.tune.tsv")

## tbl has a bunch of statistics that can be used to evaluate option combinations;
##   nhits: number of hits in unpermuted data
##   fdr: avg0 / nhits
##   max0: maximum number of hits in any permutation
##   mid0: median number of hits across permutations
##   avg0: mean number of hits across permutations
##   sd0: standard deviation of number of hits across permutations
##   norm: normalization method
##   norm_quant: normalization quantile (only used for quantile normalization)
##   impute: imputation method
##   impute_quant: imputation quantile for unif_ imputation methods
##   scale: scale parameter for rnorm_ imputation methods
##   test: method used for hypothesis testing

> head(tbl)
  nhits       fdr max0  mid0  avg0      sd0   norm norm_quant          impute imp_quant scale test
1  4366 0.1097572 2433  37.5 479.2 919.9720 TMMwsp       0.75      sample_lod      0.25   0.1 voom
2  4358 0.1150757 2538  70.5 501.5 908.2445 TMMwsp       0.75 unif_sample_lod      0.00   0.1 voom
3  4357 0.1245352 2660 134.0 542.6 911.6615 TMMwsp       0.75 unif_global_lod      0.25   0.1 voom
4  4357 0.1292862 2656 123.0 563.3 953.6246 TMMwsp       0.75 unif_global_lod      0.05   0.1 voom
5  4354 0.1289848 2648 129.0 561.6 951.6642 TMMwsp       0.75 unif_global_lod      0.00   0.1 voom
6  4352 0.1309053 2661 138.0 569.7 957.4347 TMMwsp       0.75 unif_global_lod      0.10   0.1 voom

## save results:
write.table(dat1, file="perm_results.age_strain.tsv", sep="\t", quote=F, row.names=F)
```

---
