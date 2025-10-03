# h0testR
## Hypothesis testing for ms-based proteomics

This package contains functions and workflows for analysis of non-negative 
continuous expression matrices, such as those produced by mass-spec-based 
proteomics.

In the proteomics field, there are many methods available for inter-sample 
normalization, feature aggregation, missing value imputation, and hypothesis 
testing. Many of the same methods are also applicable to, or have been 
borrowed from, other fields of expression data analysis, such as RNA-seq data. 
Many of these methods have parameters whose settings can affect the quality of 
the results.

This package was initially developed as a test platform for evaluating 
performance of different workflow configurations and parameter settings.
To facilitate achievement of this goal, it was designed to provide a uniform 
higher-level interface to a wide variety of potential methods, as well as a 
wrapper-based building-block design, which abstracts away many differences in 
native user interfaces between methods.

This building-block design facilitates use of individual methods, or 
configuring and running whole workflows. The workflow approach and the tuning 
functionality is described in the documentation below. See the `R` 
documentation for individual methods to see how to use them on their own.

---

## Install

```
## in bash: download under ~/opt/h0test:
mkdir -p ~/opt
cd ~/opt
git clone https://github.com/mkostich/h0test
R

## in R: install dependencies: 

install.packages(c("BiocManager", "glmnet", "imputeLCMD", "lmtest", 
  "missForest", "randomForest", "remotes"), dependencies=TRUE)

BiocManager::install(c("DEqMS", "edgeR", "limma", "impute", "MsCoreUtils", 
  "msqrob2", "pcaMethods", "proDA", "QFeatures", "SummarizedExperiment"))

remotes::install_github('fgcz/prolfqua', dependencies=TRUE)

## install h0testr package and view documentation:
install.packages("~/opt/h0test/h0testr", repo=NULL, type="source")
help(package="h0testr")
```

---

## Run

```
## get a default configuration object:
config <- h0testr::new_config()
str(config)                                          ## view default settings

## customize input/output files:
config$out_dir <- "."                                ## directory for output
config$in_dir <- "/path/to/input/file/dir"           ## directory to input
config$feature_file_in <- "my_feature_metadata.tsv"  ## feature metadata in config$in_dir
config$sample_file_in <- "my_sample_metadata.tsv"    ## sample metadata in config$in_dir
config$data_file_in <- "my_signal_data.tsv"          ## measurement data in config$in_dir

## customize file cross-referencing; sample_id_col becomes colnames of data_file_in
##   after aggregation of technical replicates; set obs_id_col == sample_id_col if 
##   no technical replicates:
config$feat_id_col <- "gene_id"       ## column in feature file matching rownames of data file
config$obs_id_col <- "replicate_id"   ## column in sample file matching colnames of data file
config$sample_id_col <- "sample_id"   ## column in sample file identifying samples

## customize testing configuration:
config$frm <- ~ age + gender + age:gender    ## formula to fit ('~+:' ok; '*' not tested)
config$test_term <- "age:gender"             ## term in config$frm to test
config$sample_factors=list(                  ## levels of factor variables in sample metadata
  age=c("young", "old"),    
  gender=c("Male", "Female")
)

## run workflow, generating hit table:
result <- h0testr::run(config)
head(result$standard)
```
  
---

## Dependencies

- R, with standard packages: utils, stats, methods; and add on packages: 
  DEqMS, edgeR, glmnet, impute, imputeLCMD, limma, missForest, MsCoreUtils, 
  msqrob2, pcaMethods, proDA, prolfqua, QFeatures, randomForest, and 
  SummarizedExperiment.
- Developed and tested with R 4.3.1.
- Expected to work on Linux, Mac, or Windows; tested on Rocky Linux 9.5 and 
  Windows 11.

---

## Inputs

Inputs are three tab-delimited text files in directory `config$dir_in`. By 
default, the names of these files are `features.tsv`, `samples.tsv`, and 
`expression.tsv`. The input file names can be changed using parameters 
`config$feature_file_in`, `config$sample_file_in`, and `config$data_file_in`,
respectively. You can place symbolic links to the input files, with these
default names in `config$dir_in` to reduce configuration customization.

Example input files can be found in `h0testr/inst/extdata`.

The `config$data_file_in` should contain non-negative raw (not 
log-transformed) expression data, with rows representing features (precursor, 
peptide, protein group, or gene group), and columns representing observations 
(fraction, technical replicate, or sample). This file is assumed to be saved 
using R with rownames and colnames. This results in the first row in the file 
containing the colnames, and the first column containing the rownames. 
**Row names and column names must be unique!**
The first row (containing colnames) has one less column than the rest of the 
table, since the first row has no corresponding rowname (since it does not
contain feature measurements), but all the following rows do. Row names 
correspond to unique feature identifiers and column names correspond to 
unique observation identifiers.

- `$feature_file_in`: feature metadata, with unique column name set by
     `config$feat_id_col`, and matching corresponding rownames of 
     `$data_file_in`. So, after loading of all three files completes, 
     `all(features[, feat_id_col] == rownames(data)`. In addition, this file 
     must contain a column with names specified by `config$gene_id_col`, 
     containing protein group ids, or gene group ids, which will be used to 
     aggregate feature signals to the final tested feature level (typically
     gene group or protein group). Other columns contain supplementary feature 
     metadata which will be appended as feature annotation columns in the 
     final results table. The input feature file should not have row names, 
     but must have column names.
- `$sample_file_in`: observation metadata, without row names, but with 
     column names. The names of these columns can be configured using the 
     following keys:
  - `config$obs_id_col`: unique key; technical replicate identifier matching 
    colnames of `data_file_in`.
  - `config$sample_id_col`: biosample identifier; not unique if technical 
    replication.
  - if you don't have technical replicates, set so that 
    `all(config$obs_id_col == config$sample_id_col)`.
  - all variables (column names) referred to in the formula specified by 
    configuration variable `config$frm`
- `$data_file_in`: raw (not log transformed) non-negative signal data table 
    with feature rownames and observation colnames. `0`s and `NA`s are treated 
    as missing values.
  - rownames match the `config$feat_id_col` of `config$feature_file_in`.
  - colnames match the `config$obs_id_col` of `config$sample_file_in'.

---

## Configuration
  
Runs are configured by editing the list object returned by the 
`new_config()` function, which contains the default configuration. The most 
frequently changed settings are towards the top of the returned list of 
defaults, and are briefly described below. **At a minimum**: edit the
the `frm`, `test_term`, and `sample_factors`, as the defaults are almost 
certainly wrong for your experiment.

```
## input/output paths:
feature_file_in="features.tsv"      ## feature annotation .tsv; row features
sample_file_in="samples.tsv"        ## sample annotation .tsv; row observations
data_file_in="expression.tsv"       ## quantification matrix .tsv; rowname features, colname observations
dir_in="."                          ## data_file_in, feature_file_in, and sample_file_in expected here
dir_out="."                         ## output directory

## match up expression dimnames with features metadata and samples metadata:
feat_id_col="feature_id"            ## column (scalar character) in feature_file_in matching rownames of data_file_in
gene_id_col="gene_id"               ## column (scalar character) in feature_file_in with gene group or protein group ids 
obs_id_col="observation_id"         ## column (scalar character) in sample_file_in matching colnames of data_file_in
sample_id_col="sample_id"           ## column (scalar character) in sample_file_in with sample ids; not unique if tech reps; same as obs_id_col if no tech reps
obs_col=""                          ## for internal use; leave ""; samps[, obs_col] == colnames(exprs) throughout script
feat_col=""                         ## for internal use; leave ""; feats[, feat_col] == rownames(exprs) throughout script

## formula for testing: actual formula can have '+' and ':'; not tested w/ e.g. '*' yet.
frm=~age+gender+age:gender          ## formula with variable of interest and covariates
test_term="age:gender"              ## term (scalar character) in $frm on which test is to be performed
permute_var=""                      ## name (scalar character) of variable to permute; "" for no permutation (normal execution)
sample_factors=list(                 ## set levels of factor variables in $frm
  age=c("young", "old"),             ## by default, numeric treated as numeric; if levels set here, treated as factor
  gender=c("Male", "Female")         ## by default, character treated as factor with alphabetically ordered levels
)
```

Here are the rest of the configuration options. The normalization, imputaton, 
filtering and testing options can be set further down in the list:

```
## new cols introduced into metadata data.frames by the code:
n_samples_expr_col="n_samps_expr"   ## new col (scalar character) for feature metadata; n samples expressing feature
median_raw_col="median_raw"         ## new col (scalar character) for feature metadata; median feature expression in expressing samples
n_features_expr_col="n_feats_expr"  ## new col (scalar character) for sample metadata; n features expressed

## output file naming:
log_file=""                         ## log file path (character); or "" for log to console                 
feature_mid_out=".features"         ## midfix for output feature files
sample_mid_out=".samples"           ## midfix for output samples file
data_mid_out=".expression"          ## midfix for output expression files
result_mid_out=".results"           ## midfix for output results file
suffix_out=".tsv"                   ## suffix for output files 

## tunable options: defaults are usually ok, except:
##   for dia: usually works ok: RLE:unif_sample_lod:0.05 for normalization_method:impute_method:impute_quantile
##   for dda: usually works ok: quantile:0.75:unif_sample_lod:0 for normalization_method:normalization_quantile:impute_method:impute_quantile
normalization_method="RLE"          ## normalization method; h0testr::normalize_methods() to see available choices.
normalization_quantile=0.75         ## for quantile normalization; 0.5 is median; 0.75 is upper quartile
normalization_span=0.7              ## span for normalization_method %in% "loess"
n_samples_min=2                     ## min(samples/feature w/ feature expression > 0) to keep feature
n_features_min=1000                 ## min(features/sample w/ expression > 0) to keep sample
feature_aggregation="medianPolish"  ## in c("medianPolish", "robustSummary", "none")
feature_aggregation_scaled=FALSE    ## whether to rescale peptide features prior to aggregation into protein/gene group.
impute_method="sample_lod"          ## imputation method; h0testr::impute_methods() to see available choices.
impute_quantile=0.01                ## quantile for unif_* imputation methods
impute_scale=1                      ## for rnorm_feature, adjustment on sd of distribution [1: no change];
impute_span=0.5                     ## loess span for impute_method %in% "loess_logit"
impute_k=7                          ## k for impute_method %in% c("knn", "lls")
impute_npcs=5                       ## number of PCs for PCA-based imputations "bpca", "ppca", and "svdImpute"
impute_alpha=1                      ## alpha mixing parameter for "glmnet" imputation
impute_n_pts=1e7                    ## granularity of imputed values for impute_method %in% c("glm_binom", "loess_logit")
impute_aug_steps=3                  ## data augmentation iterations for impute_rf() and impute_glmnet()
test_method="trend"                 ## hypothesis test; h0testr::test_methods() to see available choices.
test_prior_df=3                     ## prior df for test_method %in% "proda"
## run_order character vector with elements from {"normalize", "combine_replicates", "combine_features", "filter", "impute"}:
run_order=c("normalize", "combine_replicates", "combine_features", "filter", "impute")   ## order of workflow operations

## misc; 
save_state=TRUE                     ## whether to save output files; recommend FALSE for tuning/testing
probs=c(0, 0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99, 1.0)
width=110                           ## controls print width
verbose=T                           ## controls how much gets printed out during progress
```

---

## Workflow overview

The workflow always begins with loading data and ends with testing of 
hypotheses. Intermediate steps can be configured using `config$run_order`. 
The default `config$run_order` of 
`c("normalize", "combine_replicates", "filter", "impute")` yields the following
workflow:

1) Load data: read files `config$feature_file_in`, `config$sample_file_in`, 
   and `config$data_file_in` from the filesystem directory `config$dir_in`, 
   parse, reformat, [optionally permute], and check for consistency between
   all three files.
2) Inter-sample normalization specified by `config$normalization_method`: includes 
   `log2(x+1)` transform or equivalent (e.g. for `vsn`), unless 
   `config$normalization_method %in% "none"`.
3) Combine technical replicates: using median.
4) Unless `test_method %in% c("deqms", "msqrob")`, combine peptide signals 
   using method specified by `config$feature_aggregation`.
4) Filter features and samples: based on expression levels and cutoffs 
   specified by `config$n_samples_min` and `config$n_features_min`.
5) Impute missing values using method `config$impute_method`. Values `0` and 
   `NA` are treated as missing. 
6) Test null hypotheses that the effect of the term `config$test_term` in the 
   formula `config$frm` is `0` using the method `config$test_method`.

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
##   features per observation, median expresison, and observations 
##   per feature recorded:
1.initial.expression.tsv
1.initial.features.tsv
1.initial.samples.tsv

## prefiltered uninformative features and observations; 
##   permute if config$permute_var not "":
2.prepped.expression.tsv
2.prepped.features.tsv
2.prepped.samples.tsv

## inter-sample normalized and transformed (log2 or glog2):
3.normalized.expression.tsv
3.normalized.features.tsv
3.normalized.samples.tsv

## technical replicate data aggregated into sample-level data:
4.combined_replicates.expression.tsv
4.combined_replicates.features.tsv
4.combined_replicates.samples.tsv

## feature data aggregated into protein-group or gene-group-level data:
5.combined_features.expression.tsv
5.combined_features.features.tsv
5.combined_features.samples.tsv

## filtered for min samples/feature and features/sample:
6.filtered.expression.tsv
6.filtered.features.tsv
6.filtered.samples.tsv

## imputed missing values:
7.imputed.expression.tsv
7.imputed.features.tsv
7.imputed.samples.tsv

## test results:
8.results.tsv
```

The `results.tsv` with `config$test_method` set to `trend`, and 
  `config$test_term` set to `age` would look something like this:

```
results <- read.table("7.results.tsv", header=T, sep="\t", quote="", as.is=T)

> head(results)
   accession ... n_samps_expr median_raw     age12     age24  AveExpr        F      P.Value    adj.P.Val
1     K3BVX3 ...          202   7.449487 0.4403154 0.9871583 7.500060 67.11553 1.836730e-21 1.011487e-17
2     Q8JRU9 ...          202   7.106755 1.7900813 2.7292021 7.028900 58.77616 1.660120e-19 4.571139e-16
3     K9W4L5 ...          202   9.545797 0.6914303 1.2014321 9.552100 43.85464 1.118592e-15 2.053363e-12
4     E3TMK8 ...          202   9.141165 0.6385290 1.0643483 9.151548 42.81598 2.152679e-15 2.963701e-12
5     K6SUP6 ...          202   8.665571 0.6286665 0.8722942 8.642624 38.05398 4.674231e-14 5.148198e-11
6     P4TXL2 ...          202   9.486986 0.2715985 0.6570489 9.487649 37.33104 7.544069e-14 6.924198e-11

## where '...' represents feature columns from config$feature_file_in.
```

---

## NORMALIZATION_METHOD

Itersample normalization method. Values returned by each `config$normalization_method` 
are `log2(x+1)` transformed unless `config$normalization_method %in% "none"`:

**cpm**: Counts per million; for each sample: 
  `multiplier * (intensities / sum(intensities, na.rm=T))`
  
**div.mean**: Subtract an observation's mean intensity across features from 
  each feature intensity for that observation. Uses 
  `MsCoreUtils::normalize_matrix(..., method="div.mean")`.

**div.median**: Subtract an observation's median intensity across features 
  from each feature intensity for that observation. Uses 
  `MsCoreUtils::normalize_matrix(..., method="div.median")`.

**loess**: Cyclic loess normalization. Uses `limma::normalizeCyclicLoess()`.

**log2**: Values are simply `log(x+1)` transformed.

**max**: Divide each feature's intensities by the max for that feature. 
  Uses `MsCoreUtils::normalize_matrix(..., method="max")`.

**none**: Values unchanged. 

**qquantile**: Classical (from microarray days) quantile normalization, using 
  `limma::normalizeQuantiles()`.

**quantile**: For each sample: 
  `multiplier * intensities / quantile(intensities, probs=config$normalization_quantile, na.rm=T))`
  
**quantiles.robust**: Robust quantile normalization. Uses 
  `MsCoreUtils::normalize_matrix(..., method="quantiles.robust")`, which in turn is a 
  wrapper for `preprocessCore::normalize.quantiles.robust()`.

**RLE**: Relative log expression, using `edgeR::calcNormFactors()` then `edgeR::cpm()`. 
  See [https://doi.org/10.1038/npre.2010.4282.1](https://doi.org/10.1038/npre.2010.4282.1 "Anders and Huber, 2010")
  
**sum**: Divide each feature's intensities by the sum of intensities for that 
  feature. Uses `MsCoreUtils::normalize_matrix(..., method="sum")`.
  
**TMM**: Trimmed mean of medians normalization, using `edgeR::calcNormFactors()` then `edgeR::cpm()`. 
  See [https://doi.org/10.1186/gb-2010-11-3-r25](https://doi.org/10.1186/gb-2010-11-3-r25 "Robinson and Oshlack, 2010")
  
**TMMwsp**: TMM with singleton pairing, using `edgeR::calcNormFactors()` then `edgeR::cpm()`. 
  May be better than TMM when many zeros are present.

**upperquartile**: Upper quartile expression, using `edgeR::calcNormFactors()` then `edgeR::cpm()`.
  See [https://doi.org/10.1186/1471-2105-11-94](https://doi.org/10.1186/1471-2105-11-94 "Bullard et al., 2010")
  
**vsn**: Variance stabilizing transformation, using `limma::normalizeVSN()`

---

## IMPUTE_METHOD

Zero and `NA` values are both treated as missing. Methods `glmnet` and 
  `rf` are slow, and may take days to complete on large datasets.
  
**bpca**: Bayesian PCA (https://doi.org/10.1093/bioinformatics/btg287).
  Uses `pcaMethods::pca(..., method="bpca")`.

**glm_binom**: Fit a generalized linear model with logit link and binomial 
  errors to `p(missing) ~ log2(intensity)`, where `p(missing)` estimated 
  from number of missing observations for a feature, and `intensity` is
  median intensity in expressing (feature not missing) samples. Impute by
  drawing randomly from the resulting probability density.

**glmnet**: Use the `glmnet::cv.glmnet` function to model missing feature 
  intensity as a function of the intensity of other features.

**knn**: K-nearest neighbors. Uses `impute::impute.knn()`.

**lls**: Linear least squares imputation. Uses `pcaMethods::llsImpute()`. 

**loess_logit**: Fit a locally linear model to 
  `logit(p(missing)) ~ log2(intensity)`, where `p(missing)` and `intensity` 
  estimated as for imputation by `glm_binom`. Impute by randomly drawing 
  from the resulting probability density.

**min_det**: Impute using a low quantile of the observed signal 
  distribution (near putative limit of detection). 
  Uses `imputeLCMD::impute.MinDet()`.

**min_prob**: Impute using random draws from normal distribution near
  putative limit of detection. Uses `imputeLCMD::impute.MinProb()`.

**missforest**: Impute missing values with the missForest algorithm 
  (https://dx.doi.org/10.1093/bioinformatics/btr597).
  Uses `missForest::missForest()`.
  
**none**: No imputation. Depending on test method, features and/or samples
  with missing values may be dropped from the analysis. NOTE: Missing values 
  are not compatible with `config$test_method %in% c("msqrob", "voom")`.

**ppca**: Probabilistic PCA (https://dl.acm.org/doi/10.5555/3008904.3008993).
  Uses `pcaMethods::pca(..., method="ppca")`.

**qrilc**: Impute using quantile regression. Uses `imputeLCMD::impute.QRILC()`.

**rf**: Use the `randomForest::randomForest` function to model missing 
  feature intensity as a function of the intensity of other features. 

**rnorm_feature**: For each feature, calculate the mean `m` and standard 
  deviation `s`. Replace each `NA` value with a random draw from a normal 
  distribution centered at `m` with standard deviation  
  `config$impute_scale * s`: `rnorm(1, mean=m, sd=s*config$impute_scale)`.
  
**sample_lod**: For each sample, replace all missing values with the minimum 
  observed value in that sample.

**svdImpute**: Imputation with the svdImpute algorithm 
  (https://dx.doi.org/10.1093/bioinformatics/17.6.520).
  Uses `pcaMethods::pca(..., method="ppca")`.

**unif_global_lod**: For each feature remaining after filtering, calculate the 
  lowest observed value. Calculate the `config$impute_quantile` quantile `q` 
  of the lowest observed values across all features. Replace each `NA` value 
  by randomly sampling between `0` and `q`: `runif(1, min=0, max=q)`.
  
**unif_sample_lod**: For each sample, calculate the `config$impute_quantile` 
  quantile `q` of observed values. Replace each `NA` value by randomly 
  sampling between `0` and `q`: `runif(1, min=0, max=q)`.
  
---

## TEST_METHOD

Method used for hypothesis testing. 

**deqms**: Uses `DEqMS` package for peptide-based protein/gene-group analysis.
  Flow is: `limma::lmFit() -> limma::eBayes() -> DEqMS::spectraCounteBayes() 
  -> DEqMS::outputResult()`.

**msqrob**: Uses `msqrob2` package for peptide-based protein/gene-group analysis.
  Flow is: `QFeatures::readQFeatures() -> SummarizedExperiment object -> 
  QFeatures::zeroIsNA() -> QFeatures::aggregateFeatures() -> 
  msqrob2::msqrob() -> msqrob2::makeContrast() -> msqrob2::hypothesisTest()`.

**proda**: Uses `proDA::proDA()` for hypothesis testing.

**prolfqua**: Flow: `prolfqua::AnalysisTableAnnotation$new() -> 
  prolfqua::LFQData$new() -> prolfqua::strategy_lm() -> prolfqua::build_model() 
  -> model$get_anova()`.

**trend**: Flow: `limma::lmFit() -> limma::eBayes(trend=T) -> limma::topTable()`

**voom**: Flow: `limma::voom() -> limma::lmFit() -> limma::eBayes(trend=F) -> limma::topTable()`

---

## Tuning (beta)

In practice, perform one run with no permuted variable, then around 20 runs 
with permutation of a variable in the test term `config$test_term`. The 
naive code for doing so is shown below, which may be suitable for small to 
medium datasets. For larger datasets, we would normally launch a separate 
instance in its own directory, with one instance for the unpermuted data, 
and 20 instances with permuted data:

```
## R:

rm(list=ls())

## set up configuration:
config <- h0testr::new_config()   ## defaults
config$save_state <- FALSE          ## default is TRUE
config$dir_in <- system.file("extdata", package="h0testr")  ## where example data 
config$feature_file_in <- "features2.tsv"
config$sample_file_in <- "samples2.tsv"
config$data_file_in <- "expression2.tsv" 
config$feat_id_col <- "pep"
config$gene_id_col <- "gene"
config$obs_id_col <- config$sample_id_col <- "obs"
config$n_features_min <- 10         ## default 1000 too big for small demo dataset
config$frm <- ~grp
config$test_term <- "grp"
config$test_method <- "trend"
config$sample_factors <- list(grp=c("ctl", "trt"))

## one run with unpermuted data:
config$permute_var <- ""            ## no permutation
set.seed(100)
out <- h0testr::tune(config,
  normalization_methods=c("RLE", "q75", "cpm", "log2"),
  impute_methods=c("sample_lod", "unif_sample_lod", "none"),
  impute_quantiles=c(0, 0.05, 0.1),
  test_methods=c("trend", "msqrob", "proda", "prolfqua")
)
write.table(out, "0.grp.tune.tsv", quote=F, sep="\t", row.names=F)

## do another run with permuted variable in config$test_term N times;
##   since test_term is "grp", only variable to permute is "grp"; 
##   if test_term was instead "age:gender", could permute either "age" or 
##   "gender":

config$permute_var <- "grp"
N <- 20
for(idx in 1:N) {
  set.seed(100 + idx)
  out <- h0testr::tune(config,
    normalization_methods=c("RLE", "q75", "cpm", "log2"),
    impute_methods=c("sample_lod", "unif_sample_lod", "none"),
    impute_quantiles=c(0, 0.05, 0.1),
    test_methods=c("trend", "msqrob", "proda", "prolfqua")
  )
  write.table(out, file=paste0(idx, ".grp.tune.tsv"), 
    sep="\t", quote=F, row.names=F)
}
```

The results of each run will have columns:

**norm**: Normalization method. `config$normalization_method`.

**nquant**: Normalization quantile. `config$normalization_quantile`.

**impute**: Imputation method. `config$impute_method`.

**iquant**: Imputation quantile. `config$impute_quantile`.

**scale**: Imputation scale. `config$impute_scale`.

**span**: Imputation span. `config$impute_span`.

**npcs**: Number of PCs for PCA-based imputations. `config$impute_npcs`.

**k**: Number of nearest neighbors or group size for imputation. `config$impute_k`.

**test**: Hypothesis testing method. `config$test_method`.

**perm**: Variable in `config$test_term` to permute. `config$permute_var`.

**nhits**: Number of hypothesis tests yielding adjusted p-value less than `0.05`.

**ntests**: Number of hypothesis tests performed.

**time**: timestamp when iteration completed.

Example output from a more extensive parameter exploration is shown below:

```
> head(out) 
  norm nquant          impute iquant scale span npcs k  test perm nhits ntests     time
1  RLE   0.75      sample_lod   0.01     1  0.5    5 7 trend   NA    59   3940 17:00:28
2  RLE   0.75 unif_sample_lod   0.00     1  0.5    5 7 trend   NA    59   3940 17:00:29
3  RLE   0.75 unif_sample_lod   0.01     1  0.5    5 7 trend   NA    58   3940 17:00:29
4  RLE   0.75 unif_sample_lod   0.05     1  0.5    5 7 trend   NA    65   3940 17:00:29
5  RLE   0.75 unif_sample_lod   0.10     1  0.5    5 7 trend   NA    56   3940 17:00:30
6  RLE   0.75 unif_global_lod   0.00     1  0.5    5 7 trend   NA    60   3940 17:00:30

> tail(out)
     norm nquant    impute iquant scale span npcs  k test perm nhits ntests     time
4051 none   0.75 svdImpute    0.1   0.1 0.75    3 20 voom   NA  3939   3940 11:25:59
4052 none   0.75 svdImpute    0.1   0.1 0.75    5 20 voom   NA     0   3940 11:26:19
4053 none   0.75 svdImpute    0.1   0.1 0.75   10 20 voom   NA     0   3940 11:26:38
4054 none   0.75       lls    0.1   0.1 0.75   10  5 voom   NA     0   3940 11:26:40
4055 none   0.75       lls    0.1   0.1 0.75   10 10 voom   NA     0   3940 11:26:42
4056 none   0.75       lls    0.1   0.1 0.75   10 20 voom   NA     0   3940 11:26:44
```

Collect results. From parent directory, in R (only base packages required for this part):

```
rm(list=ls())

tbl <- h0testr::tune_check(dir_in="/path/to/tuning/results", prefix="", suffix=".grp.tune.tsv")

## tbl has a bunch of statistics that can be used to evaluate option combinations;
##   nhits: number of hits in unpermuted data
##   fdr: avg0 / nhits
##   max1: maximum number of hits in any permutation
##   mid1: median number of hits across permutations
##   avg1: mean number of hits across permutations
##   sd1: standard deviation of number of hits across permutations
##   norm: normalization method
##   nquant: normalization quantile (only used for quantile normalization)
##   impute: imputation method
##   iquant: imputation quantile for unif_ imputation methods
##   scale: scale parameter for rnorm_ imputation methods
##   span: span for loess-based methods
##   npcs: number of principle components for PCA/SVD-based imputations
##   k: number of nearest neighbors or groups for k-nn and related imputations
##   test: method used for hypothesis testing

> head(tbl)
  nhits         fdr max1 mid1 avg1       sd1     norm nquant    impute iquant scale span npcs  k     test
1  2011 0.000497265    1    0 0.05 0.2236068 quantile   0.50 svdImpute    0.1   0.1 0.75    3 20     voom
2   878 0.001138952    1    0 0.10 0.3077935      RLE   0.75 svdImpute    0.1   0.1 0.75    5 20     voom
3   533 0.001876173    1    0 0.05 0.2236068      RLE   0.75 svdImpute    0.1   0.1 0.75    3 20     voom
4   477 0.000000000    0    0 0.00 0.0000000 quantile   0.50      ppca    0.1   0.1 0.75    3 20 prolfqua
5   467 0.000000000    0    0 0.00 0.0000000      RLE   0.75      ppca    0.1   0.1 0.75    3 20 prolfqua
6   458 0.000000000    0    0 0.00 0.0000000     log2   0.75      ppca    0.1   0.1 0.75    3 20 prolfqua

> tail(tbl)
     nhits fdr max1 mid1   avg1          sd1 norm nquant    impute iquant scale span npcs  k   test
4051     0   1    3    0   0.15    0.6708204 none   0.75      bpca    0.1   0.1 0.75    5 20 msqrob
4052     0   1    3    0   0.15    0.6708204 none   0.75      bpca    0.1   0.1 0.75   10 20 msqrob
4053     0   1 3936    0 196.80  880.1163559 none   0.75 svdImpute    0.1   0.1 0.75    5 20   voom
4054     0   1 3938    0 393.65 1211.6292326 none   0.75 svdImpute    0.1   0.1 0.75   10 20   voom
4055     0   1 3938    0 787.00 1614.8904151 none   0.75       lls    0.1   0.1 0.75   10  5   voom
4056     0   1 3928    0 196.40  878.3275016 none   0.75       lls    0.1   0.1 0.75   10 10   voom

## save results:
write.table(tbl, file="perm_results.grp.tsv", sep="\t", quote=F, row.names=F)
```

---

## Running one step at a time:

Sometimes you might want to begin with R objects containing expression data,
feature metadata, and observation metadata, rather than reading that 
information from disk. You may also want to explicitly run individual steps 
of the workflow, rather than just calling `h0testr::run()`. These tasks are 
fairly straightforward, as long as the R objects are properly formatted as 
described and checked below: 

```
## need an expression matrix, feature metadata as a data.frame, 
##   and observation metadata as a data.frame:

> ls()
[1] "exprs" "feats"    "samps"

## expression matrix should be a numeric matrix:

> is.matrix(exprs)
[1] TRUE

> class(c(exprs))
[1] "numeric"

## rownames of expression matrix in sync with feature metadata; note that the 
##   matching column name in the feature metadata does not need to be 
##   called 'PrecursorId', and can be configured in the config list (see the
##   customization of the config list below):

> summary(rownames(exprs) == feats$PrecursorId)
   Mode    TRUE 
logical  200790 

## column names of expression matrix in sync with sample metadata; note that 
##   the matching column name in the obsvervation metadata does not need to 
##   be called 'ObservationId', and can be configured in the config list (see
##   the customization of the config list below):

> summary(colnames(exprs) == samps$ObservationId)
   Mode    TRUE 
logical      27 

## make required state list; element names must be 'expression', 'features', 
##   and 'samples':

state <- list(expression=exprs, features=feats, samples=samps)
rm(exprs, feats, samps)

## customize configuration list:

config <- h0testr::new_config()             ## default configuration
config$feat_id_col <- "PrecursorId"
config$gene_id_col <- "GeneGroup"
config$obs_id_col <- "ObservationId"
config$sample_id_col <- "SampleId"

config$frm <- ~age + sex + age:sex
config$test_term <- "age:sex"
config$permute_var <- ""
config$sample_factors <- list(
  age=c("4mo", "12mo", "24mo"),
  sex=c("female", "male")  
)

## a reasonable setup when replication is limited:

config$normalization_method <- "RLE"
config$impute_method <- "unif_global_lod"
config$impute_quantile <- 0
config$test_method <- "trend"

## run workflow step-by-step, instead of simply calling 
##   h0testr::run(state, config):

out <- h0testr::initialize(state, config)
out$state <- h0testr::add_filter_stats(out$state, out$config)
out$state <- h0testr::prefilter(out$state, out$config)
out$state <- h0testr::permute(out$state, out$config)
out <- h0testr::normalize(out$state, out$config)
out <- h0testr::combine_replicates(out$state, out$config)
out <- h0testr::combine_features(out$state, out$config)
out <- h0testr::filter(out$state, out$config)
out <- h0testr::impute(out$state, out$config, is_log_transformed=TRUE)
result <- h0testr::test(out$state, out$config)

## hits in format returned by underlying test method:
head(result$original)

## hits in standardized format for method comparison:
head(result$standard)

## most methods also return the fitted model:
summary(result$fit)
```


