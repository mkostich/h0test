# h0test.?.R
## Hypothesis testing for ms-based proteomics

---

## Quick start:

```
## bash: 

## install under ~/opt/h0test:
mkdir -p ~/opt
cd ~/opt
git clone https://github.com/mkostich/h0test

## make a working directory, just for this:
mkdir test1
cd test1

## copy h0test.1.R to working directory, and configure for this experiment:
rsync ~/opt/h0test/h0test.1.R ./
vi h0test.1.R

## invoke R instance with required packages (limma and edgeR):
container='/projects/compsci/jgeorge/kostim/resources/containers/msproteomics_r.20250325a.sif'
apptainer exec $container R

## in R:
rm(list=ls())
source("./h0test.1.R")
source("~/opt/h0test/h0test.2.R")
source("~/opt/h0test/h0test.3.R")
q()
```
  
---

## Dependencies

- R, with limma and edgeR packages
- Developed and tested with R 4.4.2, limma 3.62.2, and edgeR 4.4.2
- Dependencies available in this apptainer container:<br/>
  `/projects/compsci/jgeorge/kostim/resources/containers/msproteomics_r.20250325a.sif`
- Should work on Linux, Mac, or Windows

---

## Process

1. Prefilter for all `NA` or all `0` rows (features) or columns (observations)
2. Normalize, and if not part of normalization, transform (default: `log2`)
3. Combine technical replicates (take median)
4. Filter features and samples based on `number of values > 0`
5. Impute missing values
6. Test hypotheses, etc.

---
  
## Script summary

- h0test.0.R: collection of generic functions;
- h0test.1.R: entrypoint; configure this;
    read data -> configure sample metadata -> prefilter -> permute (if requested)
- h0test.2.R: normalize -> transform (e.g. log2) -> combine technical replicates
- h0test.3.R: filter data -> impute missing values -> h0 testing
- h0tune.1.R: used in place of h0test.1.R when tuning parameters (beta)
---

## Inputs

Inputs are three tab-delimited text files in directory `DIR_IN`. The names of 
the files can be configured using parameters `FEATURE_FILE_IN`, 
`SAMPLE_FILE_IN`, and `DATA_FILE_IN`:

- `features.tsv`: feature metadata, with unique key `FEAT_ID_COL`
- `samples.tsv`: observation metadata, with:
  - `OBS_COL`: unique key; technical replicate identifier
  - `SAMPLE_COL`: biosample identifier
  - need separate columns `OBS_COL` and `SAMPLE_COL`, even if no technical replicates
  - all variables referred to in the formula specified by configuration variable `FRM`
- `expression.tsv`: LFQ (e.g. MaxLFQ) normalized numeric expression data; 
  - feature rows; observation columns;
  - `rownames(exprs) == features[, FEAT_ID_COL]`
  - `colnames(exprs) == samples[, OBS_COL]`
  - `DIR_IN`, `FEAT_ID_COL` and `OBS_COL` configured in `h0test.1.R`
---

## Configuration

Runs are configured by editing the h0test.1.R file:

```
DIR_IN <- "."                         ## where DATA_FILE_IN, FETURE_FILE_IN, and SAMPLE_FILE_IN live
FEATURE_FILE_IN <- "features.tsv"     ## feature annotation .tsv; row features
SAMPLE_FILE_IN <- "samples.tsv"       ## sample annotation .tsv; row observations
DATA_FILE_IN <- "expression.tsv"      ## .tsv lfq normalized quant matrix; row features, column observations
DIR_OUT="."                           ## output directory

## formula for testing: actual formula can have '+' and ':'; not tested w/ e.g. '*' yet.
FRM <- ~ age + strain + gender + age:strain   ## formula with variable of interest and covariates
TEST_TERM <- "age:strain"                     ## term in FRM on which test is to be performed
PERMUTE_VAR <- ""                             ## variable to permute; "" for no permutation (normal execution)
## PERMUTE_VAR <- "age"                       ## e.g. permute the 'age' variable prior to h0 testing

## set levels of factor variables in samples; 
##   works for character and logical variables; 
##   also works for integer or numeric variables, if you want to treat as categorical;
##   first level of each factor treated as reference level:
SAMPLE_FACTORS <- list(
  age=c("4", "12", "24"),                     ## if don't factorize age, treated as continuous numeric
  strain=c("B6", "A_J", "Balbc_J", "CAST", "NZO"),
  gender=c("Male", "Female")
)  

## samples and features column names (new cols are introduced by the code):
N_SAMPLES_EXPR_COL <- "n_samps_expr"          ## new col; n samples expressing feature
MEDIAN_RAW_COL <- "median_raw"                ## new col; median feature expression in expressing samples
N_FEATURES_EXPR_COL <- "n_feats_expr"         ## new col; n features express
FEAT_ID_COL <- "accession"                    ## features[, FEAT_ID_COL] == rownames(exprs)
OBS_COL <- "assay_id"                         ## unique id for observations; samples[, OBS_COL] == colnames(exprs)
SAMPLE_COL <- "sample_id"                     ## id for samples in samples; not unique if tech reps (>1 obs/sample)

## tunable options: defaults are usually ok, except:
##   for dia: usually works ok: RLE:unif_sample_lod:0.05 for NORM_METHOD:IMPUTE_METHOD:IMPUTE_QUANTILE
##   for dda: usually works ok: quantile:0.75:unif_sample_lod:0 for NORM_METHOD:NORM_QUANTILE:IMPUTE_METHOD:IMPUTE_QUANTILE
NORM_METHOD <- "quantile"            ## in c("vsn","cpm","quantile","qquantile","TMM","TMMwsp","RLE","upperquartile")
NORM_QUANTILE <- 0.75                ## for quantile normalization; 0.5 is median; 0.75 is upper quartile;
TRANSFORM_METHOD <- "log2"           ## in c('log2', 'log10', 'none')
N_SAMPLES_MIN <- 2                   ## min samples/feature w/ feature expression > 0 to keep feature
N_FEATURES_MIN <- 1000               ## min features/sample w/ expression > 0 to keep sample
IMPUTE_METHOD <- "unif_sample_lod"   ## in c("unif_global_lod", "unif_sample_lod", "sample_lod", "rnorm_feature")
IMPUTE_QUANTILE <- 0                 ## quantile for unif_* imputation methods
IMPUTE_SCALE <- 1                    ## for rnorm_feature, adjustment on sd of distribution [1: no change];
TEST_METHOD <- "trend"               ## in c("voom", "trend")

## output file naming:
LOG_FILE <- "log.txt"                ## log file path; or "" for log to console
FEATURE_MID_OUT <- ".features"       ## midfix for output feature metadata files
SAMPLE_MID_OUT <- ".samples"         ## midfix for output sample metadata file
DATA_MID_OUT <- ".expression"        ## midfix for output expression files
RESULT_MID_OUT <- ".results"         ## prefix for output results file
SUFFIX_OUT <- ".tsv"                 ## suffix for output files

## dependency:
SRC_TEST <- "/path/to/h0test/h0test.0.R"

## misc:
PROBS <- c(0, 0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99, 1.0)
WIDTH <- 110
```

---

## Outputs

Writes a series of intermediate files and final results file to output directory 
specificed by `DIR_OUT`. Each set of intermediate files consists of an 
expression matrix file, a feature metadata file, and a sample metadata file. 
The names of the files can be customized using `DATA_MID_OUT`, `FEATURES_MID_OUT`, 
and `SAMPLE_MID_OUT`, respectively. The common file suffix can be customized by
changing `SUFFIX_OUT`.

With the default settings for `DATA_MID_OUT`, `FEATURES_MID_OUT`,  `SAMPLE_MID_OUT`, 
and `SUFFIX_OUT`, output files in `DIR_OUT` will be:

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

The `results.tsv` with `TEST_TERM` set to `age` looks something like this:

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

## where ... represents additional feature columns from FEATURE_FILE_IN.
```

---

## NORM_METHOD

Unless otherwise noted, values returned by each `NORM_METHOD` are on the original scale:

**vsn**: variance stabilizing transformation, using `limma::normalizeVSN()`; returned values on log2 scale

**cpm**: counts per million; for each sample: `multiplier * (intensities / sum(intensities, na.rm=T))`

**quantile**: for each sample: `multiplier * intensities / quantile(intensities, probs=NORM_QUANTILE, na.rm=T))`

**qquantile**: classical quantile normalization, using `limma::normalizeQuantiles()`<br/>

**TMM**: Trimmed mean of medians normalization, using `edgeR::calcNormFactors()` then `edgeR::cpm()`. 
  See [https://doi.org/10.1186/gb-2010-11-3-r25](https://doi.org/10.1186/gb-2010-11-3-r25 "Robinson and Oshlack, 2010")
  
**TMMwsp**: TMM with singleton pairing, using `edgeR::calcNormFactors()` then `edgeR::cpm()`. 
  May be better than TMM when many zeros.

**RLE**: Relative log expression, using `edgeR::calcNormFactors()` then `edgeR::cpm()`. 
  See [https://doi.org/10.1038/npre.2010.4282.1](https://doi.org/10.1038/npre.2010.4282.1 "Anders and Huber, 2010")

**upperquartile**: upper quartile expression, using `edgeR::calcNormFactors()` then `edgeR::cpm()`.
  See [https://doi.org/10.1186/1471-2105-11-94](https://doi.org/10.1186/1471-2105-11-94 "Bullard et al., 2010")


---

## TRANSFORM_METHOD

**log2**: `log2(x+1)`; assumes `all(values >= 0)`

**log10**: `log10(x+1)`; assumes `all(values >= 0)`

**none**: no transformation; pick this if normalization results already on log scale (e.g. for vsn).

---

## IMPUTE_METHOD

Typically, we replace all zero values with `NA` prior to imputing.

**unif_global_lod**: for each feature remaining after filtering, calculate the lowest observed value.
  Calculate the `IMPUTE_QUANTILE` quantile `q` of lowest observed values across all features. Replace 
  each `NA` value by randomly sampling between `0` and `q`: `runif(1, min=0, max=q)`.

**unif_sample_lod**: for each sample, calculate the `IMPUTE_QUANTILE` quantile `q` of observed
  values. Replace each `NA` value by randomly sampling between `0` and `q`: `runif(1, min=0, max=q)`.

**sample_lod**: for each sample, replace all missing values with the minimum observed value in 
  that sample.

**rnorm_feature**: for each feature, calculate the mean `m` and standard deviation `s`. Replace 
  each `NA` value with a random draw from a normal distribution centered at `m` with 
  standard deviation  `IMPUTE_SCALE * s`: `rnorm(1, mean=m, sd=IMPUTE_SCALE*s)`.

---

## TEST_METHOD

**voom**: `limma::voom() -> limma::lmFit() -> limma::eBayes(trend=F) -> limma::topTable()`

**trend**: `limma::lmFit() -> limma::eBayes(trend=T) -> limma::topTable()`

---

## Tuning (beta)

Run one instance without permutation, and `N` (`N >= 10` is recommended)
instances with permutation. Permute by setting `PERMUTE_VAR` to a variable 
(not a multivariate term, like an interaction, which will fail) that is
part of the `TEST_TERM`.

Run each instance in its own subdirectory. Here we pretend we have 
subdirectories numbered from `0` to `10`, where the unpermuted run 
will be in subdirectory `0`, and permuted runs are in subdirectories
`1` thru `10`. In the main directory, do something like the following
to get the tuning script and configure it for this run. Note that 
`h0tune.1.R` differs from `h0test.1.R` in that some configuration
parameters are commented out, and are meant to be set in the calling
environment before the scripts are sourced:

```
## bash: 

## copy h0test.1.R to working directory:
rsync /path/to/h0test/h0test.1.R ./

## configure script for this experiment:
vi h0tune.1.R            

## make subdirectories 0 thru 10:
mkdir $(seq 0 10)
```

Perform each run in a separate subdirectory:

```
## change to a subdirectory where run will be performed:
cd 0

## invoke R instance with required packages (limma and edgeR):
container='/projects/compsci/jgeorge/kostim/resources/containers/msproteomics_r.20250325a.sif'
apptainer exec $container R
```

In R, for each run:

```
## R:

rm(list=ls())

script1 <- "./h0tune.1.R"
script2 <- "~/opt/h0test/h0test.2.R"
script3 <- "~/opt/h0test/h0test.3.R"

out_file <- "0.age_strain.tune.tsv"
TEST_TERM <- "age:strain"
PERMUTE_VAR <- ""         ## to run on permuted data, set PERMUTE_VAR to e.g. "age" or "strain"
## PERMUTE_VAR <- "age"   ## e.g. to run on data where 'age' variable has been randomized/permuted

set.seed((Sys.getpid() + round((as.numeric(Sys.time()) * 100))) %% 1000)
source(script1)           ## configure vars; load data; prefilter; permute if req'd; prep exprs1, samps1, feats1; 
out <- NULL
IMPUTE_SCALE <- 1
NORM_QUANTILE <- 0.75
for(NORM_METHOD in c("vsn", "cpm", "q50", "q75", "qquantile", "TMM", "TMMwsp", "RLE", "upperquartile")) {
  if(NORM_METHOD %in% "q50") {
    NORM_METHOD <- "quantile"
    NORM_QUANTILE <- 0.5
  } else if(NORM_METHOD %in% c("q75")) {
    NORM_METHOD <- "quantile"
    NORM_QUANTILE <- 0.75
  } else if(NORM_METHOD %in% "upperquartile") NORM_QUANTILE <- 0.75
  source(script2)     ## normalize/transform; combine tech reps; exprs2, samps2, feats2 from exprs1, samps1, feats1
  for(TEST_METHOD in c("voom", "trend")) {
    for(IMPUTE_METHOD in c("unif_global_lod", "unif_sample_lod", "sample_lod", "rnorm_feature")) {
      if(IMPUTE_METHOD %in% c("unif_global_lod", "unif_sample_lod")) {
        for(IMPUTE_QUANTILE in c(0, 0.01, 0.05, 0.1, 0.25)) {
          source(script3)  ## filter -> impute -> test -> tbl; from exprs2, samps2, feats2
          out_i <- data.frame(norm=NORM_METHOD, norm_quant=NORM_QUANTILE, impute=IMPUTE_METHOD, imp_quant=IMPUTE_QUANTILE, 
            scale=IMPUTE_SCALE, test=TEST_METHOD, perm=PERMUTE_VAR, nhits=sum(tbl$adj.P.Val < 0.05), ntrys=nrow(tbl), 
            time=format(Sys.time(), "%H:%M:%S"), stringsAsFactors=F)
          out <- rbind(out, out_i)
          print(out); flush.console()
        }
      } else if(IMPUTE_METHOD %in% c("rnorm_feature")) {
        for(IMPUTE_SCALE in c(1, 0.5, 0.25, 0.1)) {
          source(script3)  ## filter -> impute -> test -> tbl; from exprs2, samps2, feats2
          out_i <- data.frame(norm=NORM_METHOD, norm_quant=NORM_QUANTILE, impute=IMPUTE_METHOD, imp_quant=IMPUTE_QUANTILE, 
            scale=IMPUTE_SCALE, test=TEST_METHOD, perm=PERMUTE_VAR, nhits=sum(tbl$adj.P.Val < 0.05), ntrys=nrow(tbl), 
            time=format(Sys.time(), "%H:%M:%S"), stringsAsFactors=F)
          out <- rbind(out, out_i)
          print(out); flush.console()
        }
      } else if(IMPUTE_METHOD %in% c("sample_lod")) {
        source(script3)  ## filter -> impute -> test -> tbl; from exprs2, samps2, feats2
        out_i <- data.frame(norm=NORM_METHOD, norm_quant=NORM_QUANTILE, impute=IMPUTE_METHOD, imp_quant=IMPUTE_QUANTILE, 
          scale=IMPUTE_SCALE, test=TEST_METHOD, perm=PERMUTE_VAR, nhits=sum(tbl$adj.P.Val < 0.05), ntrys=nrow(tbl), 
          time=format(Sys.time(), "%H:%M:%S"), stringsAsFactors=F)
        out <- rbind(out, out_i)
        print(out); flush.console()
      }
    }
  }
}
write.table(out, file=out_file, sep="\t", quote=F, row.names=F)
```

The results of each run look like this:

```
         norm norm_quant          impute imp_quant scale  test perm nhits ntrys     time
1         vsn       0.50 unif_global_lod      0.00  1.00  voom  age     0  6822 12:06:31
2         vsn       0.50 unif_global_lod      0.01  1.00  voom  age     0  6822 12:06:41
3         vsn       0.50 unif_global_lod      0.05  1.00  voom  age     1  6822 12:06:50
4         vsn       0.50 unif_global_lod      0.10  1.00  voom  age     0  6822 12:06:59
5         vsn       0.50 unif_global_lod      0.25  1.00  voom  age     1  6822 12:07:09
6         vsn       0.50 unif_sample_lod      0.00  1.00  voom  age     0  6822 12:07:19
7         vsn       0.50 unif_sample_lod      0.01  1.00  voom  age     4  6822 12:07:28
8         vsn       0.50 unif_sample_lod      0.05  1.00  voom  age    26  6822 12:07:37
9         vsn       0.50 unif_sample_lod      0.10  1.00  voom  age    48  6822 12:07:47
10        vsn       0.50 unif_sample_lod      0.25  1.00  voom  age    71  6822 12:07:56
```

Collect results. From parent directory, in R (only base packages required for this part):

```
rm(list=ls())
## main unpermuted results:
dat1 <- read.table("0/0.age_strain.tune.tsv", header=T, sep="\t", quote="", as.is=T)

## ten sets of permutation results:
obj <- list()
for(idx in 1:10) {
  f <- paste(idx, paste0(idx, ".age_strain.tune.tsv"), sep="/")
  cat("reading", f, "\n"); flush.console()
  dat_i <- read.table(f, header=T, sep="\t", quote="", as.is=T)
  dat_i$run <- idx
  obj[[idx]] <- dat_i
}
dat0 <- do.call(rbind, obj)  ## rbind the permutation results
rm(obj)

## keys with norm, norm_quant, imput, imp_quant, scale, and test:
k0 <- apply(dat0[, 1:6], 1, paste, collapse=":")
k1 <- apply(dat1[, 1:6], 1, paste, collapse=":")

## permuted results in dat0; take max, median, mean, and sd of 10 permutation results:
perm_max <- tapply(dat0$nhits, k0, max, na.rm=T)
perm_mid <- tapply(dat0$nhits, k0, median, na.rm=T)
perm_avg <- tapply(dat0$nhits, k0, mean, na.rm=T)
perm_sd <- tapply(dat0$nhits, k0, sd, na.rm=T)

## get them in the same order as dat1 (k1 made from dat1):
dat1$max0 <- perm_max[k1]
dat1$mid0 <- perm_mid[k1]
dat1$avg0 <- perm_avg[k1]
dat1$sd0 <- perm_sd[k1]
dat1$perm <- NULL

## average number of false positives == average number of hits across 10 sets of permutation results;
##   false positive rate: (average number of false positives) / (total number of positives)

if(any(dat1$nhits %in% 0)) stop("any(dat1$nhits %in% 0)")
dat1$fdr <- dat1$avg / dat1$nhits
dat1 <- dat1[order(dat1$nhits, -dat1$fdr, decreasing=T), c("nhits", "fdr", "max0", "mid0", "avg0", "sd0", 
  "norm", "norm_quant", "impute", "imp_quant", "scale", "test")]
rownames(dat1) <- NULL

## dat1 has a bunch of statistics that can be used to evaluate option combinations;
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

> head(x)
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

