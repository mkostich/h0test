## Tests for handling of continuous (numeric) covariates in h0testr, upstream of
##   hypothesis testing: covariate classification, covariate value checks,
##   factor level setting, replicate combining, and the marginal pre-screen.

usage <- function(msg=NULL) {
  if(!is.null(msg)) cat("ERROR:", msg, "\n\n", file=stderr())
  cat(
    "Test upstream (pre-hypothesis-testing) handling of continuous covariates",
    "in h0testr: covariate classification, value checks, factor levels,",
    "combine_replicates(), and filter_features_by_formula().",
    "",
    "Usage: Rscript test_covariates.R <r_dir> [--data_dir=<dir>]",
    "",
    "Required positional arguments:",
    "  <r_dir>          Path to the h0testr package R/ source directory; all",
    "                     .R files there are sourced (the installed package is",
    "                     not used).",
    "",
    "Optional named arguments:",
    "  --data_dir=<dir> Path to a directory holding pg.expression.tsv,",
    "                     pg.features.tsv, and samples.tsv; enables an",
    "                     additional smoke test on real data. Skipped if",
    "                     absent.",
    "",
    "Output: one PASS/FAIL line per assertion to stdout, then a count of",
    "  passes and failures. Messages from expected errors are written to a",
    "  temporary log file, whose path is reported at the end.",
    "",
    "Exit codes: 0 all assertions passed; 1 one or more failed; 2 usage error.",
    "",
    "Examples:",
    "  Rscript test_covariates.R C:/path/to/h0test/h0testr/R",
    "  Rscript test_covariates.R ../../h0test/h0testr/R --data_dir=/data/seer2",
    "  Rscript test_covariates.R C:/path/to/h0testr/R > test_covariates.out 2>&1",
    sep="\n", file=stderr()
  )
  quit(status=2)
}

args <- commandArgs(trailingOnly=TRUE)
if(length(args) < 1 || length(args) > 2) usage("wrong number of arguments")

r_dir <- args[1]
if(!dir.exists(r_dir)) usage(paste("r_dir not a directory:", r_dir))

data_dir <- NULL
if(length(args) %in% 2) {
  if(!grepl("^--data_dir=", args[2])) usage(paste("unrecognized argument:", args[2]))
  data_dir <- sub("^--data_dir=", "", args[2])
  if(!dir.exists(data_dir)) usage(paste("data_dir not a directory:", data_dir))
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

## returns TRUE iff expr threw an error:
threw <- function(expr) inherits(try(expr, silent=TRUE), "try-error")

mk_state <- function(samps, n_feats=12, seed=101) {
  set.seed(seed)
  exprs <- matrix(stats::rnorm(n_feats * nrow(samps), 20, 2), nrow=n_feats)
  rownames(exprs) <- paste0("f", 1:n_feats)
  colnames(exprs) <- samps$observation_id
  feats <- data.frame(feature_id=rownames(exprs), gene_id=rownames(exprs))
  list(expression=exprs, features=feats, samples=samps)
}

mk_cfg <- function(frm, test_term, reference_levels=character(0), ...) {
  cfg <- list(
    obs_id_col="observation_id", sample_id_col="observation_id",
    feat_id_col="feature_id", gene_id_col="gene_id",
    frm=frm, test_term=test_term, reference_levels=reference_levels,
    log_file=log_file, save_state=FALSE
  )
  extra <- list(...)
  for(nom in names(extra)) cfg[[nom]] <- extra[[nom]]
  cfg
}

## 8 observations; sex is a declared factor, age and wt are continuous,
##   bead is an undeclared character variable, flag is logical:

samps0 <- data.frame(
  observation_id=paste0("o", 1:8),
  sex=rep(c("F", "M"), 4),
  age=c(21, 34, 45, 52, 63, 29, 38, 47),
  wt=c(22, 25, 31, 28, 19, 33, 27, 24),
  bead=rep(c("NPA", "NPB"), each=4),
  flag=rep(c(TRUE, FALSE), 4),
  stringsAsFactors=FALSE
)

###############################################################################
section("covariate classification")

state <- mk_state(samps0)

cfg <- mk_cfg(~sex, "sex", c(sex="F"))
out <- try(initialize(state, cfg, minimal=TRUE), silent=TRUE)
report(!inherits(out, "try-error"), "0 continuous covariates accepted")
report(identical(unname(out$config$covariate_types), "factor"), "sex classified factor")

cfg <- mk_cfg(~sex + age, "age", c(sex="F"))
out <- try(initialize(state, cfg, minimal=TRUE), silent=TRUE)
report(!inherits(out, "try-error"), "1 continuous covariate accepted")
types <- out$config$covariate_types
report(identical(types[["age"]], "numeric") && identical(types[["sex"]], "factor"),
  "types: age numeric, sex factor")
report(is.numeric(out$state$samples$age), "continuous covariate left numeric")

cfg <- mk_cfg(~age + wt, "age", character(0))
out <- try(initialize(state, cfg, minimal=TRUE), silent=TRUE)
report(!inherits(out, "try-error"), "2 continuous covariates, no factors, accepted")

cfg <- mk_cfg(~age * wt, "age:wt", character(0))
out <- try(initialize(state, cfg, minimal=TRUE), silent=TRUE)
report(!inherits(out, "try-error"), "numeric:numeric interaction accepted")

cfg <- mk_cfg(~age * sex, "age:sex", c(sex="M"))
out <- try(initialize(state, cfg, minimal=TRUE), silent=TRUE)
report(!inherits(out, "try-error"), "numeric:factor interaction accepted")
report(identical(levels(out$state$samples$sex), c("M", "F")),
  "declared non-alphabetical reference level honored")

cfg <- mk_cfg(~sex + bead, "sex", c(sex="F"))
report(threw(initialize(state, cfg, minimal=TRUE)),
  "undeclared character covariate rejected")

cfg <- mk_cfg(~sex + flag, "sex", c(sex="F"))
out <- try(initialize(state, cfg, minimal=TRUE), silent=TRUE)
report(!inherits(out, "try-error"), "undeclared logical covariate accepted")
report(identical(levels(out$state$samples$flag), c("FALSE", "TRUE")),
  "logical covariate levels FALSE, TRUE")

samps <- samps0
samps$dose <- rep(c(1, 2), 4)
cfg <- mk_cfg(~sex + dose, "dose", c(sex="F", dose="2"))
out <- try(initialize(mk_state(samps), cfg, minimal=TRUE), silent=TRUE)
report(!inherits(out, "try-error"), "numeric covariate declared in reference_levels accepted")
report(identical(levels(out$state$samples$dose), c("2", "1")),
  "declared numeric covariate becomes factor with declared reference")

## a stale cached classification must not win over the current state:
cfg <- mk_cfg(~sex + age, "age", c(sex="F"))
cfg$covariate_types <- c(sex="numeric", age="factor")
out <- try(initialize(state, cfg, minimal=TRUE), silent=TRUE)
report(!inherits(out, "try-error") &&
  identical(out$config$covariate_types[["age"]], "numeric"),
  "stale config$covariate_types discarded by initialize()")

###############################################################################
section("reference levels and reported levels")

## only the reference level is declared; the remaining levels are sorted:
samps <- samps0
samps$grp <- rep(c("ctl", "trt", "xtra", "ctl"), 2)
cfg <- mk_cfg(~grp, "grp", c(grp="trt"))
out <- try(initialize(mk_state(samps), cfg, minimal=TRUE), silent=TRUE)
report(!inherits(out, "try-error"), "3-level factor with declared reference accepted")
report(identical(levels(out$state$samples$grp), c("trt", "ctl", "xtra")),
  "declared reference level first, remaining levels sorted")
report(identical(out$config$factor_levels[["grp"]], c("trt", "ctl", "xtra")),
  "config$factor_levels holds resolved level ordering")
report(setequal(names(out$config$factor_levels), "grp"),
  "config$factor_levels covers exactly the factor covariates")

## the resolved levels are also logged, so they are visible in a run log:
cfg <- mk_cfg(~grp, "grp", c(grp="trt"))
n0 <- length(readLines(log_file))
out <- try(initialize(mk_state(samps), cfg, minimal=TRUE), silent=TRUE)
txt <- readLines(log_file)
txt <- txt[(n0 + 1):length(txt)]
report(any(grepl("covariate grp : factor; levels: trt ctl xtra ; reference: trt",
  txt, fixed=TRUE)), "resolved levels and reference level logged")

## continuous covariates have no levels to report:
cfg <- mk_cfg(~sex + age, "age", c(sex="F"))
out <- try(initialize(state, cfg, minimal=TRUE), silent=TRUE)
report(setequal(names(out$config$factor_levels), "sex"),
  "continuous covariate absent from config$factor_levels")

## a typo in the one level the user does type is still caught:
cfg <- mk_cfg(~grp, "grp", c(grp="Trt"))
report(threw(initialize(mk_state(samps), cfg, minimal=TRUE)),
  "reference level absent from the data rejected")

## but a level present in the data need not be declared:
cfg <- mk_cfg(~grp, "grp", c(grp="ctl"))
out <- try(initialize(mk_state(samps), cfg, minimal=TRUE), silent=TRUE)
report(!inherits(out, "try-error") &&
  identical(levels(out$state$samples$grp), c("ctl", "trt", "xtra")),
  "undeclared non-reference levels accepted")

## a declared variable that is not in the formula is a misconfiguration:
cfg <- mk_cfg(~age, "age", c(sex="F"))
report(threw(initialize(state, cfg, minimal=TRUE)),
  "reference_levels entry absent from config$frm rejected")

###############################################################################
section("covariate value checks")

samps <- samps0
samps$age[3] <- NA
cfg <- mk_cfg(~sex + age, "age", c(sex="F"))
report(threw(initialize(mk_state(samps), cfg, minimal=TRUE)),
  "NA in continuous covariate rejected")

samps <- samps0
samps$sex[3] <- NA
cfg <- mk_cfg(~sex + age, "age", c(sex="F"))
report(threw(initialize(mk_state(samps), cfg, minimal=TRUE)),
  "NA in factor covariate rejected")

samps <- samps0
samps$age[3] <- Inf
cfg <- mk_cfg(~sex + age, "age", c(sex="F"))
report(threw(initialize(mk_state(samps), cfg, minimal=TRUE)),
  "non-finite continuous covariate rejected")

samps <- samps0
samps$age <- 40
cfg <- mk_cfg(~sex + age, "sex", c(sex="F"))
report(threw(initialize(mk_state(samps), cfg, minimal=TRUE)),
  "constant continuous covariate rejected")

samps <- samps0
samps$sex <- "F"
cfg <- mk_cfg(~sex + age, "age", c(sex="F"))
report(threw(initialize(mk_state(samps), cfg, minimal=TRUE)),
  "constant factor covariate rejected")

## warning, not error, for a continuous covariate with few distinct values:
samps <- samps0
samps$age <- rep(c(4, 12), 4)
cfg <- mk_cfg(~sex + age, "age", c(sex="F"))
n0 <- length(readLines(log_file))
out <- try(initialize(mk_state(samps), cfg, minimal=TRUE), silent=TRUE)
txt <- readLines(log_file)
txt <- txt[(n0 + 1):length(txt)]
report(!inherits(out, "try-error"),
  "continuous covariate with 2 distinct values accepted")
report(any(grepl("WARNING: numeric", txt)), "few-distinct-values warning logged")

## and the warning respects the configured cutoff:
cfg <- mk_cfg(~sex + age, "age", c(sex="F"), n_distinct_numeric_warn=1)
n0 <- length(readLines(log_file))
out <- try(initialize(mk_state(samps), cfg, minimal=TRUE), silent=TRUE)
txt <- readLines(log_file)
txt <- txt[(n0 + 1):length(txt)]
report(!any(grepl("WARNING: numeric", txt)),
  "no warning when n_distinct_numeric_warn below n distinct values")

###############################################################################
section("config parameter checks")

cfg <- new_config()
cfg$log_file <- log_file
report(isTRUE(try(check_config(cfg), silent=TRUE)), "new_config() passes check_config()")
report(identical(cfg$n_distinct_numeric_warn, 5), "n_distinct_numeric_warn default is 5")

report(identical(cfg$reference_levels, c(age="young", gender="Male")),
  "new_config() reference_levels is a named character vector")

cfg$covariate_types <- c(age="numeric", gender="factor")
report(isTRUE(try(check_config(cfg), silent=TRUE)), "covariate_types accepted by check_config()")

cfg$factor_levels <- list(age=c("young", "old"), gender=c("Male", "Female"))
report(isTRUE(try(check_config(cfg), silent=TRUE)), "factor_levels accepted by check_config()")

cfg2 <- cfg
cfg2$reference_levels <- c("young", "Male")
report(threw(check_config(cfg2)), "unnamed reference_levels rejected")

cfg2$reference_levels <- c(age="young", age="old")
report(threw(check_config(cfg2)), "duplicated names in reference_levels rejected")

cfg2$reference_levels <- c(age="young", gender=NA)
report(threw(check_config(cfg2)), "NA value in reference_levels rejected")

cfg2$reference_levels <- c(age="young", gender="")
report(threw(check_config(cfg2)), "empty value in reference_levels rejected")

cfg2$reference_levels <- list(age="young", gender="Male")
report(threw(check_config(cfg2)), "list-valued reference_levels rejected")

cfg2 <- cfg
cfg2$factor_levels <- list(c("young", "old"))
report(threw(check_config(cfg2)), "unnamed factor_levels rejected")

cfg2 <- cfg
cfg2$sample_factors <- list(age=c("young", "old"))
report(threw(check_config(cfg2)), "obsolete sample_factors parameter rejected")

cfg$n_distinct_numeric_warn <- "5"
report(threw(check_config(cfg)), "character n_distinct_numeric_warn rejected")

cfg$n_distinct_numeric_warn <- -1
report(threw(check_config(cfg)), "negative n_distinct_numeric_warn rejected")

###############################################################################
section("combine_replicates()")

## two observations per sample; covariates constant within sample:
samps <- data.frame(
  observation_id=paste0("o", 1:8),
  sample_id=rep(paste0("s", 1:4), each=2),
  sex=rep(c("F", "M"), each=4),
  age=rep(c(21, 34, 45, 52), each=2),
  stringsAsFactors=FALSE
)
cfg <- mk_cfg(~sex + age, "age", c(sex="F"), sample_id_col="sample_id")
out <- try(initialize(mk_state(samps), cfg, minimal=TRUE), silent=TRUE)
out2 <- try(combine_replicates(out$state, out$config, fn=sum), silent=TRUE)
report(!inherits(out2, "try-error"), "combine_replicates() ok with constant covariates")
report(!inherits(out2, "try-error") && ncol(out2$state$expression) %in% 4,
  "combine_replicates() collapses 8 observations to 4 samples")

## continuous covariate varying within a sample:
samps2 <- samps
samps2$age[2] <- 22
cfg <- mk_cfg(~sex + age, "age", c(sex="F"), sample_id_col="sample_id")
out <- try(initialize(mk_state(samps2), cfg, minimal=TRUE), silent=TRUE)
report(threw(combine_replicates(out$state, out$config, fn=sum)),
  "continuous covariate varying within sample rejected")

## factor covariate varying within a sample:
samps3 <- samps
samps3$sex[2] <- "M"
cfg <- mk_cfg(~sex + age, "age", c(sex="F"), sample_id_col="sample_id")
out <- try(initialize(mk_state(samps3), cfg, minimal=TRUE), silent=TRUE)
report(threw(combine_replicates(out$state, out$config, fn=sum)),
  "factor covariate varying within sample rejected")

###############################################################################
section("filter_features_by_formula()")

## continuous covariate with all-unique values must not wipe out every feature;
##   here every feature is fully observed, so all should be kept:
state <- mk_state(samps0)
cfg <- mk_cfg(~sex + age, "age", c(sex="F"))
out <- try(initialize(state, cfg, minimal=TRUE), silent=TRUE)
st <- try(filter_features_by_formula(out$state, out$config), silent=TRUE)
report(!inherits(st, "try-error") &&
  nrow(st$expression) %in% nrow(out$state$expression),
  "complete features kept under continuous covariate with unique values")

## feature with only one non-NA value is dropped; feature with three is kept.
##   Formula has the continuous covariate only, so that the factor screen
##   cannot be what drops the feature:
cfg_age <- mk_cfg(~age, "age", character(0))
out_age <- try(initialize(mk_state(samps0), cfg_age, minimal=TRUE), silent=TRUE)
st0 <- out_age$state
st0$expression[1, 2:8] <- NA
st0$expression[2, 4:8] <- NA
st <- try(filter_features_by_formula(st0, out_age$config), silent=TRUE)
report(!inherits(st, "try-error") && !("f1" %in% rownames(st$expression)),
  "feature with 1 non-NA value dropped (continuous covariate)")
report(!inherits(st, "try-error") && "f2" %in% rownames(st$expression),
  "feature with 3 non-NA values kept (continuous covariate)")

## feature observed at a single value of the continuous covariate is dropped:
samps <- samps0
samps$age <- rep(c(10, 20, 30, 40), 2)
state <- mk_state(samps)
cfg <- mk_cfg(~sex + age, "age", c(sex="F"))
out <- try(initialize(state, cfg, minimal=TRUE), silent=TRUE)
st0 <- out$state
st0$expression[3, !(samps$age %in% 10)] <- NA
st <- try(filter_features_by_formula(st0, out$config), silent=TRUE)
report(!inherits(st, "try-error") && !("f3" %in% rownames(st$expression)),
  "feature observed at a single covariate value dropped")

## factor behavior unchanged: feature observed in only one level is dropped:
state <- mk_state(samps0)
cfg <- mk_cfg(~sex, "sex", c(sex="F"))
out <- try(initialize(state, cfg, minimal=TRUE), silent=TRUE)
st0 <- out$state
st0$expression[4, samps0$sex %in% "M"] <- NA
st <- try(filter_features_by_formula(st0, out$config), silent=TRUE)
report(!inherits(st, "try-error") && !("f4" %in% rownames(st$expression)),
  "feature observed in only one factor level dropped")
report(!inherits(st, "try-error") && nrow(st$expression) %in%
  (nrow(st0$expression) - 1), "other features kept under factor covariate")

## a variable appearing only inside an interaction is screened too:
state <- mk_state(samps0)
cfg <- mk_cfg(~sex + sex:age, "sex:age", c(sex="F"))
out <- try(initialize(state, cfg, minimal=TRUE), silent=TRUE)
st0 <- out$state
st0$expression[5, 2:8] <- NA
st <- try(filter_features_by_formula(st0, out$config), silent=TRUE)
report(!inherits(st, "try-error") && !("f5" %in% rownames(st$expression)),
  "variable appearing only in an interaction is screened")

## interaction of two continuous covariates:
cfg <- mk_cfg(~age * wt, "age:wt", character(0))
out <- try(initialize(mk_state(samps0), cfg, minimal=TRUE), silent=TRUE)
st <- try(filter_features_by_formula(out$state, out$config), silent=TRUE)
report(!inherits(st, "try-error") &&
  nrow(st$expression) %in% nrow(out$state$expression),
  "numeric:numeric interaction: complete features kept")

###############################################################################
## optional smoke test on real data:

if(!is.null(data_dir)) {

  section("real data smoke test")

  exprs <- as.matrix(utils::read.table(paste0(data_dir, "/pg.expression.tsv"),
    header=TRUE, sep="\t", quote="", as.is=TRUE))
  feats <- utils::read.table(paste0(data_dir, "/pg.features.tsv"),
    header=TRUE, sep="\t", quote="", as.is=TRUE)
  samps <- utils::read.table(paste0(data_dir, "/samples.tsv"),
    header=TRUE, sep="\t", quote="", as.is=TRUE)

  obs0 <- c(
    "s448.16.A_NPB_001", "s448.16.A_NPA_001", "s448.12.A_NPA_002",
    "Seer.PC1.Control_NA_021", "Seer.DC.Control_NA_022",
    "Seer.MPE.Control_NA_023", "Seer.PC2.Control_NA_024",
    "MS.Pep.Control.extended_NA_01", "MS.Pep.Control_NA_01"
  )
  i <- samps$obs %in% obs0
  samps <- samps[!i, , drop=FALSE]
  exprs <- exprs[, !i, drop=FALSE]

  ## synthetic continuous covariate, constant within sample:
  set.seed(101)
  sids <- sort(unique(samps$sample_id))
  age <- stats::setNames(round(stats::runif(length(sids), 20, 70), 1), sids)
  samps$age <- age[samps$sample_id]

  state <- list(expression=exprs, features=feats, samples=samps)

  cfg <- new_config()
  cfg$feat_id_col <- cfg$gene_id_col <- "protein_group"
  cfg$obs_id_col <- "obs"
  cfg$sample_id_col <- "sample_id"
  cfg$dir_out <- tempdir()
  cfg$log_file <- log_file
  cfg$save_state <- FALSE
  cfg$frm <- ~sex + age
  cfg$test_term <- "age"
  cfg$permute_var <- ""
  cfg$reference_levels <- c(sex="F")
  cfg$normalization_method <- "RLE"
  cfg$impute_method <- "unif_sample_lod"
  cfg$impute_quantile <- 0
  cfg$test_method <- "trend"

  out <- try(initialize(state, cfg), silent=TRUE)
  report(!inherits(out, "try-error"), "real data: initialize() with continuous covariate")

  if(!inherits(out, "try-error")) {
    report(identical(out$config$covariate_types[["age"]], "numeric"),
      "real data: age classified numeric")
    out$state <- add_filter_stats(out$state, out$config)
    out$state <- prefilter(out$state, out$config)
    out <- normalize(out$state, out$config)
    out <- combine_replicates(out$state, out$config, fn=sum)
    out <- combine_features(out$state, out$config)
    n_before <- nrow(out$state$expression)
    st <- try(filter_features_by_formula(out$state, out$config), silent=TRUE)
    report(!inherits(st, "try-error"), "real data: filter_features_by_formula() runs")
    report(!inherits(st, "try-error") && nrow(st$expression) > 0.5 * n_before,
      paste0("real data: pre-screen keeps most features (",
        if(inherits(st, "try-error")) NA else nrow(st$expression), " of ",
        n_before, ")"))
  }

  ## an undeclared character covariate in the real metadata is an error:
  cfg$frm <- ~sex + bead
  cfg$test_term <- "sex"
  report(threw(initialize(state, cfg)), "real data: undeclared bead rejected")

  ## run holds values like 001 and 001.2, so read.table() makes it numeric, and
  ##   it is therefore treated as continuous; it has 27 distinct values, so the
  ##   default cutoff does not warn about it:
  cfg$frm <- ~sex + run
  cfg$test_term <- "sex"
  n0 <- length(readLines(log_file))
  out <- try(initialize(state, cfg), silent=TRUE)
  txt <- readLines(log_file)
  txt <- txt[(n0 + 1):length(txt)]
  report(!inherits(out, "try-error") &&
    identical(out$config$covariate_types[["run"]], "numeric"),
    "real data: numeric-looking run classified continuous")
  report(!any(grepl("WARNING: numeric", txt)),
    "real data: no warning for run at default cutoff (27 distinct values)")

  ## with the cutoff raised above its number of distinct values, it warns:
  cfg$n_distinct_numeric_warn <- 30
  n0 <- length(readLines(log_file))
  out <- try(initialize(state, cfg), silent=TRUE)
  txt <- readLines(log_file)
  txt <- txt[(n0 + 1):length(txt)]
  report(!inherits(out, "try-error") && any(grepl("WARNING: numeric", txt)),
    "real data: run warns as possible miscoded factor at raised cutoff")
}

###############################################################################

cat("\n## passes:", n_pass, "; failures:", n_fail, "; elapsed:",
  round(as.numeric(difftime(Sys.time(), t0, units="secs"))), "s\n")
cat("## log file:", log_file, "\n")

quit(save="no", status=if(n_fail > 0) 1 else 0)
