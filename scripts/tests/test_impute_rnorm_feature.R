usage <- function(msg=NULL) {
  if(!is.null(msg)) cat("ERROR:", msg, "\n\n", file=stderr())
  cat(
    "Test the scale handling of h0testr::impute_rnorm_feature(): how the",
    "is_log_transformed argument resolves against config$is_log_transformed, that a",
    "negative value is refused on the raw scale and ordinary on the log scale, that the",
    "dispersion floor is sqrt(mean) on the raw scale and a tenth of the whole matrix",
    "standard deviation on the log scale, how the scale. argument and config$impute_scale",
    "resolve, and that h0testr::impute() forwards the scale to it.",
    "",
    "Usage: Rscript test_impute_rnorm_feature.R <r_dir>",
    "",
    "Required positional arguments:",
    "  <r_dir>  Path to the h0testr package R/ source directory; all .R files",
    "             there are sourced (the installed package is not used).",
    "",
    "Output: one PASS/FAIL line per assertion to stdout, then a count of",
    "  passes and failures. Log output, including messages from expected",
    "  errors, is written to a temporary file whose path is reported at the end.",
    "",
    "Exit codes: 0 all assertions passed; 1 one or more failed; 2 usage error.",
    "",
    "Examples:",
    "  Rscript test_impute_rnorm_feature.R C:/path/to/h0test/h0testr/R",
    "  Rscript test_impute_rnorm_feature.R ../../h0test/h0testr/R",
    "  Rscript test_impute_rnorm_feature.R C:/path/to/h0testr/R > t.out 2>&1",
    sep="\n", file=stderr()
  )
  quit(status=2)
}

args <- commandArgs(trailingOnly=TRUE)
if(length(args) != 1) usage("wrong number of arguments")

r_dir <- args[1]
if(!dir.exists(r_dir)) usage(paste("r_dir not a directory:", r_dir))

for(f in list.files(r_dir, pattern="[.]R$", full.names=TRUE)) source(f)

###############################################################################

log_file <- tempfile(fileext=".log")
## logged() reads this before anything has written to it:
invisible(file.create(log_file))
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

threw <- function(expr) inherits(try(expr, silent=TRUE), "try-error")

mark <- function() length(readLines(log_file))

logged <- function(pattern, since=0) {
  txt <- readLines(log_file)
  if(since > 0) {
    if(length(txt) <= since) return(FALSE)
    txt <- txt[(since + 1):length(txt)]
  }
  any(grepl(pattern, txt, fixed=TRUE))
}

mk <- function(e) {
  if(is.null(rownames(e))) rownames(e) <- paste0("f", seq_len(nrow(e)))
  if(is.null(colnames(e))) colnames(e) <- paste0("o", seq_len(ncol(e)))
  list(expression=e, features=data.frame(feature_id=rownames(e)),
    samples=data.frame(observation_id=colnames(e)))
}

cfg_log <- list(is_log_transformed=TRUE, log_file=log_file)
cfg_raw <- list(is_log_transformed=FALSE, log_file=log_file)
cfg_bare <- list(log_file=log_file)

set.seed(101)
e_log <- sim1(n_obs=8, n_feats=60, mcar_p=0.15)$mat
e_log <- log2(e_log + 1)
e_log <- e_log[apply(e_log, 1, function(v) sum(!is.na(v)) >= 2), , drop=FALSE]
## centered, so the fixture carries values on both sides of zero:
e_log <- e_log - mean(e_log, na.rm=TRUE)
i_na <- is.na(e_log)

e_raw <- 2^(e_log - min(e_log, na.rm=TRUE) + 1) - 1

###############################################################################
section("state and scale resolution")

report(any(i_na) && any(!i_na), "fixture carries both measured and missing values")
report(any(e_log < 0, na.rm=TRUE), "and the log scale fixture carries negative values")

report(threw(impute_rnorm_feature(mk(e_log), cfg_bare)),
  "neither the argument nor the config key set is refused")

r <- try(impute_rnorm_feature(mk(e_log), cfg_bare, is_log_transformed=TRUE),
  silent=TRUE)
report(!inherits(r, "try-error"), "the argument alone resolves the scale")

r <- try(impute_rnorm_feature(mk(e_log), cfg_log), silent=TRUE)
report(!inherits(r, "try-error"), "the config key alone resolves the scale")

m <- mark()
report(threw(impute_rnorm_feature(mk(e_log), cfg_log, is_log_transformed=FALSE)),
  "an argument disagreeing with the config key is refused")
report(logged("cannot both be right", since=m), "saying they describe the same data")

report(threw(impute_rnorm_feature(list(expression=as.data.frame(e_log)), cfg_log)),
  "a state whose expression is not a matrix is refused")

e_all_na <- e_log
e_all_na[1, ] <- NA
report(threw(impute_rnorm_feature(mk(e_all_na), cfg_log)),
  "a feature with no measured value at all is refused")

###############################################################################
section("log scale")

set.seed(11)
r_log <- impute_rnorm_feature(mk(e_log), cfg_log)

report(!any(is.na(r_log$expression)), "every missing value is filled in")
report(all(r_log$expression[!i_na] == e_log[!i_na]), "measured values are left alone")
report(any(r_log$expression[i_na] < 0), "imputed values below zero are kept")
report(all(is.finite(r_log$expression)), "no imputed value is infinite")

m_row <- rowMeans(e_log, na.rm=TRUE)
d_log <- r_log$expression[i_na] - rep(m_row, ncol(e_log))[i_na]
report(abs(mean(d_log)) < 0.5, "imputed values sit around the feature mean")

###############################################################################
section("raw scale")

m <- mark()
report(threw(impute_rnorm_feature(mk(e_log), cfg_raw)),
  "a negative value on the raw scale is refused")
report(logged("contains negative values on the", since=m),
  "naming the values rather than the scale as the problem")

set.seed(11)
r_raw <- try(impute_rnorm_feature(mk(e_raw), cfg_raw), silent=TRUE)
report(!inherits(r_raw, "try-error"), "a non-negative raw matrix is accepted")

if(!inherits(r_raw, "try-error")) {
  report(!any(is.na(r_raw$expression)), "every missing value is filled in")
  report(all(r_raw$expression > 0), "and no imputed value is at or below zero")
}

e_zero <- e_raw
e_zero[!is.na(e_zero)][1] <- 0
report(threw(impute_rnorm_feature(mk(e_zero), cfg_raw)),
  "a zero on the raw scale is refused by the f.pos_mat check")

###############################################################################
section("dispersion floor")

## one measured value per feature, so every draw uses the floor:
n_f <- 60
n_o <- 12
set.seed(3)
obs1 <- stats::rnorm(n_f, mean=10, sd=2)
e_one <- matrix(NA_real_, nrow=n_f, ncol=n_o)
for(k in seq_len(n_f)) e_one[k, sample(n_o, 1)] <- obs1[k]
i_one <- is.na(e_one)

set.seed(5)
r_one <- impute_rnorm_feature(mk(e_one), cfg_log)
d_one <- (r_one$expression - rep(obs1, n_o))[i_one]
s_want <- 0.1 * stats::sd(c(e_one), na.rm=TRUE)

report(abs(stats::sd(d_one) / s_want - 1) < 0.2,
  "one measured value on the log scale draws at a tenth of the matrix spread")

obs2 <- exp(stats::rnorm(n_f, mean=log(1000), sd=0.3))
e_one_raw <- matrix(NA_real_, nrow=n_f, ncol=n_o)
for(k in seq_len(n_f)) e_one_raw[k, sample(n_o, 1)] <- obs2[k]
i_one_raw <- is.na(e_one_raw)

set.seed(5)
r_one_raw <- impute_rnorm_feature(mk(e_one_raw), cfg_raw)
d_one_raw <- ((r_one_raw$expression - rep(obs2, n_o)) /
  rep(sqrt(obs2), n_o))[i_one_raw]

report(abs(stats::sd(d_one_raw) - 1) < 0.2,
  "one measured value on the raw scale draws at sqrt(mean)")

set.seed(13)
e_flat <- matrix(rep(stats::rnorm(n_f, mean=10, sd=2), n_o), nrow=n_f)
e_flat[cbind(seq_len(n_f), sample(n_o, n_f, replace=TRUE))] <- NA
i_flat <- is.na(e_flat)

set.seed(5)
r_flat <- impute_rnorm_feature(mk(e_flat), cfg_log)
d_flat <- (r_flat$expression - rowMeans(e_flat, na.rm=TRUE))[i_flat]
s_flat <- 0.1 * stats::sd(c(e_flat), na.rm=TRUE)

report(stats::sd(d_flat) > 0.5 * s_flat,
  "a feature whose measured values have no spread still gets the floor")

###############################################################################
section("scale resolution")

set.seed(17)
e_wide <- matrix(stats::rnorm(n_f * n_o, mean=10, sd=3), nrow=n_f)
e_wide[sample(length(e_wide), 0.3 * length(e_wide))] <- NA
e_wide <- e_wide[apply(e_wide, 1, function(v) sum(!is.na(v)) >= 2), , drop=FALSE]
i_wide <- is.na(e_wide)
m_wide <- rowMeans(e_wide, na.rm=TRUE)

spread <- function(scale., cfg) {
  set.seed(23)
  r <- impute_rnorm_feature(mk(e_wide), cfg, scale.=scale.)
  stats::sd((r$expression - m_wide)[i_wide])
}

s1 <- spread(1, cfg_log)
s3 <- spread(3, cfg_log)
report(abs(s3 / s1 - 3) < 0.3, "the scale. argument multiplies the dispersion")

set.seed(23)
r_arg <- impute_rnorm_feature(mk(e_wide), cfg_log, scale.=3)
cfg_s3 <- cfg_log
cfg_s3$impute_scale <- 3
set.seed(23)
r_cfg <- impute_rnorm_feature(mk(e_wide), cfg_s3)
report(identical(r_arg$expression, r_cfg$expression),
  "config$impute_scale gives what the argument gives")

cfg_s001 <- cfg_log
cfg_s001$impute_scale <- 0.01
set.seed(23)
r_both <- impute_rnorm_feature(mk(e_wide), cfg_s001, scale.=3)
report(identical(r_both$expression, r_arg$expression),
  "and the argument wins where both are set")

###############################################################################
section("impute() forwards the scale")

cfg_i <- new_config()
cfg_i$log_file <- log_file
cfg_i$impute_method <- "rnorm_feature"
cfg_i$is_log_transformed <- TRUE
cfg_i$feat_id_col <- cfg_i$gene_id_col <- cfg_i$feat_col <- "feature_id"
cfg_i$obs_id_col <- cfg_i$sample_id_col <- cfg_i$obs_col <- "observation_id"

r_i <- try(suppressMessages(impute(mk(e_log), cfg_i)), silent=TRUE)
report(!inherits(r_i, "try-error"),
  "impute(method='rnorm_feature') runs on a log scale matrix with negative values")
if(!inherits(r_i, "try-error")) {
  e_i <- r_i$state$expression
  report(is.matrix(e_i) && identical(dim(e_i), dim(e_log)),
    "returning a matrix of the same shape")
  report(!any(is.na(e_i)), "and fills in every missing value")
  report(any(e_i[i_na] < 0), "keeping the imputed values below zero")
}

cfg_i$is_log_transformed <- FALSE
report(threw(suppressMessages(impute(mk(e_log), cfg_i))),
  "and the same call on the raw scale is refused, so the scale reaches the imputer")

###############################################################################

cat("\n## log file:", log_file, "\n")
cat("\n## passes:", n_pass, "; failures:", n_fail, "; elapsed:",
  round(as.numeric(difftime(Sys.time(), t0, units="secs"))), "s\n")

quit(save="no", status=if(n_fail > 0) 1 else 0)
