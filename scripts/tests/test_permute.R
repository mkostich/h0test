## Tests for h0testr::permute(), which shuffles one observation covariate so a
##   run against the shuffled labels estimates the null.

usage <- function(msg=NULL) {
  if(!is.null(msg)) cat("ERROR:", msg, "\n\n", file=stderr())
  cat(
    "Test h0testr::permute(): the conditions under which it does nothing, that it",
    "moves only the named column of state$samples and leaves the expression",
    "matrix and feature table alone, that it shuffles at the sample level so",
    "technical replicates of one sample keep one label, that it is reproducible",
    "under a seed, how the variable argument resolves against config$permute_var,",
    "the errors it raises, and that it destroys the association it is meant to",
    "destroy.",
    "",
    "Usage: Rscript test_permute.R <r_dir>",
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
    "  Rscript test_permute.R C:/path/to/h0test/h0testr/R",
    "  Rscript test_permute.R ../../h0test/h0testr/R",
    "  Rscript test_permute.R C:/path/to/h0testr/R > test_permute.out 2>&1",
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

## two technical replicates per sample, so sample level and observation level
##   are distinguishable:

set.seed(101)
e <- sim1(n_obs=12, n_feats=20)$mat
sid <- rep(paste0("s", 1:6), each=2)
st <- list(
  expression=e,
  features=data.frame(feature_id=rownames(e), gene=rownames(e)),
  samples=data.frame(observation_id=colnames(e), sid=sid,
    grp=rep(c("ctl", "ctl", "ctl", "trt", "trt", "trt"), each=2),
    age=rep(c(4, 9, 14, 19, 24, 29), each=2),
    batch=rep(c("b1", "b2"), 6), stringsAsFactors=FALSE)
)
cfg <- list(feat_col="feature_id", obs_col="observation_id",
  sample_id_col="sid", log_file=log_file)

###############################################################################
section("permute() does nothing unless asked")

report(identical(permute(st, cfg), st),
  "neither the argument nor config$permute_var set leaves the state untouched")

m <- mark()
cfg_off <- cfg
cfg_off$permute_var <- ""
report(identical(permute(st, cfg_off), st),
  "config$permute_var set to the empty string leaves the state untouched")
report(logged("skipping permutation", since=m),
  "and says so, rather than passing silently")

report(identical(permute(st, cfg, variable=""), st),
  "the variable argument set to the empty string leaves the state untouched")

###############################################################################
section("only the named column moves")

cfg_grp <- cfg
cfg_grp$permute_var <- "grp"

set.seed(7)
p <- permute(st, cfg_grp)

report(identical(p$expression, st$expression),
  "permute() leaves the expression matrix byte for byte")
report(identical(p$features, st$features),
  "and leaves the feature table alone")
report(identical(dim(p$samples), dim(st$samples)) &&
  identical(names(p$samples), names(st$samples)),
  "and returns a samples table of the same shape")

for(nom in c("observation_id", "sid", "age", "batch")) {
  report(identical(p$samples[[nom]], st$samples[[nom]]),
    paste0("and leaves samples$", nom, " alone"))
}

report(identical(sort(p$samples$grp), sort(st$samples$grp)),
  "the permuted column holds the same labels, in some order")
report(identical(class(p$samples$grp), class(st$samples$grp)),
  "and keeps the type it had")

###############################################################################
section("the shuffle is at the sample level")

## a technical replicate is the same biological sample, so both replicates have
##   to carry the same permuted label or the design becomes incoherent:

one_per_sample <- function(s) {
  tapply(s$samples$grp, s$samples$sid, function(v) length(unique(v)))
}
report(all(one_per_sample(p) == 1),
  "both replicates of a sample carry one permuted label")

set.seed(3)
ok <- TRUE
for(k in 1:20) {
  ok <- ok && all(one_per_sample(permute(st, cfg_grp)) == 1)
}
report(ok, "and that holds over twenty shuffles")

## six samples, three labels each, so the sample level multiset is preserved:
tab0 <- table(st$samples$grp[!duplicated(st$samples$sid)])
tab1 <- table(p$samples$grp[!duplicated(p$samples$sid)])
report(identical(tab0, tab1), "the labels are dealt out to samples without replacement")

###############################################################################
section("reproducibility")

set.seed(31)
a <- permute(st, cfg_grp)
set.seed(31)
b <- permute(st, cfg_grp)
report(identical(a, b), "the same seed gives the same permutation")

orders <- character(0)
for(k in 1:20) {
  set.seed(k)
  orders <- c(orders, paste(permute(st, cfg_grp)$samples$grp, collapse=""))
}
report(length(unique(orders)) > 1, "and different seeds give different permutations")
report(any(orders != paste(st$samples$grp, collapse="")),
  "at least one of which differs from the labels it started with")

###############################################################################
section("the variable argument")

set.seed(5)
p_arg <- permute(st, cfg, variable="grp")
set.seed(5)
p_cfg <- permute(st, cfg_grp)
report(identical(p_arg, p_cfg),
  "the variable argument gives what config$permute_var gives")

cfg_age <- cfg
cfg_age$permute_var <- "age"
set.seed(5)
p_over <- permute(st, cfg_age, variable="grp")
report(identical(p_over$samples$age, st$samples$age),
  "the argument wins where both are set, leaving config$permute_var alone")
report(identical(p_over$samples$grp, p_arg$samples$grp),
  "and permutes the column the argument named")

set.seed(9)
p_num <- permute(st, cfg_age)
report(is.numeric(p_num$samples$age) &&
  identical(sort(p_num$samples$age), sort(st$samples$age)),
  "a numeric covariate permutes as readily as a factor one")

###############################################################################
section("errors")

m <- mark()
report(threw(permute(st, cfg, variable="no_such_column")),
  "a variable naming no column of samples is refused")
report(logged("no_such_column", since=m),
  "and the message names the variable that was asked for")

cfg_bad <- cfg
cfg_bad$sample_id_col <- "no_such_column"
report(threw(permute(st, cfg_bad, variable="grp")),
  "a config$sample_id_col naming no column of samples is refused")

st_bad <- st
st_bad$features <- st_bad$features[-1, , drop=FALSE]
report(threw(permute(st_bad, cfg_grp)),
  "permute() checks the state it is handed before returning it")

###############################################################################
section("permutation destroys the association it targets")

## the point of the exercise: a planted effect is detectable in the real
##   labels and not in permuted ones. Counted with an ordinary t-test rather
##   than through h0testr::test(), so this measures permute() and nothing else.

set.seed(202)
samps <- sim_samples(factors=list(grp=c("ctl", "trt")), n_per_cell=6)
sim <- sim_design(samps, frm=~grp, test_term="grp", n_genes=200,
  n_genes_signif=100, effects=3, p_drop=0, mcar_p=0)
st2 <- sim$state
cfg2 <- sim$config
cfg2$log_file <- log_file
cfg2$permute_var <- cfg2$test_term
## init_state() would set these from the id columns; sim_design() leaves them empty:
cfg2$feat_col <- cfg2$feat_id_col
cfg2$obs_col <- cfg2$obs_id_col

x <- log2(st2$expression + 1)

n_hits <- function(labels) {
  p <- apply(x, 1, function(v) {
    if(length(unique(labels[!is.na(v)])) < 2) return(NA_real_)
    tryCatch(stats::t.test(v ~ labels)$p.value, error=function(e) NA_real_)
  })
  sum(p < 0.05, na.rm=TRUE)
}

hits0 <- n_hits(st2$samples[[cfg2$test_term]])
report(hits0 > 50, "the planted effect is detectable in the real labels")

set.seed(11)
hits1 <- integer(0)
for(k in 1:10) {
  sp <- permute(st2, cfg2)
  hits1 <- c(hits1, n_hits(sp$samples[[cfg2$test_term]]))
}
report(max(hits1) < hits0 / 2,
  "and every one of ten permutations finds fewer than half as many")
report(stats::median(hits1) < 0.15 * nrow(x),
  "leaving a hit rate near what a null of this size should give")

###############################################################################

cat("\n## log file:", log_file, "\n")
cat("\n## passes:", n_pass, "; failures:", n_fail, "; elapsed:",
  round(as.numeric(difftime(Sys.time(), t0, units="secs"))), "s\n")

quit(save="no", status=if(n_fail > 0) 1 else 0)
