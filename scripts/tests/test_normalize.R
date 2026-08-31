## Tests for h0testr::normalize() and the normalize_*() functions behind it:
##   what happens when normalization has nothing to work with (an observation
##   with no measured value; no feature measured in every observation, which is
##   where edgeR's RLE takes its reference from), validation of the n_pts
##   argument of normalize_vsn(), and agreement between the documentation and
##   the code. Not a test of what each method computes; that is each method's
##   own business.

usage <- function(msg=NULL) {
  if(!is.null(msg)) cat("ERROR:", msg, "\n\n", file=stderr())
  cat(
    "Test h0testr::normalize() and the normalize_*() functions behind it:",
    "refusal and reporting when normalization has nothing to work with, the",
    "n_pts argument of h0testr::normalize_vsn(), and agreement between the",
    "documentation and the code.",
    "",
    "Usage: Rscript test_normalize.R <r_dir>",
    "",
    "Required positional arguments:",
    "  <r_dir>  Path to the h0testr package R/ source directory; all .R files",
    "             there are sourced (the installed package is not used). The",
    "             man/ directory beside it and readme.md above it are read for",
    "             the documentation assertions, and skipped if not found.",
    "",
    "Output: one PASS/FAIL line per assertion to stdout, then a count of",
    "  passes and failures. Log output, including messages from expected",
    "  errors, is written to a temporary file whose path is reported at the end.",
    "",
    "Exit codes: 0 all assertions passed; 1 one or more failed; 2 usage error.",
    "",
    "Examples:",
    "  Rscript test_normalize.R C:/path/to/h0test/h0testr/R",
    "  Rscript test_normalize.R ../../h0test/h0testr/R",
    "  Rscript test_normalize.R C:/path/to/h0testr/R > test_normalize.out 2>&1",
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

## TRUE iff the log holds the text, which is where f.err() and f.msg() put
##   their detail; since= skips lines written before a mark(), so that a message
##   left by an earlier assertion is not read as evidence for a later one:

mark <- function() {
  if(!file.exists(log_file)) return(0)     ## nothing has written to it yet
  length(readLines(log_file))
}

logged <- function(pattern, since=0) {
  if(!file.exists(log_file)) return(FALSE)
  txt <- readLines(log_file)
  if(since > 0) {
    if(length(txt) <= since) return(FALSE)
    txt <- txt[(since + 1):length(txt)]
  }
  any(grepl(pattern, txt, fixed=TRUE))
}

## the value of an expression, or NULL if it threw:
val <- function(expr) {
  out <- try(expr, silent=TRUE)
  if(inherits(out, "try-error")) NULL else out
}

###############################################################################
## shared data; sim1() draws strictly positive values with NA for dropouts, so
##   mnar_c0=-Inf gives a matrix with nothing missing:

set.seed(101)
m0 <- sim1(n_obs=6, n_feats=60)$mat
m0[m0 %in% 0] <- NA
mfull <- sim1(n_obs=6, n_feats=60, mnar_c0=-Inf)$mat

feats <- data.frame(pep=rownames(m0), gene=rownames(m0))
samps <- data.frame(obs=colnames(m0), grp=rep(c("ctl", "trt"), 3))

mk_state <- function(e) list(expression=e, features=feats, samples=samps)

cfg0 <- list(
  obs_id_col="obs", sample_id_col="obs", feat_id_col="pep", gene_id_col="gene",
  frm=~grp, test_term="grp", reference_levels=c(grp="trt"), log_file=log_file,
  feat_col="pep", obs_col="obs", save_state=FALSE
)

## an observation that measured nothing:
me <- m0
me[, 3] <- NA

## no feature measured in every observation: one missing value per feature,
##   spread around the columns so that every observation still measures plenty:
mz <- mfull
for(i in 1:nrow(mz)) mz[i, ((i - 1) %% ncol(mz)) + 1] <- NA

## the methods that cannot proceed on an observation measuring nothing, and the
##   ones that carry it through as the all-NA column it already was:
strict <- c("TMM", "TMMwsp", "RLE", "upperquartile", "qquantile")
tolerant <- c("cpm", "quantile", "div.mean", "div.median", "sum", "max",
  "loess", "vsn", "log2")

## and the one that refuses an empty observation for a different reason: it
##   refuses any missing value at all, whatever the pattern. Its own section
##   further down; here it is only kept out of the two lists above:

no_missing <- c("quantiles.robust")

report(sum(apply(me, 2, function(v) all(is.na(v)))) == 1,
  "the test matrix has exactly one observation measuring nothing")
report(sum(apply(mz, 1, function(v) !anyNA(v))) == 0,
  "and the other has no feature measured in every observation")
report(sum(apply(m0, 1, function(v) !anyNA(v))) > 0,
  "while the reference matrix has some of both")
report(identical(sort(c(strict, tolerant, no_missing, "none")),
  sort(normalize_methods())),
  "every method normalize_methods() names is covered by one list or the other")

###############################################################################
section("an observation measuring nothing: the methods that cannot proceed")

for(meth in strict) {
  i <- mark()
  report(threw(normalize(mk_state(me), cfg0, method=meth)),
    paste0(meth, ": an observation measuring nothing is refused"))
  report(logged("cannot normalize", i) && logged("no measured value", i),
    paste0("  ", meth, ": the refusal says what is wrong"))
  report(logged("observations: obs_3", i) && logged("filter_observations", i),
    paste0("  ", meth, ": and names the observation and what to do"))
  report(!logged("library sizes", i) && !logged("lengths differ", i) &&
      !logged("TRUE/FALSE needed", i),
    paste0("  ", meth, ": and no third-party message is what stopped it"))
}

###############################################################################
section("an observation measuring nothing: the methods that carry it through")

for(meth in tolerant) {
  i <- mark()
  q <- if(meth %in% "quantile") 0.75 else NULL
  out <- val(normalize(mk_state(me), cfg0, method=meth,
    normalization_quantile=q))
  report(!is.null(out) && logged("WARNING", i) &&
      logged("no measured value", i),
    paste0(meth, ": still normalizes, and says what it was handed"))
  e <- out$state$expression
  report(!is.null(out) && all(is.na(e[, 3])) && sum(is.nan(e)) == 0 &&
      sum(is.infinite(e)) == 0,
    paste0("  ", meth, ": the empty observation stays empty, no NaN, no Inf"))
}

i <- mark()
report(!is.null(val(normalize(mk_state(me), cfg0, method="none"))) &&
    !logged("no measured value", i),
  "none: normalizes nothing, so makes no claim about the data")

###############################################################################
section("a matrix with no such observation is not warned about")

for(meth in c("RLE", "upperquartile", "cpm", "quantile", "loess", "none")) {
  i <- mark()
  q <- if(meth %in% "quantile") 0.75 else NULL
  out <- val(normalize(mk_state(m0), cfg0, method=meth,
    normalization_quantile=q))
  report(!is.null(out) && !logged("no measured value", i),
    paste0(meth, ": runs without a word about missing observations"))
}

###############################################################################
section("RLE needs a feature measured in every observation")

i <- mark()
report(threw(normalize(mk_state(mz), cfg0, method="RLE")),
  "RLE is refused when no feature is measured in every observation")
report(logged("cannot be used on this matrix", i) &&
    logged("no feature is measured in every observation: 0 of 60", i),
  "  and the refusal names the condition and counts it")
report(logged("takes its reference", i),
  "  and says where RLE's reference comes from")
report(logged("library sizes should be finite and non-negative", i),
  "  and names the edgeR message it stands in for")
report(logged("'TMM', 'TMMwsp' and 'upperquartile'", i),
  "  and names the methods that do not need such a feature")

for(meth in c("TMM", "TMMwsp", "upperquartile", "cpm", "loess")) {
  out <- val(normalize(mk_state(mz), cfg0, method=meth))
  e <- out$state$expression
  report(!is.null(out) && sum(is.nan(e)) == 0 && sum(is.infinite(e)) == 0,
    paste0("  ", meth, ": normalizes that matrix anyway, no NaN, no Inf"))
}

m1 <- mz
m1[1, ] <- mfull[1, ]
report(!threw(normalize(mk_state(m1), cfg0, method="RLE")),
  "  RLE works again as soon as one feature is measured everywhere")

## NA is the only indicator of a missing value, so the count is on NA alone;
##   the zeros edgeR excludes are the ones normalize_edger() writes in itself:

m1[1, 2] <- 0
i <- mark()
threw(normalize_edger(mk_state(m1), cfg0, method="RLE"))
report(!logged("no feature is measured in every observation", i),
  "  a zero in an otherwise complete feature does not trip the count")

i <- mark()
report(threw(normalize_edger(mk_state(mz), cfg0, method="RLE")) &&
    logged("takes its reference", i),
  "a direct call to normalize_edger() is refused the same way")
i <- mark()
report(threw(normalize_edger(mk_state(me), cfg0, method="TMM")) &&
    logged("cannot normalize", i),
  "and so is a direct call handed an observation measuring nothing")
i <- mark()
report(threw(normalize_qquantile(mk_state(me), cfg0)) &&
    logged("cannot normalize", i),
  "the same for a direct call to normalize_qquantile()")
report(!threw(normalize_qquantile(mk_state(m0), cfg0)),
  "  which is otherwise unaffected")

###############################################################################
section("normalize_vsn() checks n_pts before vsn2 sees it")

for(np in list("42", c(3, 4), NA, NA_integer_, TRUE, NULL)) {
  i <- mark()
  report(threw(normalize_vsn(mk_state(mfull), cfg0, n_pts=np)) &&
      logged("has to be a single number", i),
    paste0("n_pts of class ", class(np), ", length ", length(np),
      ", is refused as not a single number"))
  report(!logged("missing value where", i) && !logged("length = 2", i),
    "  and no raw R message leaks out")
}

for(np in list(-5, 0, 3.7)) {
  i <- mark()
  report(threw(normalize_vsn(mk_state(mfull), cfg0, n_pts=np)) &&
      logged("whole number of 1 or more", i),
    paste0("n_pts of ", np, " is refused as not a whole number of 1 or more"))
}

n_max <- round(sqrt(nrow(mfull)))
i <- mark()
out <- val(normalize_vsn(mk_state(mfull), cfg0))
report(!is.null(out) && logged("n_pts > n_pts_max", i),
  "the default n_pts is still reduced to n_pts_max, with a message")
report(!is.null(out) && isTRUE(all.equal(out$expression,
  limma::normalizeVSN(mfull, minDataPointsPerStratum=n_max))),
  "  and the result is still limma's, at the reduced n_pts")
i <- mark()
out <- val(normalize_vsn(mk_state(mfull), cfg0, n_pts=3))
report(!is.null(out) && !logged("n_pts > n_pts_max", i) &&
    isTRUE(all.equal(out$expression,
      limma::normalizeVSN(mfull, minDataPointsPerStratum=3))),
  "a valid n_pts below n_pts_max is passed through untouched")
report(!threw(normalize_vsn(mk_state(mfull), cfg0, n_pts=8L)),
  "an integer n_pts is accepted")

###############################################################################
section("normalize_loess() checks method, and refuses affy/pairs on gaps")

for(meth in list("bogus", c("fast", "affy"), NA, NA_character_, 42, NULL)) {
  i <- mark()
  report(threw(normalize_loess(mk_state(mfull), cfg0, method=meth)) &&
      logged("method has to be one of fast affy pairs", i),
    paste0("loess method of class ", class(meth), ", length ", length(meth),
      ", is refused"))
  report(!logged("missing value where", i) && !logged("length = 2", i) &&
      !logged("should be one of", i),
    "  and no raw R or limma message leaks out")
}

## the refusal is not decoration: limma's affy and pairs normalize each pair of
##   observations on their difference, so a feature the pair does not both
##   measure comes back missing from both. Measured here rather than asserted,
##   so that this stops holding if limma ever stops doing it:

n_na0 <- sum(is.na(m0))
n_na_affy <- sum(is.na(limma::normalizeCyclicLoess(m0, span=0.7, method="affy")))

report(n_na_affy > n_na0,
  paste0("limma's affy really does spread missingness: ", n_na0, " missing in, ",
    n_na_affy, " out"))

for(meth in c("affy", "pairs")) {
  i <- mark()
  report(threw(normalize_loess(mk_state(m0), cfg0, method=meth)) &&
      logged("cannot be used on data with", i),
    paste0("loess method ", meth, " is refused on data with missing values"))
  report(logged(paste("missing values in state$expression:", n_na0), i),
    "  and the refusal counts them")
  report(logged("use method='fast'", i),
    "  and names the method that does not have the problem")
  report(!threw(normalize_loess(mk_state(mfull), cfg0, method=meth)),
    paste0("  while ", meth, " is accepted on a matrix with nothing missing"))
}

out <- val(normalize_loess(mk_state(m0), cfg0, method="fast"))
report(!is.null(out) && sum(is.na(out$expression)) == n_na0,
  "the default method 'fast' takes the gaps and leaves them where they were")

###############################################################################
section("loess is log transformed before the fit, not after it")

## cyclic loess is an additive correction, so on raw intensities it drives the
##   smallest measurements negative, and log2(x + 1) of anything at or below -1
##   is NaN. Measured on the same matrix the function is handed, so that the
##   assertions below are about a real hazard rather than a hypothetical one:

raw_fit <- limma::normalizeCyclicLoess(m0, span=0.7, method="fast")
n_nan_old <- sum(is.nan(log2(raw_fit + 1)))

report(n_nan_old > 0,
  paste0("normalizing before transforming would delete measurements: ",
    n_nan_old, " NaN"))

i <- mark()
out <- val(normalize(mk_state(m0), cfg0, method="loess"))

report(!is.null(out) && logged("transforming data before the loess fit", i),
  "normalize() says it transformed the data before the loess fit")
report(!is.null(out) && !any(is.nan(out$state$expression)),
  "and there is no NaN in the result")
report(!is.null(out) && sum(is.na(out$state$expression)) == n_na0,
  "and no measurement was lost: the gaps are the ones it started with")
report(!is.null(out) && isTRUE(all.equal(out$state$expression,
  limma::normalizeCyclicLoess(log2(m0 + 1), adaptive.span=TRUE,
    method="fast"))),
  "and the result is limma's fit of the transformed matrix")

i <- mark()
out <- val(normalize(mk_state(m0), cfg0, method="quantile"))
report(!is.null(out) && logged("transforming data", i) &&
    !logged("before the loess fit", i),
  "while the other methods are still transformed after normalizing")

###############################################################################
section("a single feature keeps its shape, or is refused by name")

m1f <- mfull[1, , drop=FALSE]           ## one feature, nothing missing
m1f_na <- m1f
m1f_na[1, 2] <- NA                       ## with one feature, a gap empties an obs

## the shared mk_state() carries 60 features, so a one-feature matrix needs its
##   own state, or f.check_state() refuses it before the method is reached:

mk1 <- function(e) {
  list(expression=e, features=data.frame(pep=rownames(e), gene=rownames(e)),
    samples=samps)
}

## the hazard: apply() drops the dim attribute of a single-row matrix, which is
##   how normalize_quantile(), normalize_cpm() and normalize_mscoreutils() each
##   lost the shape of a one-feature matrix. Measured rather than assumed, so
##   this notices if apply() ever stops doing it:

report(is.null(dim(apply(m1f, 2, function(v) v / 2))) &&
    !is.null(dim(apply(mfull, 2, function(v) v / 2))),
  "apply(m, 2, f) drops the dim of a 1-row matrix, but not of a 60-row one")

shape_ok <- function(x, m) {
  is.matrix(x) && identical(dim(x), dim(m)) && identical(dimnames(x), dimnames(m))
}

## through normalize(), where f.check_state() would notice a lost shape:

for(meth in c("quantile", "cpm", "sum", "max", "div.mean", "div.median",
    "quantiles.robust", "loess", "log2", "none")) {
  out <- val(normalize(mk1(m1f), cfg0, method=meth))
  report(!is.null(out) && shape_ok(out$state$expression, m1f),
    paste0("normalize(method='", meth, "') returns a 1 x 6 matrix, dimnames and all"))
}

## and by a direct call, where nothing would: normalize_quantile() returned a
##   named vector of one value per observation, and the two mscoreutils paths
##   handed a vector to a third party that reported "'x' must be an array of at
##   least two dimensions" or "dim(X) must have a positive length" instead:

report(shape_ok(val(normalize_quantile(mk1(m1f), cfg0))$expression, m1f),
  "normalize_quantile() called directly keeps the shape too")
report(shape_ok(val(normalize_cpm(mk1(m1f), cfg0))$expression, m1f),
  "normalize_cpm() called directly keeps the shape too")

for(meth in c("sum", "max", "div.mean", "div.median", "quantiles.robust")) {
  out <- val(normalize_mscoreutils(mk1(m1f), cfg0, method=meth))
  report(!is.null(out) && shape_ok(out$expression, m1f),
    paste0("normalize_mscoreutils(method='", meth, "') keeps it as well"))
}

## the shape fix must not have changed any number: same arithmetic as the
##   apply() version, on the 60-feature matrix with gaps:

f_q <- function(v) 1e3 * v / stats::quantile(v, probs=0.75, na.rm=TRUE)
report(isTRUE(all.equal(normalize_quantile(mk_state(m0), cfg0)$expression,
  apply(m0, 2, f_q), tolerance=0)),
  "normalize_quantile() returns exactly what the apply() version returned")

f_c <- function(v) 1e6 * (v / sum(v, na.rm=TRUE))
report(isTRUE(all.equal(normalize_cpm(mk_state(m0), cfg0)$expression,
  apply(m0, 2, f_c), tolerance=0)),
  "normalize_cpm() likewise")

for(meth in c("sum", "max", "div.mean", "div.median")) {
  old <- MsCoreUtils::normalize_matrix(apply(m0, 2, as.numeric), method=meth)
  rownames(old) <- rownames(m0)
  new <- val(normalize_mscoreutils(mk_state(m0), cfg0, method=meth))
  report(!is.null(new) && isTRUE(all.equal(new$expression, old, tolerance=0)),
    paste0("normalize_mscoreutils(method='", meth, "') likewise"))
}

## quantiles.robust is refused on m0, which has gaps, so the same comparison
##   runs on the complete matrix:

old <- MsCoreUtils::normalize_matrix(apply(mfull, 2, as.numeric),
  method="quantiles.robust")
rownames(old) <- rownames(mfull)
new <- val(normalize_mscoreutils(mk_state(mfull), cfg0,
  method="quantiles.robust"))
report(!is.null(new) && isTRUE(all.equal(new$expression, old, tolerance=0)),
  "normalize_mscoreutils(method='quantiles.robust') likewise, on mfull")

## the as.numeric() that the apply() call used to do still has to happen:

mi <- m0
mi[is.na(mi)] <- 0L
storage.mode(mi) <- "integer"
out <- val(normalize_mscoreutils(mk_state(mi), cfg0, method="div.mean"))
report(!is.null(out) && storage.mode(out$expression) %in% "double",
  "integer input still reaches MsCoreUtils as doubles")

## where one feature cannot be normalized at all, it is refused by name rather
##   than by whatever the third party has to say about it:

report(threw(limma::normalizeQuantiles(m1f)),
  "limma cannot quantile normalize a single feature")

i <- mark()
report(threw(normalize_qquantile(mk1(m1f), cfg0)) &&
    logged("normalize_qquantile: needs at least 2 features", i),
  "  so normalize_qquantile() refuses one")
report(logged("nrow(state$expression): 1", i),
  "  and says how many features there were")
report(logged("need at least two non-NA values to interpolate", i),
  "  and quotes what limma would have said")
report(!is.null(val(normalize_qquantile(mk1(mfull[1:2, , drop=FALSE]), cfg0))),
  "  while two features are enough for it")

## vsn is refused on one feature whatever the values, because what it does with
##   one depends on them: 12 of the first 20 single rows of mfull come back (all
##   six observations equal, so every difference between them erased) and the
##   other 8 fail with "L-BFGS-B needs finite values of 'fn'":

n_vsn_ok <- 0
for(k in 1:20) {
  if(!is.null(val(limma::normalizeVSN(mfull[k, , drop=FALSE],
    minDataPointsPerStratum=1)))) n_vsn_ok <- n_vsn_ok + 1
}
report(n_vsn_ok > 0 && n_vsn_ok < 20,
  paste0("vsn takes some single features and not others: ", n_vsn_ok, " of 20"))

v <- val(limma::normalizeVSN(mfull[1:20, , drop=FALSE][1, , drop=FALSE],
  minDataPointsPerStratum=1))
report(is.null(v) || length(unique(as.vector(v))) == 1,
  "  and where it takes one, every observation comes back with the same value")

i <- mark()
for(m in list(m1f, m1f_na)) {
  report(threw(normalize_vsn(mk1(m), cfg0)) &&
      logged("normalize_vsn: needs at least 2 features", i),
    paste0("so normalize_vsn() refuses a single feature (", sum(is.na(m)),
      " missing)"))
}
report(logged("nrow(state$expression): 1", i),
  "  and says how many features there were")
report(!is.null(val(normalize_vsn(mk_state(mz), cfg0))),
  "  while a 60-feature matrix with a gap in every observation goes through")

###############################################################################
section("normalize() names the method it refuses, before touching the data")

## the old message named config$normalization_method rather than the value it
##   had rejected, so it could name a method that is perfectly good and then
##   list that same method among the allowed values:

cfg_rle <- cfg0
cfg_rle$normalization_method <- "RLE"

i <- mark()
report(threw(normalize(mk_state(mfull), cfg_rle, method="bogus")) &&
    logged("method: bogus", i),
  "normalize() names the method it rejected, not config$normalization_method")
report(!logged("normalization_method: RLE", i),
  "  so a good config$normalization_method is not reported as the bad value")
report(logged("method has to be one of", i) && logged("quantiles.robust", i),
  "  and the refusal lists the methods normalize_methods() returns")

## and the name is checked before anything happens to the data: with
##   config$is_log_transformed TRUE the report used to be the log-scale error,
##   naming the unknown method as the thing that would transform the data twice:

cfg_lt <- cfg_rle
cfg_lt$is_log_transformed <- TRUE
i <- mark()
report(threw(normalize(mk_state(mfull), cfg_lt, method="bogus")) &&
    logged("method has to be one of", i) &&
    !logged("is_log_transformed is TRUE", i),
  "an unknown method is refused before the is_log_transformed check")

## every shape of bad method argument, none of them leaking a raw R message:

for(meth in list("bogus", c("RLE", "TMM"), NA, NA_character_, 42, TRUE)) {
  i <- mark()
  report(threw(normalize(mk_state(mfull), cfg_rle, method=meth)) &&
      logged("method has to be one of", i),
    paste0("normalize() method of class ", class(meth), ", length ",
      length(meth), ", is refused"))
  report(!logged("missing value where", i) && !logged("length = 2", i) &&
      !logged("condition has length", i),
    "  and no raw R message leaks out")
}

## NULL, character(0) and "" all mean unset, so they fall back to config and are
##   an error only when config does not name a method either:

for(meth in list(NULL, character(0), "")) {
  report(!threw(normalize(mk_state(mfull), cfg_rle, method=meth)),
    paste0("method of class ", class(meth), ", length ", length(meth),
      ", falls back to config$normalization_method"))
  i <- mark()
  report(threw(normalize(mk_state(mfull), cfg0, method=meth)) &&
      logged("both unset", i),
    "  and is 'both unset' when config does not name one either")
}

## with the check up front, normalize()'s final else is unreachable, which holds
##   only while every name normalize_methods() returns has a branch:

i <- mark()
n_ok <- 0
for(meth in normalize_methods()) {
  if(!is.null(val(normalize(mk_state(mfull), cfg0, method=meth)))) {
    n_ok <- n_ok + 1
  }
}
report(n_ok == length(normalize_methods()),
  paste0("every name normalize_methods() returns has a branch: ", n_ok, " of ",
    length(normalize_methods())))
report(!logged("no branch for method", i),
  "  so the guard left in the final else never fires")

###############################################################################
section("an argument that is not a single usable value is refused by name")

## no raw R message from any of these: each of the four below is what one of
##   these arguments used to produce when it reached the arithmetic or the
##   comparison unchecked:

no_raw <- function(i) {
  !logged("missing value where", i) && !logged("length = 2", i) &&
    !logged("condition has length", i) && !logged("non-numeric argument", i)
}

## a quantile or a span: two values, an NA of either class, a character and a
##   logical:

bad_prop <- list(c(0.5, 0.75), NA, NA_real_, "0.75", TRUE)

for(x in bad_prop) {
  i <- mark()
  report(threw(normalize_quantile(mk_state(mfull), cfg0,
        normalization_quantile=x)) &&
      logged("normalize_quantile: normalization_quantile has to be a single", i),
    paste0("normalize_quantile() normalization_quantile of class ", class(x),
      ", length ", length(x), ", is refused"))
  report(no_raw(i), "  and no raw R message leaks out")

  i <- mark()
  report(threw(normalize_edger(mk_state(mfull), cfg0, method="upperquartile",
        normalization_quantile=x)) &&
      logged("normalize_edger: normalization_quantile has to be a single", i),
    paste0("normalize_edger() normalization_quantile of class ", class(x),
      ", length ", length(x), ", is refused"))
  report(no_raw(i), "  and no raw R message leaks out")

  i <- mark()
  report(threw(normalize_loess(mk_state(mfull), cfg0, span=x)) &&
      logged("normalize_loess: span has to be a single number", i),
    paste0("normalize_loess() span of class ", class(x), ", length ",
      length(x), ", is refused"))
  report(no_raw(i), "  and no raw R message leaks out")
}

## the range checks that used to be reached with these values are still there:

for(x in c(-0.5, 1.5)) {
  i <- mark()
  report(threw(normalize_quantile(mk_state(mfull), cfg0,
        normalization_quantile=x)) && logged("< 0 || normalization_quantile > 1", i),
    paste0("a single normalization_quantile of ", x, " is still out of range"))
  i <- mark()
  report(threw(normalize_loess(mk_state(mfull), cfg0, span=x)) &&
      logged("span < 0 || span > 1", i),
    paste0("a single span of ", x, " is still out of range"))
}

## multiplier had no check at all, and it is the one of these arguments that
##   nothing downstream refuses either. Measured rather than asserted, so that
##   this stops holding if R ever stops recycling:

f_bad <- function(v) c(1e3, 1e6) * v / stats::quantile(v, probs=0.75, na.rm=TRUE)
m_bad <- mfull
for(j in seq_len(ncol(m_bad))) m_bad[, j] <- f_bad(mfull[, j])
sum_ok <- sum(normalize_quantile(mk_state(mfull), cfg0)$expression, na.rm=TRUE)

report(sum(m_bad, na.rm=TRUE) > 10 * sum_ok,
  paste0("a length-2 multiplier recycles into a different result, not an error: ",
    signif(sum_ok, 6), " vs ", signif(sum(m_bad, na.rm=TRUE), 6)))

bad_mult <- list(c(1e3, 1e6), NA, NA_real_, "1e3", TRUE, -1, 0, Inf, NaN,
  numeric(0))

for(x in bad_mult) {
  i <- mark()
  report(threw(normalize_quantile(mk_state(mfull), cfg0, multiplier=x)) &&
      logged("normalize_quantile: multiplier has to be a single finite", i),
    paste0("normalize_quantile() multiplier ", paste(x, collapse=","),
      " (class ", class(x), ", length ", length(x), ") is refused"))
  report(no_raw(i), "  and no raw R message leaks out")

  i <- mark()
  report(threw(normalize_cpm(mk_state(mfull), cfg0, multiplier=x)) &&
      logged("normalize_cpm: multiplier has to be a single finite", i),
    paste0("normalize_cpm() multiplier ", paste(x, collapse=","),
      " (class ", class(x), ", length ", length(x), ") is refused"))
  report(no_raw(i), "  and no raw R message leaks out")
}

## and a valid multiplier still does what it did, which is scale the result:

out1 <- val(normalize_quantile(mk_state(mfull), cfg0, multiplier=1e3))
out2 <- val(normalize_quantile(mk_state(mfull), cfg0, multiplier=2e3))
report(!is.null(out1) && !is.null(out2) &&
    isTRUE(all.equal(out2$expression, 2 * out1$expression)),
  "a valid multiplier still scales the result and nothing else")
report(!is.null(out1) && isTRUE(all.equal(out1$expression,
  normalize_quantile(mk_state(mfull), cfg0)$expression, tolerance=0)),
  "  and 1e3 is still the default")

## a method name, in the two functions that pick one from a list of their own:

bad_meth <- list("bogus", c("RLE", "TMM"), NA, NA_character_, 42, TRUE)

for(x in bad_meth) {
  i <- mark()
  report(threw(normalize_edger(mk_state(mfull), cfg0, method=x)) &&
      logged("normalize_edger: method has to be one of RLE upperquartile", i),
    paste0("normalize_edger() method of class ", class(x), ", length ",
      length(x), ", is refused"))
  report(no_raw(i), "  and no raw R message leaks out")

  i <- mark()
  report(threw(normalize_mscoreutils(mk_state(mfull), cfg0, method=x)) &&
      logged("normalize_mscoreutils: method has to be one of sum max", i),
    paste0("normalize_mscoreutils() method of class ", class(x), ", length ",
      length(x), ", is refused"))
  report(no_raw(i), "  and no raw R message leaks out")
}

## NULL and a zero-length value mean unset, which is what NULL already meant:
##   the default, or config, or an error naming both:

for(x in list(NULL, numeric(0))) {
  out <- val(normalize_quantile(mk_state(mfull), cfg0, normalization_quantile=x))
  report(!is.null(out) && isTRUE(all.equal(out$expression,
      normalize_quantile(mk_state(mfull), cfg0,
        normalization_quantile=0.75)$expression, tolerance=0)),
    paste0("normalization_quantile of length ", length(x),
      " means unset, so 0.75 is used"))
  out <- val(normalize_loess(mk_state(mfull), cfg0, span=x))
  report(!is.null(out) && isTRUE(all.equal(out$expression,
      limma::normalizeCyclicLoess(mfull, adaptive.span=TRUE, method="fast"),
      tolerance=0)),
    paste0("span of length ", length(x), " means unset, so the span is adaptive"))
}

for(x in list(NULL, character(0))) {
  out <- val(normalize_edger(mk_state(mfull), cfg0, method=x))
  report(!is.null(out) && isTRUE(all.equal(out$expression,
      normalize_edger(mk_state(mfull), cfg0, method="RLE")$expression,
      tolerance=0)),
    paste0("edger method of length ", length(x),
      " means unset, so RLE is used"))
  i <- mark()
  report(threw(normalize_mscoreutils(mk_state(mfull), cfg0, method=x)) &&
      logged("normalize_mscoreutils: method and config$normalization_method", i),
    paste0("mscoreutils method of length ", length(x),
      " with nothing in config is 'both unset'"))
  cfg_m <- cfg0
  cfg_m$normalization_method <- "sum"
  out <- val(normalize_mscoreutils(mk_state(mfull), cfg_m, method=x))
  report(!is.null(out), "  and comes from config when config names one")
}

## the config path of the two functions that read a value from config: these
##   two call check_config() as well now (see the section below), so a bad
##   config value is named there, as a config problem, before either guard sees
##   it. The guards above are what covers the same value passed as an argument:

cfg_bad <- cfg0
cfg_bad$normalization_span <- c(0.3, 0.7)
i <- mark()
report(threw(normalize_loess(mk_state(mfull), cfg_bad)) &&
    logged("check_config: param not scalar proportion", i),
  "a config$normalization_span of length 2 is refused, as a config problem")
cfg_bad$normalization_span <- NA_real_
i <- mark()
report(threw(normalize_loess(mk_state(mfull), cfg_bad)) &&
    logged("check_config: param is NA or NaN", i),
  "and a config$normalization_span of NA_real_ likewise")
cfg_bad <- cfg0
cfg_bad$normalization_method <- c("sum", "max")
i <- mark()
report(threw(normalize_mscoreutils(mk_state(mfull), cfg_bad)) &&
    logged("check_config: param not scalar character", i),
  "a config$normalization_method of length 2 is refused, as a config problem")

###############################################################################
section("every normalize_*() checks the config it was handed")

nrm <- list(
  normalize_edger=function(st, cfg) normalize_edger(st, cfg, method="RLE"),
  normalize_quantile=function(st, cfg) normalize_quantile(st, cfg),
  normalize_cpm=function(st, cfg) normalize_cpm(st, cfg),
  normalize_vsn=function(st, cfg) normalize_vsn(st, cfg),
  normalize_loess=function(st, cfg) normalize_loess(st, cfg),
  normalize_qquantile=function(st, cfg) normalize_qquantile(st, cfg),
  normalize_mscoreutils=function(st, cfg) normalize_mscoreutils(st, cfg,
    method="sum"),
  normalize=function(st, cfg) normalize(st, cfg, method="cpm")
)

report(length(nrm) == 8,
  "all eight entry points are covered: the five that did not check, and the three that did")

## a config carrying a name check_config() does not know is refused by every one
##   of them, not only by the three that always called it:

cfg_junk <- cfg0
cfg_junk$bogus_param <- 1

for(nom in names(nrm)) {
  i <- mark()
  report(threw(nrm[[nom]](mk_state(mfull), cfg_junk)) &&
      logged("check_config: unexpected parameter name: bogus_param", i),
    paste0(nom, "() refuses a config carrying an unrecognized parameter name"))
}

## and it is the first thing each of them does: handed a bad config and an
##   expression matrix that is not a matrix, every one of them reports the
##   config:

st_bad <- mk_state(as.data.frame(mfull))

for(nom in names(nrm)) {
  i <- mark()
  report(threw(nrm[[nom]](st_bad, cfg_junk)) &&
      logged("check_config: unexpected parameter name", i) &&
      !logged("is.matrix", i),
    paste0(nom, "() checks the config before it looks at the state"))
  report(threw(nrm[[nom]](st_bad, cfg0)) && logged("is.matrix", i),
    paste0("  while with a good config ", nom, "() gets as far as the state"))
}

## the value checks come with it, which is what the two functions that read a
##   value from config get out of this: a wrong class is named as a config
##   problem before the argument guard sees it at all:

cfg_s <- cfg0
cfg_s$normalization_span <- "0.7"
i <- mark()
report(threw(normalize_loess(mk_state(mfull), cfg_s)) &&
    logged("check_config: param not scalar proportion", i),
  "a config$normalization_span of the wrong class is caught by check_config()")
report(!logged("normalize_loess: span has to be", i),
  "  rather than by the span guard, which would have named the argument")

cfg_m <- cfg0
cfg_m$normalization_method <- 42
i <- mark()
report(threw(normalize_mscoreutils(mk_state(mfull), cfg_m)) &&
    logged("check_config: param not scalar character", i),
  "a config$normalization_method of the wrong class is caught by check_config()")

## and an empty config still works where the documentation says it can:

for(nom in c("normalize_cpm", "normalize_vsn", "normalize_qquantile")) {
  report(!is.null(val(nrm[[nom]](mk_state(mfull), list()))),
    paste0(nom, "() still takes an empty config, as its @param config says"))
}

###############################################################################
section("one feature: the methods that collapse are warned about")

## the list in f.one_feature_methods() is checked against what the methods
##   actually do rather than against itself: for each method, whether
##   normalize() warned, and whether the result came back with one distinct
##   value for the whole matrix:

warn_txt <- "every observation will come back with the same value"

mk1b <- function(e) {
  list(expression=e, features=data.frame(pep=rownames(e), gene=rownames(e)),
    samples=samps)
}

warned <- character(0)
collapsed <- character(0)
refused <- character(0)

for(meth in normalize_methods()) {
  i <- mark()
  q <- if(meth %in% "quantile") 0.75 else NULL
  out <- val(normalize(mk1b(mfull[1, , drop=FALSE]), cfg0, method=meth,
    normalization_quantile=q))
  if(logged(warn_txt, i)) warned <- c(warned, meth)
  if(is.null(out)) {
    refused <- c(refused, meth)
  } else if(length(unique(signif(as.vector(out$state$expression), 10))) == 1) {
    collapsed <- c(collapsed, meth)
  }
}

report(identical(sort(warned), sort(collapsed)),
  paste0("normalize() warns about exactly the methods that collapse: ",
    length(warned), " warned, ", length(collapsed), " collapsed"))
report(identical(sort(warned), sort(f.one_feature_methods())),
  "and f.one_feature_methods() is that same list, not a list of its own")
report(length(warned) == 10,
  paste0("ten of the sixteen collapse a single feature: ",
    paste(sort(warned), collapse=" ")))
report(identical(sort(refused), sort(c("qquantile", "vsn"))),
  "the two that refuse a single feature are still qquantile and vsn")
report(setequal(setdiff(normalize_methods(), c(warned, refused)),
    c("sum", "max", "log2", "none")),
  "sum, max, log2 and none keep a single feature's six distinct values")

## why those two are exempt, measured rather than asserted: MsCoreUtils
##   normalizes sum and max per feature, so nothing about observations was
##   equalized and one feature keeps its differences:

m60 <- mfull
storage.mode(m60) <- "double"
out_sum <- MsCoreUtils::normalize_matrix(m60, method="sum")
report(isTRUE(all.equal(out_sum, m60 / rowSums(m60))) &&
    !isTRUE(all.equal(out_sum, t(t(m60) / colSums(m60)))),
  "MsCoreUtils 'sum' is m / rowSums(m), per feature and not per observation")
report(isTRUE(all.equal(MsCoreUtils::normalize_matrix(m60, method="max"),
  m60 / apply(m60, 1, max))),
  "  and 'max' is m / apply(m, 1, max), likewise")
report(isTRUE(all.equal(MsCoreUtils::normalize_matrix(m60, method="div.mean"),
  t(t(m60) / colMeans(m60)))),
  "  while 'div.mean' is per observation, which is why one feature collapses")

## it is a warning, not a refusal, and it says what it is about:

i <- mark()
out <- val(normalize(mk1b(mfull[1, , drop=FALSE]), cfg0, method="quantile",
  normalization_quantile=0.75))
report(!is.null(out), "the report is a warning: normalize() still returns a state")
report(logged("normalize: WARNING:", i) &&
    logged("state$expression has 1 feature", i),
  "  and names the feature count")
report(logged("sum and max come through", i),
  "  and says which methods do not have the problem")
report(!is.null(out) &&
    length(unique(signif(as.vector(out$state$expression), 10))) == 1,
  "  and the thing it warns about is what happens")

## and nothing is said from two features up, where the methods have something
##   to work with:

for(nf in c(2, 3, 60)) {
  i <- mark()
  out <- val(normalize(mk1b(mfull[1:nf, , drop=FALSE]), cfg0, method="quantile",
    normalization_quantile=0.75))
  report(!is.null(out) && !logged(warn_txt, i),
    paste0("no warning at ", nf, " features"))
}

## quantiles.robust is degenerate at two features as well, which this warning
##   does not cover: measured here so that it is on the record rather than
##   discovered again. Every row comes back constant:

out <- val(normalize(mk1b(mfull[1:2, , drop=FALSE]), cfg0,
  method="quantiles.robust"))
report(!is.null(out) &&
    all(apply(out$state$expression, 1,
      function(v) length(unique(signif(v, 10))) == 1)),
  "quantiles.robust still collapses each row at two features, and is not warned about")

###############################################################################
section("quantiles.robust is refused on any missing value")

## what the guard is for, measured against MsCoreUtils rather than against the
##   guard: one gap in an otherwise complete matrix, at feat_7 of obs_3.

mg <- mfull
mg[7, 3] <- NA
raw <- MsCoreUtils::normalize_matrix(mg, method="quantiles.robust")

report(sum(is.na(raw)) > sum(is.na(mg)),
  paste0("MsCoreUtils returns more missing values than it was given: ",
    sum(is.na(mg)), " in, ", sum(is.na(raw)), " out"))
report(all(colSums(is.na(raw)) > 0),
  "  one in every observation, not only in the one the gap was in")
report(sum(is.na(raw) & !is.na(mg)) > 0,
  "  and in cells that were measured, so values are relocated, not only lost")
report(length(unique(which(is.na(raw), arr.ind=TRUE)[, "row"])) > 1,
  "  spread over more than one feature, whichever sorts to the vacated rank")

## so it is refused, whatever the pattern of the missingness:

for(nm in c("one_gap", "empty_obs", "scattered")) {
  mm <- switch(nm, one_gap=mg, empty_obs=me, scattered=mz)
  i <- mark()
  report(threw(normalize_mscoreutils(mk_state(mm), cfg0,
    method="quantiles.robust")),
    paste0("quantiles.robust is refused on the ", nm, " matrix"))
  report(logged("cannot be used on data with missing values", i),
    paste0("  ", nm, ": and says what is wrong"))
  report(logged(paste("missing values in state$expression:", sum(is.na(mm))), i),
    paste0("  ", nm, ": and counts them"))
}

## and through normalize(), which is how the workflow reaches it:

i <- mark()
report(threw(normalize(mk_state(m0), cfg0, method="quantiles.robust")),
  "refused through normalize() as well")
report(logged("assigns values by rank within each observation", i),
  "  the refusal explains the mechanism")
report(logged("in whatever feature sorts", i),
  "  including where the missing value ends up")
report(logged("sum, max, div.mean, div.median", i) &&
    logged("normalize_qquantile", i),
  "  and names what to use instead")
report(!logged("INTEGER() can only be applied", i) &&
    !logged("must be an array", i),
  "  and no third-party message is what stopped it")

out <- val(normalize(mk_state(mfull), cfg0, method="quantiles.robust"))
report(!is.null(out) && sum(is.na(out$state$expression)) == 0,
  "a matrix with nothing missing still normalizes, and comes back complete")

## the other four methods here are not refused, which is the reason the guard
##   names one method rather than the function: they return the missing values
##   they were given, in the cells they were given them in:

for(meth in c("sum", "max", "div.mean", "div.median")) {
  out <- val(normalize_mscoreutils(mk_state(m0), cfg0, method=meth))
  report(!is.null(out) &&
      identical(which(is.na(out$expression)), which(is.na(m0))),
    paste0(meth, ": not refused, and leaves the missing values where they were"))
}

## tune() normalizes before imputing and does not wrap that call in try(), so a
##   method that refuses missing data would stop a whole sweep:

tune_meths <- eval(formals(tune)$normalization_methods)
report(!("quantiles.robust" %in% tune_meths),
  "tune()'s default normalization_methods no longer offers quantiles.robust")
report(all(c("qquantile", "RLE") %in% tune_meths),
  "  while the tolerant methods it does offer are still swept")

###############################################################################
section("the documentation agrees with the code")

man_dir <- file.path(dirname(r_dir), "man")
readme <- file.path(dirname(dirname(r_dir)), "readme.md")

src <- unlist(lapply(list.files(r_dir, pattern="[.]R$", full.names=TRUE),
  readLines, warn=FALSE))
docs <- src
if(dir.exists(man_dir)) {
  docs <- c(docs, unlist(lapply(list.files(man_dir, pattern="[.]Rd$",
    full.names=TRUE), readLines, warn=FALSE)))
} else {
  cat("## note: man/ not found at", man_dir, "-- Rd files not read\n")
}
if(file.exists(readme)) {
  docs <- c(docs, readLines(readme, warn=FALSE))
} else {
  cat("## note: readme.md not found at", readme, "-- not read\n")
}

## regression guards: a documented name that does not exist sends the reader
##   looking for a function that was never there:

report(exists("normalize_methods") && length(normalize_methods()) == 16,
  "normalize_methods() is the name, and it returns all 16 methods")
report(!any(grepl("normalization_methods()", docs, fixed=TRUE)),
  "nothing names h0testr::normalization_methods(), which does not exist")
report(!any(grepl("normalizaiton", docs, fixed=TRUE)),
  "the 'normalizaiton' typo has not come back")
report(!any(grepl('c("vsn","cpm","quantile","qquantile"', docs, fixed=TRUE)),
  "no config row lists a subset of the methods in place of the function")
report(!any(grepl("unexpected config$normalization_method", src, fixed=TRUE)),
  "nothing reports config$normalization_method as the method it rejected")
report(!any(grepl("inter-observation normalization of one feature is degenerate",
  docs, fixed=TRUE)),
  "the paragraph that named three of the ten collapsing methods is gone")
report(any(grepl("ten of the sixteen methods return the same", docs,
  fixed=TRUE)),
  "and the one that replaced it says how many there are")
report(any(grepl("assigns values by rank within each observation", docs,
  fixed=TRUE)),
  "the docs record why quantiles.robust cannot be used on data with gaps")

## normalize_loess() reads config$normalization_span, so its config row cannot
##   say it uses no keys:

i <- grep("^normalize_loess <- function", src)
j <- max(grep("^#' @param config", src[1:i]))
report(length(i) == 1 && !grepl("Does not use any keys", src[j]) &&
    grepl("normalization_span", paste(src[j:(j + 3)], collapse=" ")),
  "normalize_loess()'s config row names the key it reads")

cfg_s <- cfg0
cfg_s$normalization_span <- 0.2
report(!isTRUE(all.equal(
  normalize_loess(mk_state(mfull), cfg_s)$expression,
  normalize_loess(mk_state(mfull), cfg0)$expression)),
  "  and that key really is read")
report(isTRUE(all.equal(
  normalize_loess(mk_state(mfull), cfg0, span=0.2)$expression,
  normalize_loess(mk_state(mfull), cfg_s)$expression)),
  "  with the span argument overriding it")

if(dir.exists(man_dir)) {
  rd <- readLines(file.path(man_dir, "normalize_loess.Rd"), warn=FALSE)
  k <- grep("item{config}", rd, fixed=TRUE)
  report(length(k) == 1 && !grepl("Does not use any keys", rd[k]) &&
      grepl("normalization_span", paste(rd[k:(k + 3)], collapse=" ")),
    "  and normalize_loess.Rd says the same")

  rd <- readLines(file.path(man_dir, "normalize_mscoreutils.Rd"), warn=FALSE)
  report(any(grepl("any missing value", rd, fixed=TRUE)),
    "normalize_mscoreutils.Rd says the method is refused on data with gaps")

  rd <- readLines(file.path(man_dir, "tune.Rd"), warn=FALSE)
  report(!any(grepl("quantiles.robust", rd, fixed=TRUE)),
    "and tune.Rd no longer shows it in the default sweep")
}

###############################################################################

cat("\n## passes:", n_pass, "; failures:", n_fail, "; elapsed:",
  round(as.numeric(difftime(Sys.time(), t0, units="secs"))), "s\n")
cat("## log file:", log_file, "\n")

quit(save="no", status=if(n_fail > 0) 1 else 0)
