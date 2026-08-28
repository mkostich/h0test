## Tests for the four faults found in simulate.R on Aug 19 and fixed the same day.
##   f.sim_rnorm_pos() looped forever whenever a non-positive draw landed anywhere but the
##   first position and m and s were scalars, which is how sim1() and f.sim_tech_reps() call
##   it: the retry indexed a length 1 vector with a length n logical, redrew with mean NA,
##   got NaN, and never satisfied its own exit condition. f.pep_drop() lacked drop=F, so a
##   single surviving peptide became a vector and lost the rownames the peptide level
##   simulator reads gene labels from, and every peptide could be dropped. mnar_off was
##   accepted and documented by both simulators but never passed through to f.mnar(). And a
##   one row matrix came back from f.mnar()'s apply() as a vector, so sim1(n_feats=1) died on
##   nrow(mat) being NULL.
## Later the same day: sim1() validated none of its arguments and the peptide level simulator
##   only two of its, so both now check every argument through f.chk_num(); sim1()'s final
##   round() became ceiling(), round() being able to send a strictly positive draw to zero;
##   and sim1()'s inlined MCAR step became a call to f.mcar(), which already did the same
##   three lines.
## Aug 19, later still: sim2() was retired in favor of sim_design(), which carries the same
##   intensity, noise and missingness layers, so the assertions below that reached those
##   layers through sim2() now reach them through sim_design(). The ones that were about
##   sim2() itself, its return shape and its own refusals, went with it; the four argument
##   checks it was the only caller of moved to test_sim_design.R.

usage <- function(msg=NULL) {
  if(!is.null(msg)) cat("ERROR:", msg, "\n\n", file=stderr())
  cat(
    "Test that h0testr's simulators terminate, keep their matrices matrices, and honor",
    "every argument they document: that f.sim_rnorm_pos() redraws correctly with scalar",
    "parameters however heavy the rejection rate, that f.pep_drop() returns a labelled",
    "matrix in which every gene keeps at least one of its peptides, that mnar_off",
    "reaches f.mnar() from sim1() and sim_design(),",
    "and that a single feature or a single observation is a workable simulation. Also",
    "that every documented argument constraint is enforced with a message naming the",
    "argument and the offending value, and that the returned matrix is whole numbers of",
    "at least 1.",
    "",
    "Usage: Rscript test_simulate.R <r_dir>",
    "",
    "Required positional arguments:",
    "  <r_dir>  Path to the h0testr package R/ source directory; all .R files",
    "             there are sourced (the installed package is not used).",
    "",
    "Optional named arguments: none.",
    "",
    "Needs no packages beyond stats; the simulators draw from stats::rnorm() and",
    "  stats::rbinom() only. The rejection sampling assertions are wrapped in an elapsed",
    "  time limit, so a return of the infinite loop fails rather than hanging the suite.",
    "",
    "Output: one PASS/FAIL line per assertion to stdout, then a count of passes and",
    "  failures and the elapsed time.",
    "",
    "Exit codes: 0 all assertions passed; 1 one or more failed; 2 usage error.",
    "",
    "Examples:",
    "  Rscript test_simulate.R C:/path/to/h0test/h0testr/R",
    "  Rscript test_simulate.R ../../h0test/h0testr/R",
    "  Rscript test_simulate.R C:/path/to/h0testr/R > simulate.out 2>&1",
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

## the fault under test was non-termination, so every call that used to hang is run under an
##   elapsed limit; setTimeLimit() interrupts an R level loop, which turns a regression into a
##   failing assertion instead of a suite that never finishes:

bounded <- function(expr, secs=20) {
  setTimeLimit(elapsed=secs, transient=TRUE)
  on.exit(setTimeLimit(elapsed=Inf, transient=FALSE))
  return(try(expr, silent=TRUE))
}

###############################################################################
section("f.sim_rnorm_pos() with scalar parameters, where the loop used to be")

set.seed(3)
rslt <- bounded(sim1(n_obs=200, n_feats=200, log_m_mean=0.2, log_m_sd=0.5,
  mnar_c0=-Inf, mnar_c1=0, mcar_p=0))
report(!inherits(rslt, "try-error") && is.matrix(rslt$mat),
  "a mean well inside the rejection region returns instead of spinning")
report(!inherits(rslt, "try-error") && all(rslt$feat_mean > 0) && !any(is.nan(rslt$feat_mean)),
  "and every feature mean is a real positive number, not a NaN from a mean of NA")

set.seed(11)
v <- bounded(f.sim_rnorm_pos(n=2000, m=0.1, s=1))     ## roughly 46% of draws need a retry
report(!inherits(v, "try-error") && length(v) == 2000 && all(v > 0) && !any(is.nan(v)),
  "scalar m and s recycle, so a retry at any position redraws with that position's mean")

set.seed(12)
v <- bounded(f.sim_rnorm_pos(n=50, m=rep(0.1, 50), s=rep(1, 50)))
report(!inherits(v, "try-error") && all(v > 0),
  "vector parameters, the case that always worked, still work")

set.seed(13)
v <- bounded(f.sim_rnorm_pos(n=6, m=c(10, 1000), s=c(1, 10)))
report(!inherits(v, "try-error") && length(v) == 6 && all(v > 0),
  "and parameters shorter than n recycle rather than running off the end as NA")

why <- function(expr) {
  rslt <- try(expr, silent=TRUE)
  if(!inherits(rslt, "try-error")) return("")
  return(conditionMessage(attr(rslt, "condition")))
}

report(grepl("no positive value is possible", why(f.sim_rnorm_pos(3, 0, 0))),
  "m and s both non-positive is refused: no draw could ever satisfy the loop")
report(grepl("must be finite", why(f.sim_rnorm_pos(3, NA, 1))),
  "a non-finite mean is named as such, not left to if(any(m < 0)) as an NA condition")
report(grepl("must be finite", why(f.sim_rnorm_pos(3, 10, Inf))), "and a non-finite sd")
report(grepl("m < 0", why(f.sim_rnorm_pos(c(3), c(1, -1, 1), c(1, 1, 1)))),
  "a negative mean is still refused, as before, and its position reported")
report(grepl("2:", why(f.sim_rnorm_pos(c(3), c(1, -1, 1), c(1, 1, 1)))),
  "by position rather than by dumping every parameter into the message")
report(threw(f.sim_rnorm_pos(3, c(10, 20), 1)),
  "and mismatched parameter lengths, also as before")

###############################################################################
section("f.pep_drop(), which has to hand back a matrix with rownames")

mat <- matrix(1:12, nrow=4, dimnames=list(paste0("gene1_pep", 1:4), paste0("obs", 1:3)))

set.seed(7)
out <- f.pep_drop(mat, peps_per_gene=2, p_drop=0.99)
report(is.matrix(out) && !is.null(rownames(out)) && nrow(out) >= 1,
  "a lone survivor stays a one row matrix and keeps its label")
report(all(rownames(out) %in% rownames(mat)) && ncol(out) == 3,
  "the surviving rows are the original rows, columns untouched")

out <- f.pep_drop(mat, peps_per_gene=2, p_drop=1)
report(is.matrix(out) && nrow(out) == 1,
  "one gene dropped entirely keeps one peptide rather than emptying the matrix or refusing")

out <- f.pep_drop(mat, peps_per_gene=1, p_drop=0.99)
report(identical(out, mat), "one peptide per gene means nothing to drop, whatever p_drop says")

out <- f.pep_drop(mat, peps_per_gene=4, p_drop=0)
report(identical(out, mat), "and p_drop of zero drops nothing")

## the floor is per gene, not per matrix. It used to be the latter, so a gene could lose every
##   peptide and drop out of the simulation and out of the truth table it is scored against,
##   and only a matrix with nothing left at all was patched: p_drop of 1 over several genes
##   returned a single feature of a single gene:

gene4 <- rep(paste0("gene", 1:4), each=3)
mat4 <- matrix(1, nrow=12, ncol=3, dimnames=list(paste0(gene4, "_pep", 1:3),
  paste0("obs", 1:3)))

out <- f.pep_drop(mat4, peps_per_gene=3, p_drop=1)
report(nrow(out) %in% 4 &&
  identical(sort(unique(sub("_.*", "", rownames(out)))), sort(unique(gene4))),
  "p_drop of 1 leaves one peptide of every gene, not one peptide overall")

n_bad <- 0
for(k in 1:40) {
  set.seed(200 + k)
  cnt <- table(factor(sub("_.*", "", rownames(f.pep_drop(mat4, peps_per_gene=3,
    p_drop=0.75))), levels=unique(gene4)))
  if(any(cnt < 1) || any(cnt > 3)) n_bad <- n_bad + 1
}
report(n_bad %in% 0,
  "and every gene keeps between one and peps_per_gene peptides over 40 seeds")

## the gene of each row may be named instead of read off the labels, which is how
##   sim_design() calls it: a feature id is not the gene id with a suffix in general. Without
##   the argument these six rows are six genes and nothing can be dropped at all:

mat6 <- matrix(1, nrow=6, ncol=2, dimnames=list(paste0("f", 1:6), c("o1", "o2")))

report(nrow(f.pep_drop(mat6, peps_per_gene=3, p_drop=1)) %in% 6 &&
  nrow(f.pep_drop(mat6, peps_per_gene=3, p_drop=1,
    genes=rep(c("gA", "gB"), each=3))) %in% 2,
  "a genes argument names the gene of each row where the labels do not")
report(threw(f.pep_drop(mat6, peps_per_gene=3, p_drop=1, genes=c("gA", "gB"))),
  "and one that does not line up with the rows is refused")

## the fallback strips the _pep<n> suffix rather than cutting at the first underscore, which made
##   a single gene of every gene id containing one: these six rows are two genes, not one, so
##   p_drop of 1 leaves two features rather than one:

mat_u <- matrix(1, nrow=6, ncol=2, dimnames=list(paste0(rep(c("gene_a", "gene_b"), each=3),
  "_pep", 1:3), c("o1", "o2")))
set.seed(11)
out <- f.pep_drop(mat_u, peps_per_gene=3, p_drop=1)
report(nrow(out) %in% 2 &&
  identical(sort(sub("_pep[0-9]+$", "", rownames(out))), c("gene_a", "gene_b")),
  "the label fallback splits on the peptide suffix, so an underscore in a gene id is safe")

###############################################################################
section("f.sim_tech_reps() gets no NAs, and no longer pretends otherwise")

## it ran before f.mnar() and f.mcar() and its result is checked, so the is.na(val) branch that
##   widened an NA into a column of NAs was unreachable. Removed: an NA arriving here is named
##   and refused by f.sim_rnorm_pos() rather than quietly spread across the replicates:

mat_r <- matrix(c(100, NA), nrow=2, dimnames=list(c("f1", "f2"), "o1"))
report(grepl("must be finite", why(f.sim_tech_reps(mat_r, reps_per_sample=2, cv_reps=0.1))),
  "an NA value is refused rather than widened, nothing in the package sending one")
report(identical(f.sim_tech_reps(mat_r, reps_per_sample=1, cv_reps=0.1), mat_r),
  "and one replicate per sample is still a passthrough, NA or not")

mat_r <- matrix(c(100, 200), nrow=2, dimnames=list(c("f1", "f2"), "o1"))
set.seed(12)
out <- f.sim_tech_reps(mat_r, reps_per_sample=3, cv_reps=0.1)
report(is.matrix(out) && nrow(out) == 2 && ncol(out) == 3 && all(out > 0) &&
  identical(colnames(out), paste0("o1_rep", 1:3)),
  "and the replicates it does make are positive and labelled by sample and replicate")

n_bad <- 0
for(k in 1:40) {
  set.seed(100 + k)
  samps <- sim_samples(factors=list(grp=c("ctl", "trt")), n_per_cell=2)
  rslt <- try(sim_design(samps, frm=~grp, test_term="grp", n_genes=2, n_genes_signif=1,
    peps_per_gene=2), silent=TRUE)
  mat <- if(inherits(rslt, "try-error")) NULL else rslt$state$expression
  if(is.null(mat) || !is.matrix(mat) || is.null(rownames(mat)) || nrow(mat) < 1) {
    n_bad <- n_bad + 1
  }
}
report(n_bad == 0,
  "the smallest peptide level simulation, four rows against the default p_drop, survives 40 seeds")

###############################################################################
section("mnar_off, accepted and documented by both simulators")

## the offset only matters where it is an appreciable fraction of the intensity, so these use
##   feature means near 1; the count of missing values then has to move when the offset does.
##   Both offsets are inside the documented open interval, which the checks now enforce:

set.seed(21)
a <- sim1(n_obs=40, n_feats=300, log_m_mean=0, log_m_sd=0.2, mnar_c0=0, mnar_c1=-3,
  mnar_off=1e-6, mcar_p=0)
set.seed(21)
b <- sim1(n_obs=40, n_feats=300, log_m_mean=0, log_m_sd=0.2, mnar_c0=0, mnar_c1=-3,
  mnar_off=0.99, mcar_p=0)
report(sum(is.na(a$mat)) != sum(is.na(b$mat)),
  "sim1 passes it through to f.mnar() rather than leaving the helper default in place")

set.seed(22)
samps <- sim_samples(factors=list(grp=c("ctl", "trt")), n_per_cell=6)
a <- sim_design(samps, frm=~grp, test_term="grp", n_genes=200, log_m_mean=0, log_m_sd=0.2,
  mnar_c0=0, mnar_c1=-3, mnar_off=1e-6, mcar_p=0, p_drop=0)
set.seed(22)
samps <- sim_samples(factors=list(grp=c("ctl", "trt")), n_per_cell=6)
b <- sim_design(samps, frm=~grp, test_term="grp", n_genes=200, log_m_mean=0, log_m_sd=0.2,
  mnar_c0=0, mnar_c1=-3, mnar_off=0.99, mcar_p=0, p_drop=0)
report(sum(is.na(a$state$expression)) != sum(is.na(b$state$expression)),
  "and so does sim_design")

###############################################################################
section("a single feature or a single observation")

set.seed(31)
rslt <- try(sim1(n_obs=4, n_feats=1), silent=TRUE)
report(!inherits(rslt, "try-error") && identical(dim(rslt$mat), c(1L, 4L)),
  "sim1(n_feats=1) returns a one row matrix, as its own @param allows")
report(!inherits(rslt, "try-error") && identical(rownames(rslt$mat), "feat_1") &&
  identical(colnames(rslt$mat), paste0("obs_", 1:4)), "with both sets of labels intact")

set.seed(32)
rslt <- try(sim1(n_obs=1, n_feats=1), silent=TRUE)
report(!inherits(rslt, "try-error") && identical(dim(rslt$mat), c(1L, 1L)),
  "and the 1x1 corner of that is a matrix too")

m1 <- matrix(c(100, 200, 300), nrow=1, dimnames=list("f1", c("o1", "o2", "o3")))
out <- f.mnar(m1, mnar_c0=-Inf, mnar_c1=0)
report(is.matrix(out) && identical(dimnames(out), dimnames(m1)),
  "f.mnar() restores the shape apply() drops on a one row matrix")
report(all(out == m1), "leaving the values alone where p(mnar) is zero")

out <- f.mnar(matrix(c(100, 200), nrow=2, dimnames=list(c("f1", "f2"), "o1")),
  mnar_c0=-Inf, mnar_c1=0)
report(is.matrix(out) && identical(dim(out), c(2L, 1L)), "and a one column matrix as well")

###############################################################################
section("f.mnar() at the ends of the logistic, where exp(resp) / (1 + exp(resp)) overflowed")

warned <- function(expr) {               ## the value and the warning, for asserting on both
  msg <- ""
  val <- withCallingHandlers(expr, warning=function(w) {
    msg <<- conditionMessage(w)
    invokeRestart("muffleWarning")
  })
  return(list(val=val, msg=msg))
}

## the ratio is NaN once mnar_c0 + mnar_c1 * log(intensity) is above about 710, and
##   stats::rbinom(prob=NaN) is NA, whose all NA logical index assigns nothing: a cell certain
##   to be missing was left in place, with only an "NAs produced" warning to show for it.
##   stats::plogis() saturates instead. Both ends are asserted, the low one being how the
##   documented mnar_c0=-Inf turns MNAR off:

big <- matrix(1e5, nrow=4, ncol=3, dimnames=list(paste0("f", 1:4), paste0("o", 1:3)))

set.seed(33)
out <- warned(f.mnar(big, mnar_c0=4.65, mnar_c1=62))
report(all(is.na(out$val)) && !nzchar(out$msg),
  "a response above the overflow of the ratio drops every cell, and does so quietly")

set.seed(34)
out <- warned(f.mnar(big, mnar_c0=-1e4, mnar_c1=0))
report(!any(is.na(out$val)) && !nzchar(out$msg),
  "and a response far below it drops none, p(mnar) underflowing to 0 rather than to NaN")

set.seed(35)
out <- warned(sim1(n_obs=3, n_feats=6, mnar_c1=62, mcar_p=0))
report(all(is.na(out$val$mat)) && !nzchar(out$msg),
  "the same through sim1, which used to hand back a matrix with no missing values at all")

set.seed(36)
rslt <- sim1(n_obs=4, n_feats=20, mnar_c0=-Inf, mnar_c1=0, mcar_p=0)
report(!any(is.na(rslt$mat)), "while mnar_c0 of -Inf still turns MNAR off, as it did before")

###############################################################################
section("argument checking, which sim1() had none of")

## every argument of sim1() is checked by f.chk_num(), so each is probed once out of range
##   and, where the boundary is legal, once on it. The message has to name the argument and
##   the offending value, since a bad simulation parameter is otherwise noticed only much
##   later, as a strange result. sim_design()'s own checks are in test_sim_design.R:

report(grepl("sim1: n_obs", why(sim1(0, 5))) && grepl("got: 0", why(sim1(0, 5))),
  "sim1 names the argument and the value it refused")
report(grepl("whole number", why(sim1(2.5, 5))), "n_obs must be a whole number")
report(grepl("scalar", why(sim1(c(2, 3), 5))), "and a scalar")
report(threw(sim1(NA, 5)) && threw(sim1("4", 5)), "not NA and not a character string")
report(threw(sim1(4, 0)), "n_feats below 1 is refused")
report(!threw(sim1(1, 1)), "while the smallest legal simulation is allowed")

report(threw(sim1(4, 5, log_m_mean=-1)), "log_m_mean below 0 is refused")
report(threw(sim1(4, 5, log_m_sd=-1)) && threw(sim1(4, 5, log_cv_sd=-1)),
  "as is a negative sd for either parameter")
report(!threw(sim1(4, 5, log_m_mean=0, log_cv_sd=0)) && !threw(sim1(4, 5, log_m_sd=0)),
  "zero is on the legal side of each of the three")
report(grepl("log_m_mean and log_m_sd cannot both be 0", why(sim1(4, 5, log_m_mean=0,
  log_m_sd=0))) && grepl("log_m_sd: 0", why(sim1(4, 5, log_m_mean=0, log_m_sd=0))),
  "though zero for both is refused by sim1 itself, naming both arguments and their values")
report(!threw(sim1(4, 5, log_cv_mean=-2)),
  "log_cv_mean is unconstrained: the default of -0.75 is itself negative")

report(grepl("must be finite", why(sim1(4, 5, mnar_c0=Inf))),
  "mnar_c0 of +Inf is refused, every cell being missing there")
report(!threw(sim1(4, 5, mnar_c0=-Inf, mnar_c1=0)),
  "but -Inf is allowed, that being the documented way to turn MNAR off")
report(threw(sim1(4, 5, mnar_c1=NA)), "mnar_c1 must be present and finite")
report(threw(sim1(4, 5, mnar_off=0)) && threw(sim1(4, 5, mnar_off=1)),
  "mnar_off is refused at both ends of its documented open interval")
report(threw(sim1(4, 5, mcar_p=-0.1)) && threw(sim1(4, 5, mcar_p=1.2)),
  "and mcar_p outside 0 to 1")

rslt <- sim1(4, 5, mnar_c0=-Inf, mnar_c1=0, mcar_p=1)
report(all(is.na(rslt$mat)),
  "mcar_p of 1 blanks the matrix, so sim1 still drops cells now that f.mcar() does it")

###############################################################################
section("ceiling() in sim1(), where round() could undo the positivity guarantee")

## with a feature mean near 1, round() sent draws below 0.5 to zero, which is exactly what
##   f.sim_rnorm_pos() exists to prevent, and what the peptide level simulator already
##   avoided with ceiling():

set.seed(51)
rslt <- sim1(n_obs=30, n_feats=200, log_m_mean=0, log_m_sd=0.1, mnar_c0=-Inf, mnar_c1=0,
  mcar_p=0)
report(all(rslt$mat >= 1), "no cell rounds down to zero when the means sit near 1")
report(all(rslt$mat %% 1 == 0), "and the matrix is still whole numbers")

set.seed(52)
rslt <- sim1(n_obs=6, n_feats=40)
report(all(rslt$mat %% 1 == 0, na.rm=TRUE) && all(rslt$mat >= 1, na.rm=TRUE),
  "at the defaults too, where rounding never mattered")

###############################################################################
section("what the fixes were not meant to change")

## the peptide level assertions that were here, on the five elements sim2() returned and on
##   its own refusals, went with it; the corresponding claims about sim_design()'s return and
##   its refusals are in test_sim_design.R, as is technical replication:

set.seed(41)
rslt <- sim1(n_obs=6, n_feats=8)
report(identical(names(rslt), c("mat", "feat_mean", "feat_cv")),
  "sim1 returns the same three elements")
report(all(rownames(rslt$mat) == names(rslt$feat_mean)) &&
  all(rownames(rslt$mat) == names(rslt$feat_cv)),
  "still aligned row for row with the parameters that generated them")
report(all(rslt$mat >= 1, na.rm=TRUE), "and still strictly positive where not missing")

set.seed(43)
rslt <- sim1(n_obs=6, n_feats=40, mnar_c0=-Inf, mnar_c1=0, mcar_p=0)
report(!any(is.na(rslt$mat)), "the documented no-missing-values settings still give none")

###############################################################################
cat("\n## passes:", n_pass, "; failures:", n_fail, "; elapsed:",
  round(as.numeric(difftime(Sys.time(), t0, units="secs"))), "s\n")

quit(save="no", status=if(n_fail > 0) 1 else 0)
