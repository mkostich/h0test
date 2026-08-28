## Tests five things about h0testr::tune_check(): that it joins the unpermuted results to
##   the permuted ones on named columns rather than on the first nine, that it says so when
##   the two sides do not cover the same set of parameter combinations, that it does not
##   count a permuted run which tested nothing as a permutation that found no false
##   positives, that it reads only the files belonging to the sweep it was given the prefix
##   of, and that it refuses a directory with no permuted file at the point the count of
##   them is taken rather than several steps downstream.
##   The key was apply(dat0[, 1:9], 1, paste, collapse=":") on both sides, the first nine
##   columns being the parameter combination as f.tune2() happens to emit it. Add a column
##   to that data.frame, or reorder it, and the key silently changes meaning: it still
##   pastes nine values together, and if both sides changed together it still matches, just
##   on the wrong nine. Every fdr in the table is then computed against the wrong permuted
##   rows and nothing errors. Naming the nine makes a missing name a stop and makes an
##   inserted or reordered column a no-op, and both are checked here.
##   The join itself was also silent about combinations it could not match. perm_max[k0] is
##   NA for an unpermuted combination with no permuted counterpart, so that combination got
##   an NA fdr, which the sort files with the combinations that never ran: it did run, and
##   there was simply nothing to measure it against. A permuted combination missing from
##   the unpermuted table went the other way and was dropped without trace. Both directions
##   are now counted and warned about, and both counts are checked here.
##   The fixture is the package's own inst/extdata/tune: one unpermuted file and five
##   permuted ones, 170 parameter combinations each, the six of them in exact agreement, so
##   every variant below is a named departure from a table known to line up.
##   The fixture's own rows all have ntests above 0, so the permuted rows that tested
##   nothing are made here rather than found, and the 30 msqrob rows it does have with
##   nhits NA and ntests 30 are left alone: those are tests that ran.
##   Elsewhere in this suite, test_marginality.R covers what tune_check() does with a
##   combination that never ran, and test_tune.R covers the sweep that writes these files.

usage <- function(msg=NULL) {
  if(!is.null(msg)) cat("ERROR:", msg, "\n\n", file=stderr())
  cat(
    "Test how h0testr::tune_check() joins unpermuted tuning results to permuted ones:",
    "that the join key is built from named columns, so that inserting or reordering a",
    "column in either file changes no result and dropping or renaming one of the nine",
    "key columns is refused with a message naming it, and that a combination present in",
    "one of the two sets and absent from the other is counted and warned about in both",
    "directions rather than silently given no fdr or silently dropped. Also tests that a",
    "permuted result row with ntests 0 is dropped from the permutation summaries rather",
    "than counted as a permutation that found no false positives, so that a combination",
    "skipped in every permutation gets no fdr instead of an fdr of 0, and that only files",
    "beginning with the given prefix are read, so that two sweeps sharing a directory and",
    "a suffix do not pool their results. Also tests that a directory holding no permuted",
    "file is refused with a message naming the missing files, rather than reaching an",
    "empty table several steps later and failing with an R level message about it.",
    "",
    "Usage: Rscript test_tune_check.R <r_dir>",
    "",
    "Required positional arguments:",
    "  <r_dir>  Path to the h0testr package R/ source directory; all .R files there",
    "             are sourced (the installed package is not used). The tuning results",
    "             in ../inst/extdata/tune alongside it are the fixture.",
    "",
    "Output: one 'PASS: <what>' or 'FAIL: <what>' line per assertion on stdout,",
    "  section headers with elapsed seconds, and a final count of passes and",
    "  failures. Exit code 0 if every assertion passed, 1 if any failed, 2 on a",
    "  usage error.",
    "",
    "Examples:",
    "  Rscript test_tune_check.R C:/path/to/h0testr/R",
    "  Rscript test_tune_check.R ../../h0test/h0testr/R",
    "  Rscript test_tune_check.R C:/Users/me/projects/h0testr/R",
    "",
    sep="\n", file=stderr()
  )
  quit(save="no", status=2)
}

args <- commandArgs(trailingOnly=TRUE)
if(length(args) != 1) usage("wrong number of arguments")

r_dir <- args[1]
if(!dir.exists(r_dir)) usage(paste("r_dir not a directory:", r_dir))

for(f in list.files(r_dir, pattern="[.]R$", full.names=TRUE)) source(f)

###############################################################################
## harness:

log_file <- tempfile(fileext=".log")
invisible(file.create(log_file))   ## mark() reads it before anything here has logged
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

mark <- function() length(readLines(log_file))
logged <- function(pat, since=0) {
  txt <- readLines(log_file)
  if(since >= length(txt)) return(FALSE)
  any(grepl(pat, txt[(since + 1):length(txt)], fixed=TRUE))
}

###############################################################################
## fixture:

## reached from the source tree the way test_check_config.R reaches the demo data, so
##   that a stale installed copy of the package cannot be what gets tested:

tune_src <- file.path(dirname(r_dir), "inst", "extdata", "tune")
if(!dir.exists(tune_src)) tune_src <- system.file("extdata/tune", package="h0testr")
if(!dir.exists(tune_src)) usage(paste("no tuning fixture found at:", tune_src))

suffix <- ".condition.tune.tsv"
cfg <- list(log_file=log_file)

key_cols <- c("norm", "nquant", "impute", "iquant", "scale", "span", "npcs", "k", "test")
out_cols <- c("nhits", "ntests", "fdr", "max1", "mid1", "avg1", "sd1", "norm", "nquant",
  "impute", "iquant", "scale", "span", "npcs", "k", "test")

## the fixture read once; every variant below is written out from these:

fix <- list()
for(f in list.files(tune_src, pattern=paste0(suffix, "$"))) {
  fix[[f]] <- utils::read.table(file.path(tune_src, f), header=T, sep="\t", quote="",
    as.is=T)
}

unperm_file <- paste0("0", suffix)
perm_files <- setdiff(names(fix), unperm_file)

## write a variant of the fixture into a fresh directory and return its path. fn is
##   applied to each data.frame, and is given the filename so that a variant can change
##   one side of the join and not the other:

mk_dir <- function(tag, fn=function(d, f) d) {
  dir_i <- file.path(tempdir(), paste0("tune_check_", tag))
  unlink(dir_i, recursive=TRUE)
  dir.create(dir_i, showWarnings=FALSE, recursive=TRUE)
  for(f in names(fix)) {
    utils::write.table(fn(fix[[f]], f), file=file.path(dir_i, f), sep="\t", quote=FALSE,
      row.names=FALSE)
  }
  return(dir_i)
}

chk <- function(dir_i) {
  try(suppressMessages(tune_check(dir_i, "", suffix, cfg)), silent=TRUE)
}

###############################################################################
section("baseline: the fixture as shipped")

dir_base <- mk_dir("base")
base <- chk(dir_base)

report(!inherits(base, "try-error"), "tune_check() runs on the shipped tuning fixture")

n_key <- length(unique(apply(fix[[unperm_file]][, key_cols], 1, paste, collapse=":")))

if(!inherits(base, "try-error")) {
  report(nrow(base) %in% nrow(fix[[unperm_file]]),
    "one output row per unpermuted input row")
  report(n_key %in% nrow(fix[[unperm_file]]),
    "the nine key columns identify each input row uniquely, so the join is one to one")
  report(identical(names(base), out_cols), "the output columns are the documented ones")
  report(!any(is.na(base$fdr)),
    "every combination gets an fdr, the two sides of the fixture covering the same grid")
  report(!any(is.na(base$max1)), "and every combination is matched to permuted results")
}

###############################################################################
section("item 4: the join key is built from named columns")

## a column inserted ahead of the nine, in both files. Under the positional key this
##   shifts the window by one, so "test" drops out of it and the five test methods of a
##   combination collide into one group: max1 becomes the largest permuted hit count of
##   any of the five, and every fdr that touches is wrong, with nothing said. Under the
##   named key the extra column is not part of the join and changes nothing:

dir_ins <- mk_dir("insert", function(d, f) cbind(run_id="r1", d, stringsAsFactors=FALSE))
res_ins <- chk(dir_ins)

report(!inherits(res_ins, "try-error"),
  "a column inserted ahead of the key columns does not stop the join")
report(!inherits(res_ins, "try-error") && !inherits(base, "try-error") &&
  identical(res_ins, base), "and changes no value in the result")

## the same, but reordering rather than inserting, and on the unpermuted side only.
##   Under the positional key the two sides then paste the same nine values in different
##   orders, so nothing matches at all and every max1 is NA:

shuffle <- function(d, f) {
  if(!(f %in% unperm_file)) return(d)
  return(d[, c(rev(key_cols), setdiff(names(d), key_cols)), drop=FALSE])
}

dir_ord <- mk_dir("reorder", shuffle)
res_ord <- chk(dir_ord)

report(!inherits(res_ord, "try-error"),
  "reordering the key columns of one file does not stop the join")
report(!inherits(res_ord, "try-error") && !inherits(base, "try-error") &&
  identical(res_ord, base), "and changes no value in the result")

## a key column that is not there at all. Renamed rather than dropped, so that the file
##   is otherwise the shape tune_check() expects and the refusal is about the name:

m <- mark()
dir_bad0 <- mk_dir("bad_unperm", function(d, f) {
  if(f %in% unperm_file) names(d)[names(d) %in% "norm"] <- "normalisation"
  return(d)
})
res_bad0 <- chk(dir_bad0)

report(inherits(res_bad0, "try-error"),
  "a key column missing from the unpermuted file is refused")
report(logged("is missing join column(s): norm", since=m),
  "and the message names the column that is missing")
report(logged(unperm_file, since=m), "and names the file it is missing from")
report(logged("normalisation", since=m), "and lists the columns that are there instead")

m <- mark()
dir_bad1 <- mk_dir("bad_perm", function(d, f) {
  if(!(f %in% unperm_file)) names(d)[names(d) %in% "test"] <- "test_method"
  return(d)
})
res_bad1 <- chk(dir_bad1)

report(inherits(res_bad1, "try-error"),
  "a key column missing from the permuted files is refused")
report(logged("permuted results are missing join column(s): test", since=m),
  "and the message says it was the permuted side and names the column")

## the two sides are checked separately, so a name present on one side and not the other
##   is caught rather than joined on whatever the two happen to share:

m <- mark()
dir_bad2 <- mk_dir("bad_both", function(d, f) {
  names(d)[names(d) %in% "iquant"] <- if(f %in% unperm_file) "impute_quantile" else "iq"
  return(d)
})
report(inherits(chk(dir_bad2), "try-error"),
  "a key column renamed differently on the two sides is refused")
report(logged("iquant", since=m), "naming iquant, the column neither side has")

###############################################################################
section("item 5: grids that do not line up are counted and warned about")

## drop one test method from the permuted files only. Every unpermuted combination using
##   it then has nothing to compare against, gets an NA fdr, and sorts in with the
##   combinations that never ran:

drop_perm <- function(d, f) {
  if(f %in% unperm_file) return(d)
  return(d[!(d$test %in% "proda"), , drop=FALSE])
}

n_proda <- sum(fix[[unperm_file]]$test %in% "proda")

m <- mark()
dir_m0 <- mk_dir("miss_perm", drop_perm)
res_m0 <- chk(dir_m0)

report(!inherits(res_m0, "try-error"), "a permuted set missing combinations still runs")
report(logged(paste(n_proda, "of", n_key, "unpermuted combination(s)"), since=m),
  "and says how many unpermuted combinations have no permuted counterpart")
report(logged("have no permuted counterpart", since=m),
  "and says which direction the mismatch is in")
report(!logged("permuted combination(s) are absent from", since=m),
  "and does not report a mismatch in the other direction, there being none")

if(!inherits(res_m0, "try-error")) {
  i_out <- res_m0$test %in% "proda"
  report(sum(i_out) %in% n_proda && all(is.na(res_m0$fdr[i_out])),
    "the unmatched combinations are the ones given no fdr")
  report(!any(is.na(res_m0$fdr[!i_out])), "and every matched combination still has one")
  report(min(which(i_out)) > max(which(!i_out)),
    "and the unmatched ones sort below every combination that has an fdr")
  ## matched by key rather than by position, the unmatched rows having moved:
  kb <- apply(base[, key_cols], 1, paste, collapse=":")
  km <- apply(res_m0[, key_cols], 1, paste, collapse=":")
  j <- match(km[!i_out], kb)
  report(!any(is.na(j)) && isTRUE(all.equal(res_m0[!i_out, out_cols], base[j, out_cols],
    check.attributes=FALSE)),
    "and dropping the permuted rows of one method leaves the other methods untouched")
}

## the other direction: drop a method from the unpermuted file only, so that permuted
##   rows are read, grouped and then thrown away with nothing said:

drop_unperm <- function(d, f) {
  if(!(f %in% unperm_file)) return(d)
  return(d[!(d$test %in% "deqms"), , drop=FALSE])
}

n_deqms <- sum(fix[[unperm_file]]$test %in% "deqms")

m <- mark()
dir_m1 <- mk_dir("miss_unperm", drop_unperm)
res_m1 <- chk(dir_m1)

report(!inherits(res_m1, "try-error"),
  "an unpermuted set missing combinations still runs")
report(logged(paste(n_deqms, "of", n_key, "permuted combination(s) are absent from"),
  since=m), "and says how many permuted combinations have no unpermuted counterpart")
report(logged(unperm_file, since=m), "and names the file they are absent from")
report(!logged("have no permuted counterpart", since=m),
  "and does not report a mismatch in the other direction, there being none")

if(!inherits(res_m1, "try-error")) {
  report(nrow(res_m1) %in% (nrow(fix[[unperm_file]]) - n_deqms),
    "the output covers the unpermuted combinations and no others")
  report(!any(res_m1$test %in% "deqms"),
    "and the permuted combinations with no unpermuted row are absent from it")
  report(!any(is.na(res_m1$fdr)),
    "and every combination that is there is matched, so all of them have an fdr")
}

## both at once, which is what an interrupted or a re-parameterised sweep looks like:

m <- mark()
dir_m2 <- mk_dir("miss_both", function(d, f) drop_perm(drop_unperm(d, f), f))
res_m2 <- chk(dir_m2)

report(!inherits(res_m2, "try-error"), "a mismatch in both directions still runs")
report(logged(paste(n_proda, "of", n_key - n_deqms, "unpermuted combination(s)"),
  since=m), "and counts the unpermuted side against the combinations actually present")
report(logged(paste(n_deqms, "of", n_key - n_proda,
  "permuted combination(s) are absent from"), since=m),
  "and counts the permuted side the same way")

## and the fixture as shipped, so that the warnings are known not to fire on a sweep
##   whose two sides agree:

m <- mark()
invisible(chk(dir_base))
report(!logged("have no permuted counterpart", since=m),
  "neither warning fires on the fixture as shipped")
report(!logged("permuted combination(s) are absent from", since=m),
  "in either direction")

###############################################################################
section("item 1: permuted runs that tested nothing are not counted as zero hits")

## tune() writes a row for a combination it skipped, and for one that lost every feature
##   to filtering, with ntests 0 and nhits NA. On the permuted side that NA became 0 and
##   was summarized beside the runs that ran, which says the permutation looked for false
##   positives and found none when it never looked: a combination skipped in every
##   permuted file came out with max1 0 and an fdr of 0, the best score in the table, and
##   sorted near the top of it. The shipped fixture has no such row, so every row below is
##   made into one here. The 30 msqrob rows it does have, with nhits NA and ntests 30, are
##   a different thing and are deliberately left alone by this: a test that ran.

kk0 <- apply(fix[[unperm_file]][, key_cols], 1, paste, collapse=":")
skip_keys <- head(kk0[!is.na(fix[[unperm_file]]$nhits)], 7)

blank <- function(d) {
  i <- apply(d[, key_cols], 1, paste, collapse=":") %in% skip_keys
  d$ntests[i] <- 0
  d$nhits[i] <- NA
  return(d)
}

## skipped in every permuted file, which is the case that scored best of all:

m <- mark()
dir_s0 <- mk_dir("skip_all", function(d, f) if(f %in% unperm_file) d else blank(d))
res_s0 <- chk(dir_s0)

report(!inherits(res_s0, "try-error"),
  "permuted rows that tested nothing do not stop the check")
report(logged(paste(length(skip_keys) * length(perm_files), "of",
  length(perm_files) * nrow(fix[[unperm_file]]), "permuted result row(s) have ntests 0"),
  since=m), "and are counted in a warning, against the number of permuted rows read")
report(logged("dropped rather than counted as permutations that found no hits", since=m),
  "which says what was done with them")

if(!inherits(res_s0, "try-error")) {
  ks0 <- apply(res_s0[, key_cols], 1, paste, collapse=":")
  i_skip <- ks0 %in% skip_keys
  report(sum(i_skip) %in% length(skip_keys),
    "the combinations they belong to are still in the output")
  report(all(is.na(res_s0$fdr[i_skip])),
    "and get no fdr, rather than an fdr of 0 from permutations that never ran")
  report(all(is.na(res_s0$max1[i_skip])) && all(is.na(res_s0$avg1[i_skip])),
    "and no max1 or avg1 either, there being nothing to take them over")
  n_cut <- sum(!is.na(res_s0$fdr) & res_s0$fdr < 0.05)
  report(all(which(i_skip) > n_cut),
    "and sort below every combination that met the fdr cutoff, not at the head of them")
  report(logged("have no permuted counterpart that tested anything", since=m),
    "and the mismatch warning says they have no permuted counterpart that ran")
  ## everything else is untouched, matched by key rather than position:
  kb <- apply(base[, key_cols], 1, paste, collapse=":")
  j <- match(ks0[!i_skip], kb)
  report(!any(is.na(j)) && isTRUE(all.equal(res_s0[!i_skip, out_cols], base[j, out_cols],
    check.attributes=FALSE)),
    "and blanking some permuted rows leaves the other combinations untouched")
}

## skipped in all but one permuted file, which is the quieter form: the runs that never
##   happened used to enter mid1, avg1 and sd1 as permutations that found no hits, pulling
##   all three toward 0. The one surviving file is given a distinctive hit count so that
##   an average over five files is telling apart from an average over one:

keep_file <- perm_files[1]
keep_hits <- 4

m <- mark()
dir_s1 <- mk_dir("skip_some", function(d, f) {
  if(f %in% unperm_file) return(d)
  if(!(f %in% keep_file)) return(blank(d))
  i <- apply(d[, key_cols], 1, paste, collapse=":") %in% skip_keys
  d$nhits[i] <- keep_hits
  return(d)
})
res_s1 <- chk(dir_s1)

report(!inherits(res_s1, "try-error"), "a partial skip still runs")
report(logged(paste(length(skip_keys) * (length(perm_files) - 1), "of"), since=m),
  "and counts only the rows that tested nothing")

if(!inherits(res_s1, "try-error")) {
  i_skip <- apply(res_s1[, key_cols], 1, paste, collapse=":") %in% skip_keys
  report(all(res_s1$max1[i_skip] %in% keep_hits),
    "the summaries are taken over the permuted runs that ran")
  report(all(res_s1$avg1[i_skip] %in% keep_hits),
    "so avg1 is the average over those, not over those plus the ones that did not")
  report(all(res_s1$mid1[i_skip] %in% keep_hits), "and mid1 likewise")
  report(all(is.na(res_s1$sd1[i_skip])),
    "and sd1 is NA, one run being one number, rather than an sd of mostly zeroes")
  report(!any(is.na(res_s1$fdr[i_skip])),
    "and the combinations still get an fdr, one permutation being enough for one")
}

## and the degenerate case, where no permuted row anywhere tested anything:

m <- mark()
dir_s2 <- mk_dir("skip_every", function(d, f) {
  if(f %in% unperm_file) return(d)
  d$ntests <- 0
  d$nhits <- NA
  return(d)
})
report(inherits(chk(dir_s2), "try-error"),
  "permuted results in which nothing at all was tested are refused")
report(logged("no permuted result row has ntests above 0", since=m),
  "with a message saying there is nothing to estimate an fdr from")

## the fixture as shipped has no such row, so none of this fires on it:

m <- mark()
invisible(chk(dir_base))
report(!logged("have ntests 0", since=m),
  "and none of the shipped permuted rows is dropped as untested")

###############################################################################
section("item 2: prefix picks the files, so two sweeps in one directory stay apart")

## prefix was taken and then ignored: the file pattern was built from suffix alone, so
##   every file in dir_in ending in suffix was read as part of this sweep. A second sweep
##   writing to the same directory with the same suffix had its files pooled into this
##   one's null without a word. Here the fixture is written out under prefix "a" beside
##   one file from another sweep under prefix "b", whose hit counts are far larger than
##   anything in the fixture, so pooling it is visible in every fdr:

other_hits <- 999

mk_dir_pfx <- function(tag) {
  dir_i <- file.path(tempdir(), paste0("tune_check_", tag))
  unlink(dir_i, recursive=TRUE)
  dir.create(dir_i, showWarnings=FALSE, recursive=TRUE)
  for(f in names(fix)) {
    utils::write.table(fix[[f]], file=file.path(dir_i, paste0("a", f)), sep="\t",
      quote=FALSE, row.names=FALSE)
  }
  d <- fix[[unperm_file]]
  d$nhits <- other_hits
  utils::write.table(d, file=file.path(dir_i, paste0("b", unperm_file)), sep="\t",
    quote=FALSE, row.names=FALSE)
  return(dir_i)
}

m <- mark()
dir_p <- mk_dir_pfx("prefix")
res_p <- try(suppressMessages(tune_check(dir_p, "a", suffix, cfg)), silent=TRUE)

report(!inherits(res_p, "try-error"), "a sweep reads its own files when prefixed")
report(!inherits(res_p, "try-error") && !inherits(base, "try-error") &&
  identical(res_p, base),
  "and gets the same answer it gets from those files unprefixed and alone")
report(!inherits(res_p, "try-error") && !any(res_p$max1 %in% other_hits),
  "so the other sweep's hit counts are nowhere in the summaries")
report(logged("do not begin with the prefix 'a'", since=m),
  "and the file passed over is reported rather than silently skipped")
report(logged(paste("1 file(s) in", dir_p), since=m), "with a count and the directory")
report(logged(paste0("b", unperm_file), since=m), "and the name of the file itself")

## and from the other side: prefix "b" is one file, an unpermuted one, so the sweep has
##   no permutations of its own. It used to find five, belonging to prefix "a":

m <- mark()
res_pb <- try(suppressMessages(tune_check(dir_p, "b", suffix, cfg)), silent=TRUE)

report(inherits(res_pb, "try-error"),
  "a prefix with no permuted files of its own does not borrow another sweep's")
report(logged("found 1 unpermuted file and 0 permuted files", since=m),
  "and says it found none, rather than the other prefix's five")
report(logged("no permuted results", since=m),
  "and is refused for having no null rather than for a downstream symptom of it")

## the unprefixed fixture, so the warning is known not to fire when there is nothing to
##   pass over:

m <- mark()
invisible(chk(dir_base))
report(!logged("do not begin with the prefix", since=m),
  "no files are reported as passed over when every file belongs to the sweep")

###############################################################################
section("item 3: a directory with no permuted file is refused where it is noticed")

## the count of permuted files was reported and then not acted on. With none of them,
##   do.call(rbind, list()) is NULL, the ntests fill turns that into a one-element list,
##   and the drop of rows that tested nothing dies on it with "incorrect number of
##   dimensions"; before that drop existed the same case reached the join column check and
##   came back as a missing column. Both name something other than the missing files.
##   f.err() stops with "Stopping" and puts its message in the log, so the message is
##   checked there and the absence of the R level symptom in the condition itself:

dir_np <- file.path(tempdir(), "tune_check_no_perm")
unlink(dir_np, recursive=TRUE)
dir.create(dir_np, showWarnings=FALSE, recursive=TRUE)
utils::write.table(fix[[unperm_file]], file=file.path(dir_np, unperm_file), sep="\t",
  quote=FALSE, row.names=FALSE)

m <- mark()
res_np <- chk(dir_np)

report(inherits(res_np, "try-error"),
  "a directory holding only the unpermuted file is refused")
report(logged("no permuted results", since=m),
  "for having no permuted results to estimate an fdr from")
report(logged(unperm_file, since=m) && logged(dir_np, since=m),
  "naming the unpermuted file it did find and the directory it looked in")
report(logged(paste0("ending in '", suffix, "'"), since=m),
  "and the suffix it looked for, that being what a typo lands in")
report(inherits(res_np, "try-error") &&
  !grepl("number of dimensions", as.character(res_np), fixed=TRUE),
  "rather than dying downstream on an empty table with an R level message")

###############################################################################
for(tag in c("base", "insert", "reorder", "bad_unperm", "bad_perm", "bad_both",
    "miss_perm", "miss_unperm", "miss_both", "skip_all", "skip_some", "skip_every",
    "prefix", "no_perm")) {
  unlink(file.path(tempdir(), paste0("tune_check_", tag)), recursive=TRUE)
}

cat("\n", strrep("=", 70), "\n", sep="")
cat("\n## passes:", n_pass, "; failures:", n_fail, "; elapsed:",
  round(as.numeric(difftime(Sys.time(), t0, units="secs"))), "s\n")

quit(save="no", status=if(n_fail > 0) 1 else 0)
