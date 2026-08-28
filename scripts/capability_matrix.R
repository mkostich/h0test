## Tabulate which of the formula/test_term combinations that h0testr accepts upstream
##   (f.parse_frm(), init_state(), the filters) each test_*() method can actually run,
##   with both categorical and continuous covariates. One row per (frm, test_term) case,
##   one column per method.

usage <- function(msg=NULL) {
  if(!is.null(msg)) cat("ERROR:", msg, "\n\n", file=stderr())
  cat(
    "Tabulate the formula/test_term combinations accepted upstream by h0testr",
    "against what each test_*() method can run, over categorical (2 and 3 level),",
    "continuous, interaction, no-intercept and intercept-test cases.",
    "",
    "Usage: Rscript capability_matrix.R <r_dir> [--outfile=<path>]",
    "",
    "Required positional arguments:",
    "  <r_dir>          Path to the h0testr package R/ source directory; all .R",
    "                     files there are sourced (the installed package is not used).",
    "",
    "Optional named arguments:",
    "  --outfile=<path> Path to write the table to, wrapped to 100 characters.",
    "                     Printed to stdout in either case.",
    "",
    "Output: the capability table to stdout (and to --outfile if given), followed by",
    "  the first line of every refusal message. Progress lines go to stdout.",
    "",
    "Exit codes: 0 success; 1 a case could not be set up upstream; 2 usage error.",
    "",
    "Examples:",
    "  Rscript capability_matrix.R C:/path/to/h0test/h0testr/R",
    "  Rscript capability_matrix.R ../../h0test/h0testr/R --outfile=answer.txt",
    "  Rscript capability_matrix.R C:/path/to/h0testr/R > caps.log 2>&1",
    sep="\n", file=stderr()
  )
  quit(status=2)
}

args <- commandArgs(trailingOnly=TRUE)
if(length(args) < 1 || length(args) > 2) usage("wrong number of arguments")

r_dir <- args[1]
if(!dir.exists(r_dir)) usage(paste("r_dir not a directory:", r_dir))

outfile <- NULL
if(length(args) %in% 2) {
  if(!grepl("^--outfile=", args[2])) usage(paste("unrecognized argument:", args[2]))
  outfile <- sub("^--outfile=", "", args[2])
}

for(f in list.files(r_dir, pattern="[.]R$", full.names=TRUE)) source(f)

log_file <- "capability_matrix.err.log"
cat("", file=log_file)
t0 <- Sys.time()

## f.err() writes the informative message to config$log_file and then stops with a
##   bare "Stopping", so the message has to be read back out of the log:

log_len <- function() length(readLines(log_file, warn=FALSE))

log_since <- function(n0, max_chars=400) {
  txt <- readLines(log_file, warn=FALSE)
  if(length(txt) > n0) txt <- txt[-seq_len(n0)] else txt <- character(0)
  txt <- trimws(txt)
  txt <- txt[nzchar(txt)]
  if(!length(txt)) return("")

  ## the failing call may have logged a whole block first (test_deqms() runs
  ##   combine_features(), which reports the state), so take the error itself:

  i <- grep("ERROR", txt)
  if(length(i)) txt <- txt[i[length(i)]:length(txt)]

  txt <- paste(txt, collapse=" ")
  txt <- sub("[[:space:]]*at:[[:space:]]*[0-9].*$", "", txt)     ## trailing timestamp
  if(nchar(txt) > max_chars) txt <- paste0(substr(txt, 1, max_chars), " [...]")
  return(txt)
}

why <- function(res, n0) {
  msg <- sub("\n.*", "", conditionMessage(attr(res, "condition")))
  if(grepl("^Stopping", msg)) {
    alt <- log_since(n0)
    if(nzchar(alt)) return(alt)
  }
  return(msg)
}

progress <- function(...) {
  cat("[", round(as.numeric(difftime(Sys.time(), t0, units="secs"))), "s ]", ..., "\n")
  utils::flush.console()
}

###############################################################################
## fixture: 24 observations; grp three levels, sex and batch two each, age
##   continuous. Forty genes with one to six peptides each, values already on a log
##   scale. The counts have to span more than a couple of values: DEqMS fits its variance
##   prior as a loess of the log residual variance on the log peptide count, and with two
##   distinct counts, which is what this fixture had, that loess is singular and its
##   posterior variance comes back NaN for some genes. That went unnoticed while nothing
##   read those columns; test_method "deqms" reports them now:

set.seed(101)

nobs <- 24
ngene <- 40
npep_max <- 6

## fully crossed, so that even ~grp*sex*batch is of full rank (12 parameters, 12
##   residual degrees of freedom) and no cell of the design is empty; a layout that
##   confounds two of the factors makes the higher-order cases non-estimable, which
##   is a property of the data rather than of the methods:

samps <- expand.grid(rep=1:2, batch=c("b1", "b2"), sex=c("F", "M"),
  grp=c("a", "b", "c"), stringsAsFactors=FALSE)
samps <- data.frame(
  obs=paste0("o", sprintf("%02d", 1:nobs)),
  grp=samps$grp,
  sex=samps$sex,
  batch=samps$batch,
  age=round(rnorm(nobs, 50, 10), 1),
  stringsAsFactors=FALSE
)

npep <- ((seq_len(ngene) - 1) %% npep_max) + 1          ## 1..6 peptides per gene
gene <- rep(paste0("g", sprintf("%02d", seq_len(ngene))), times=npep)
exprs <- matrix(rnorm(length(gene) * nobs, 20, 1.5), nrow=length(gene))
rownames(exprs) <- paste0(gene, ".p", unlist(lapply(npep, seq_len)))
colnames(exprs) <- samps$obs
feats <- data.frame(pep=rownames(exprs), gene=gene)

cfg0 <- new_config()
cfg0$obs_id_col <- cfg0$obs_col <- "obs"
cfg0$sample_id_col <- "obs"
cfg0$feat_id_col <- cfg0$feat_col <- "pep"
cfg0$gene_id_col <- "gene"
cfg0$reference_levels <- c(grp="a", sex="F", batch="b1")
cfg0$estimability <- "test"
cfg0$df_resid_min <- 2
cfg0$save_state <- FALSE
cfg0$permute_var <- ""
cfg0$is_log_transformed <- TRUE
cfg0$impute_method <- "none"
cfg0$n_features_min <- 5
cfg0$log_file <- log_file
cfg0$frm <- ~grp + sex + batch + age
cfg0$test_term <- "sex"

state_raw <- list(expression=exprs, features=feats, samples=samps)

###############################################################################
## cases: every construct the upstream formula checker accepts (main effects,
##   interactions of any order, mixed factor/continuous, intercept suppression,
##   intercept as the tested term), with a factor of more than two levels and a
##   continuous covariate each appearing alone and inside an interaction:

cases <- list(
  list(frm=~sex,             term="sex",       lab="2-level factor"),
  list(frm=~grp,             term="grp",       lab="3-level factor"),
  list(frm=~age,             term="age",       lab="continuous"),
  list(frm=~sex + age,       term="age",       lab="continuous, factor adjusted"),
  list(frm=~sex + age,       term="sex",       lab="factor, continuous adjusted"),
  list(frm=~grp + age,       term="grp",       lab="3-level, continuous adjusted"),
  list(frm=~sex * batch,     term="sex",       lab="factor in interaction"),
  list(frm=~sex * batch,     term="sex:batch", lab="factor:factor term"),
  list(frm=~grp * sex,       term="grp",       lab="3-level in interaction"),
  list(frm=~grp * sex,       term="grp:sex",   lab="3-level:factor term"),
  list(frm=~age * sex,       term="age",       lab="continuous in interaction"),
  list(frm=~age * sex,       term="age:sex",   lab="continuous:factor term"),
  list(frm=~grp * sex * batch, term="grp:sex:batch", lab="3-way term"),
  list(frm=~grp * sex * batch, term="sex",     lab="factor in 3-way"),
  list(frm=~0 + sex,         term="sex",       lab="no intercept"),
  list(frm=~sex,             term="1",         lab="intercept tested")
)

methods <- c("lm", "trend", "voom", "deqms", "msqrob", "msqrob_agg", "proda",
  "prolfqua", "prolfqua_lmer")

## deqms and msqrob take peptide-level input, since both moderate on peptide
##   counts, and prolfqua_lmer and msqrob_agg take it because they model the peptides of
##   a gene as a random effect; the rest are run on the gene-level aggregate, as tune()
##   runs them. The predicate is the package's own, so this cannot drift from tune():

run_test <- function(nm, state, cfg) {
  fn <- get(paste0("test_", sub("_(lmer|agg)$", "", nm)))
  if(f.gene_level_method(nm)) {
    st <- state
    cf <- cfg
  } else {
    agg <- suppressMessages(combine_features(state, cfg))
    st <- agg$state
    cf <- agg$config
  }
  args <- list(st, cf)
  if(nm %in% "proda") args <- list(st, cf, is_log_transformed=TRUE, prior_df=5)
  if(nm %in% "prolfqua") args <- list(st, cf, is_log_transformed=TRUE)
  if(nm %in% "prolfqua_lmer") {
    args <- list(st, cf, is_log_transformed=TRUE, mixed=TRUE)
  }
  if(nm %in% "msqrob_agg") args <- list(st, cf, aggregate=TRUE)
  return(try(suppressMessages(do.call(fn, args)), silent=TRUE))
}

## what the method reported: 't' a single-coefficient test with a fold change,
##   'F' a joint test over several coefficients, 'F1' an F on one coefficient
##   (test_prolfqua() always reports F, being a nested model comparison):

describe <- function(nm, res, n_cols) {

  h <- res$hits
  if(is.null(h) || !nrow(h)) return("empty")

  ## test_lm() and test_prolfqua() both report the F of a nested model comparison,
  ##   whether one column carries the test or several; test_lm() reported the
  ##   likelihood ratio chi-square of that comparison until it was replaced by the
  ##   exact F. test_lm() reports the coefficients under test alongside it:

  if(nm %in% c("lm", "prolfqua", "prolfqua_lmer")) {
    return(if(n_cols > 1) "F" else "F1")
  }

  if(n_cols > 1) return("F")
  return(if(any(c("logFC", "diff") %in% names(h))) "t" else "F1")
}

###############################################################################
## one row per case, one cell per method:

tab <- NULL
notes <- character(0)

for(k in seq_along(cases)) {

  cs <- cases[[k]]
  cfg <- cfg0
  cfg$frm <- cs$frm
  cfg$test_term <- cs$term

  ## only the factor variables that this formula actually uses: init_state() rejects
  ##   a reference level declared for a variable that is not in config$frm:

  vars <- all.vars(cs$frm)
  cfg$reference_levels <- cfg0$reference_levels[names(cfg0$reference_levels) %in% vars]

  progress("case", k, "of", length(cases), ":", deparse(cs$frm), "/", cs$term)

  ## upstream: init_state() classifies the covariates and resolves factor levels,
  ##   the two filters screen features against this same formula, and
  ##   f.design_test_cols() is the shared derivation of the tested columns:

  n0 <- log_len()
  ini <- try(suppressMessages(init_state(state_raw, cfg, minimal=TRUE)), silent=TRUE)

  if(inherits(ini, "try-error")) {
    notes <- c(notes, paste(k, "init_state()", why(ini, n0), sep="\t"))
    next
  }

  st <- ini$state
  cf <- ini$config

  up <- c(
    frm=class(try(suppressMessages(filter_features_by_formula(st, cf)), silent=TRUE))[1],
    est=class(try(suppressMessages(filter_features_by_estimability(st, cf)),
      silent=TRUE))[1]
  )
  up_ok <- !any(up %in% "try-error")

  n0 <- log_len()
  des <- try(suppressMessages(f.design_test_cols(st, cf)), silent=TRUE)

  if(inherits(des, "try-error")) {
    notes <- c(notes, paste(k, "f.design_test_cols()", why(des, n0), sep="\t"))
    next
  }

  n_cols <- length(des$cols_test)

  row <- c(
    case=as.character(k),
    frm=gsub("[[:space:]]", "", paste(deparse(cs$frm), collapse="")),
    term=cs$term,
    what=cs$lab,
    cols=as.character(n_cols),
    df=as.character(des$df_intend),
    up=if(up_ok) "ok" else "ERR"
  )

  for(nm in methods) {
    n0 <- log_len()
    res <- run_test(nm, st, cf)
    if(inherits(res, "try-error")) {
      row[nm] <- "no"
      notes <- c(notes, paste(k, nm, why(res, n0), sep="\t"))
    } else {
      row[nm] <- describe(nm, res, n_cols)
    }
  }

  tab <- rbind(tab, row)
}

tab <- as.data.frame(tab, stringsAsFactors=FALSE)
rownames(tab) <- NULL

###############################################################################
## render, wrapped within 100 characters:

pad <- function(x, w, right=FALSE) {
  x <- as.character(x)
  if(right) return(formatC(x, width=w, flag=" "))
  return(formatC(x, width=-w, flag=" "))
}

## the numerator df of the test is the number of design columns carrying it in every
##   case here, so it is reported once rather than as a column of its own:

df_is_cols <- all(tab$df %in% tab$cols)
up_all_ok <- all(tab$up %in% "ok")

## the table itself is allowed past 100 characters: nine method columns do not fit
##   otherwise, and abbreviating the header costs more in legibility than the width
##   saves. The surrounding prose is still wrapped to 100.

mw <- pmax(8, nchar(methods) + 1)          ## keep a gap ahead of the widest name

## the df and up columns are dropped from the table when they say the same thing for
##   every case, which the prose above then states once, so the table stays inside
##   100 characters; they reappear as soon as a case differs:

cell <- function(i, nm) {
  txt <- tab[[nm]][i]
  if(nm %in% "cols" && !df_is_cols) txt <- paste0(txt, "/", tab$df[i])
  return(txt)
}

lines <- character(0)
lines <- c(lines,
  paste0(pad("case", 5), pad("formula", 15), pad("test_term", 14),
    pad(if(df_is_cols) "cols" else "cols/df", 8, TRUE),
    if(up_all_ok) "" else pad("up", 5, TRUE), " ",
    paste(vapply(seq_along(methods), function(j) pad(methods[j], mw[j], TRUE),
      character(1)), collapse="")))
lines <- c(lines, strrep("-", nchar(lines[1])))

for(i in seq_len(nrow(tab))) {
  lines <- c(lines,
    paste0(pad(tab$case[i], 5), pad(tab$frm[i], 15), pad(tab$term[i], 14),
      pad(cell(i, "cols"), 8, TRUE),
      if(up_all_ok) "" else pad(tab$up[i], 5, TRUE), " ",
      paste(vapply(seq_along(methods), function(j) pad(tab[[methods[j]]][i], mw[j], TRUE),
        character(1)), collapse="")))
}

out <- c(
  "h0testr: what each test_*() method can run, over the formula/test_term",
  "  combinations accepted upstream. Fully crossed fixture: 24 observations, grp (3",
  "  levels), sex and batch (2 levels each), age continuous; 40 genes, 1 to 6",
  "  peptides each, no missing values, values already on a log scale.",
  "",
  "cols = design matrix columns carrying the test of test_term, which is also the",
  "         numerator df of the test in every case below.",
  if(up_all_ok) c(
    "Every case below runs init_state(), filter_features_by_formula() and",
    "  filter_features_by_estimability() without error, so all of them are formulas",
    "  h0testr accepts upstream of the test.") else
    "up   = whether init_state() and the two formula-aware filters ran (ok) or not.",
  "",
  "Cells: t  = single-coefficient test (moderated t or Wald);",
  "       F  = joint test over all the columns carrying the test;",
  "       F1 = nested-model F on a single column, which is what test_lm() and",
  "              test_prolfqua() report there (the same p-value as the corresponding",
  "              moderated t);",
  "       no = refused with an informative error (quoted below).",
  "Every cell that runs carries an effect size in the standardized logfc column: the",
  "  signed coefficient where one column carries the test, and the total swing over the",
  "  tested columns where several do; see f.logfc_effect().",
  "",
  lines,
  "",
  strwrap(paste0("Cases run, of ", nrow(tab), " accepted upstream: ",
    paste(vapply(methods, function(m) paste0(m, " ", sum(!tab[[m]] %in% "no")),
      character(1)), collapse="; "), "."), width=98, exdent=2),
  ""
)

if(length(notes)) {

  ## the same refusal recurs across cases, differing only in the term and coefficient
  ##   names it quotes, so group by the shape of the message:

  parts <- do.call(rbind, strsplit(notes, "\t", fixed=TRUE))
  shape <- gsub("[0-9]+", "N", gsub("\\([^)]*\\)", "(...)",
    gsub("'[^']*'", "'...'", parts[, 3])))

  ## the coefficient-count refusal quotes the formula as well, so collapse it to one
  ##   entry per method rather than one per formula:

  shape[grepl("can test at most", parts[, 3])] <- "max_cols"
  key <- paste(parts[, 2], shape)

  out <- c(out,
    "Refusals, one entry per distinct reason, quoting the message from the first case",
    "  listed:", "")

  for(kk in unique(key)) {
    i <- which(key %in% kk)
    out <- c(out, strwrap(paste0(parts[i[1], 2], ", case",
      if(length(i) > 1) "s" else "", " ", paste(parts[i, 1], collapse=", "), ": ",
      parts[i[1], 3]), width=98, indent=2, exdent=6), "")
  }
}

cat(paste(out, collapse="\n"), "\n")

if(!is.null(outfile)) {
  writeLines(out, outfile)
  progress("wrote", outfile)
}

progress("done; expected-error messages in", log_file)

quit(status=0)
