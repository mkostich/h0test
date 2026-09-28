## Tabulate what each test_*() method reports for a config$contrast run, over the shapes of
##   contrast that config$contrast accepts, and check the identities that must hold between
##   certain of those runs. Companion to capability_matrix.R, which does the same for
##   config$test_term; same fixture, so the two tables are comparable side by side.

usage <- function(msg=NULL) {
  if(!is.null(msg)) cat("ERROR:", msg, "\n\n", file=stderr())
  cat(
    "Tabulate what each h0testr test_*() method reports for a config$contrast run, over",
    "pairwise, averaged, scaled, reversed, intercept, interaction, no-intercept and",
    "non-marginal contrasts, and verify the identities that must hold between runs.",
    "",
    "Usage: Rscript contrast_matrix.R <r_dir> [--outfile=<path>]",
    "",
    "Required positional arguments:",
    "  <r_dir>          Path to the h0testr package R/ source directory; all .R",
    "                     files there are sourced (the installed package is not used).",
    "",
    "Optional named arguments:",
    "  --outfile=<path> Path to write the table to, prose wrapped to 100 characters.",
    "                     Printed to stdout in either case.",
    "",
    "Output: the contrast table to stdout (and to --outfile if given), followed by the",
    "  identity checks and the first line of every refusal message. Progress lines go to",
    "  stdout; expected-error messages go to contrast_matrix.err.log.",
    "",
    "Exit codes: 0 success; 1 a case could not be set up upstream; 2 usage error.",
    "",
    "Examples:",
    "  Rscript contrast_matrix.R C:/path/to/h0test/h0testr/R",
    "  Rscript contrast_matrix.R ../../h0test/h0testr/R --outfile=contrasts.txt",
    "  Rscript contrast_matrix.R C:/path/to/h0testr/R > contrasts.log 2>&1",
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

log_file <- "contrast_matrix.err.log"
cat("", file=log_file)
t0 <- Sys.time()

## f.err() writes the informative message to config$log_file and then stops with a bare
##   "Stopping", so the message has to be read back out of the log; f.msg() writes there
##   too, which is how the marginality warning is detected below:

log_len <- function() length(readLines(log_file, warn=FALSE))

log_slice <- function(n0) {
  txt <- readLines(log_file, warn=FALSE)
  if(length(txt) > n0) return(txt[-seq_len(n0)])
  return(character(0))
}

log_since <- function(n0, max_chars=400) {
  txt <- trimws(log_slice(n0))
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
## fixture: identical to capability_matrix.R's, so a cell here can be read against the
##   corresponding cell there. 24 observations; grp three levels, sex and batch two each,
##   age continuous. Forty genes with one to six peptides each, values already on a log
##   scale. The counts have to span more than a couple of values, or DEqMS's variance
##   prior, a loess of the log residual variance on the log peptide count, is singular:

set.seed(101)

nobs <- 24
ngene <- 40
npep_max <- 6

## fully crossed, so that even ~grp*sex*batch is of full rank and no cell of the design
##   is empty; a contrast within an interaction needs every cell it names to be occupied:

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
cfg0$test_prior_df <- 5
cfg0$frm <- ~grp + sex + batch + age
cfg0$test_term <- "sex"
cfg0$contrast <- ""

state_raw <- list(expression=exprs, features=feats, samples=samps)

###############################################################################
## cases: one per shape of contrast that config$contrast accepts, each probing a
##   different mechanism rather than a different formula. f.contrast_vector() binds each
##   coefficient name literally and lets R's own arithmetic build the weights, so the
##   expression has to name each coefficient the way R names it, and a coefficient whose
##   name is not a syntactic R name has to be backquoted. That covers both the intercept,
##   whose name is "(Intercept)", and every interaction coefficient, whose name contains
##   a colon: written bare, grpb:sexM parses as the ':' operator applied to two symbols
##   and is refused, since ':' is not one of the operators a weighted sum may contain.

cases <- list(
  list(frm=~sex,        con="sexM",                  lab="one coefficient"),
  list(frm=~grp,        con="grpb",                  lab="one coefficient of 3"),
  list(frm=~grp,        con="grpb - grpc",           lab="pairwise, non-reference"),
  list(frm=~grp,        con="(grpb + grpc)/2",       lab="two levels averaged"),
  list(frm=~grp,        con="2*grpb",                lab="scaled by a constant"),
  list(frm=~grp,        con="grpc - grpb",           lab="case 3 reversed"),
  list(frm=~sex + age,  con="age",                   lab="continuous coefficient"),
  list(frm=~0 + sex,    con="sexM - sexF",           lab="no intercept"),
  list(frm=~grp,        con="`(Intercept)` + grpb",  lab="a cell mean"),
  list(frm=~grp * sex,  con="`grpb:sexM` - `grpc:sexM`",
    lab="within an interaction"),
  list(frm=~grp * sex,  con="grpb - grpc",           lab="non-marginal"),
  list(frm=~age * sex,  con="age + `age:sexM`",      lab="slope in one group"),
  list(frm=~sex * batch, con="sexM + `sexM:batchb2`", lab="effect at one level")
)

## the reference run every identity below is measured against: the same hypothesis as
##   case 1, asked the way config$test_term asks it. Run through the same code path as a
##   case, but not tabulated, since capability_matrix.R already reports it:

ref_case <- list(frm=~sex, term="sex", lab="test_term reference")

methods <- c("lm", "trend", "voom", "deqms", "msqrob", "msqrob_agg", "proda",
  "prolfqua", "prolfqua_lmer")

## deqms and msqrob take peptide-level input, since both moderate on peptide counts, and
##   prolfqua_lmer and msqrob_agg take it because they model the peptides of a gene as a
##   random effect; the rest are run on the gene-level aggregate, as tune() runs them. The
##   predicate is the package's own, so this cannot drift from tune(). test_h0() rather than
##   test_*() is called, because the identities below are between the standardized pval and
##   logfc columns, which is where the contrast's own effect size is reported:

run_test <- function(nm, state, cfg) {
  if(f.gene_level_method(nm)) {
    st <- state
    cf <- cfg
  } else {
    agg <- suppressMessages(combine_features(state, cfg))
    st <- agg$state
    cf <- agg$config
  }
  return(try(suppressMessages(test_h0(st, cf, method=nm, is_log_transformed=TRUE,
    prior_df=5)), silent=TRUE))
}

## what the method reported. Every config$contrast run is one degree of freedom, whatever
##   config$frm looks like, so the question is which one-df statistic each engine reports:
##   't' a single-coefficient test with a signed effect, or 'F1' a nested-model F on one
##   degree of freedom, which is what test_lm() reports for any hypothesis, and what
##   test_prolfqua() reports here. test_lm() reported the likelihood ratio chi-square of
##   that same comparison until it was replaced by the exact F:

describe <- function(nm, res) {

  h <- res$original
  if(is.null(h) || !nrow(h)) return("empty")

  if(nm %in% c("lm", "prolfqua", "prolfqua_lmer")) return("F1")

  return(if(any(c("logFC", "diff") %in% names(h))) "t" else "F1")
}

## The answer computed outside h0testr, so that a row is anchored to something rather
##   than only to its neighbours. Same algebra and the same tolerances as
##   1/test_contrast.R's ind_effect() and ind_nested_p(), which do this for one contrast;
##   here it is done for every case. The identities further down are invariant to a
##   misaligned weight vector, since scaling a permuted L still scales the effect and
##   negating it still negates it, so this is the check that a contrast was assembled over
##   the coefficients it names:

ind_effect <- function(mat, X, L) {
  return(apply(mat, 1, function(y) {
    b <- stats::lm.fit(X, y)$coefficients
    b[is.na(b)] <- 0
    return(sum(L * b))
  }))
}

ind_nested_p <- function(mat, X, X_red) {
  df1 <- qr(X)$rank - qr(X_red)$rank
  df2 <- ncol(mat) - qr(X)$rank
  return(apply(mat, 1, function(y) {
    rss1 <- sum(stats::lm.fit(X, y)$residuals^2)
    rss0 <- sum(stats::lm.fit(X_red, y)$residuals^2)
    return(stats::pf(((rss0 - rss1) / df1) / (rss1 / df2), df1, df2,
      lower.tail=FALSE))
  }))
}

## Which methods must reproduce the least squares weighted sum exactly. The other five
##   fit something else to something else and are reported rather than checked: voom
##   reweights the observations, proDA models dropout with its own regularization, msqrob2
##   fits robustly, and msqrob_agg and prolfqua_lmer fit a mixed model to the peptides
##   instead of to the aggregate this reference is computed on. deqms belongs in the exact
##   group although it takes peptide-level input, because it aggregates internally with
##   combine_features(), which is the same function and configuration the reference uses:

exact_ref <- c("lm", "trend", "deqms", "prolfqua")

pad <- function(x, w, right=FALSE) {
  x <- as.character(x)
  if(right) return(formatC(x, width=w, flag=" "))
  return(formatC(x, width=-w, flag=" "))
}

## relative comparison, since a p-value can be arbitrarily small; the effect sizes are on
##   a log2 scale of order 1, so those are compared against the larger of the two and 1:

rel_diff <- function(a, b) {
  i <- intersect(names(a), names(b))
  if(!length(i)) return(NA_real_)
  a <- a[i]
  b <- b[i]
  ok <- is.finite(a) & is.finite(b)
  if(!any(ok)) return(NA_real_)
  return(max(abs(a[ok] - b[ok]) / pmax(abs(b[ok]), 1e-300)))
}

cor_with <- function(a, b) {
  i <- intersect(names(a), names(b))
  if(length(i) < 3) return(NA_real_)
  a <- a[i]
  b <- b[i]
  ok <- is.finite(a) & is.finite(b)
  if(sum(ok) < 3 || stats::sd(a[ok]) %in% 0 || stats::sd(b[ok]) %in% 0) {
    return(NA_real_)
  }
  return(stats::cor(a[ok], b[ok]))
}

abs_diff <- function(a, b) {
  i <- intersect(names(a), names(b))
  if(!length(i)) return(NA_real_)
  a <- a[i]
  b <- b[i]
  ok <- is.finite(a) & is.finite(b)
  if(!any(ok)) return(NA_real_)
  return(max(abs(a[ok] - b[ok])) / max(1, max(abs(b[ok]))))
}

## the standardized p-value and effect size, keyed by feature, for the identity checks:

standard_of <- function(res, col) {
  tbl <- res$standard
  if(is.null(tbl) || !nrow(tbl) || !all(c("feature", col) %in% names(tbl))) return(NULL)
  out <- tbl[[col]]
  names(out) <- as.character(tbl$feature)
  return(out)
}

###############################################################################
## run: one row per case, one cell per method, plus the reference run:

pvals <- list()          ## pvals[[case key]][[method]], named by feature
lfcs <- list()
refs <- list()           ## refs[[case key]] = the least squares effect and nested F p
unmod <- list()          ## unmod[[case key]] = prolfqua with moderation off
tab <- NULL
notes <- character(0)

setup <- function(cfg, what) {

  ## upstream: init_state() classifies the covariates and resolves factor levels, the two
  ##   filters screen features against this same formula, and f.design_test_cols() is the
  ##   shared derivation, dispatching to f.design_contrast() when config$contrast is set:

  n0 <- log_len()
  ini <- try(suppressMessages(init_state(state_raw, cfg, minimal=TRUE)), silent=TRUE)

  if(inherits(ini, "try-error")) {
    return(list(err=paste(what, "init_state()", why(ini, n0), sep="\t")))
  }

  st <- ini$state
  cf <- ini$config

  up <- c(
    frm=class(try(suppressMessages(filter_features_by_formula(st, cf)), silent=TRUE))[1],
    est=class(try(suppressMessages(filter_features_by_estimability(st, cf)),
      silent=TRUE))[1]
  )

  n0 <- log_len()
  des <- try(suppressMessages(f.design_test_cols(st, cf)), silent=TRUE)

  if(inherits(des, "try-error")) {
    return(list(err=paste(what, "f.design_test_cols()", why(des, n0), sep="\t")))
  }

  ## f.warn_contrast_marginality() reports through f.msg(), which appends to
  ##   config$log_file, so whether it fired is read back out of the lines this call added:

  marg <- any(grepl("f.warn_contrast_marginality|WARNING: f.design_contrast",
    log_slice(n0)))

  return(list(state=st, config=cf, design=des, up_ok=!any(up %in% "try-error"),
    marg=marg))
}

## the reference run first, so a failure of it is reported before 13 cases of comparisons
##   against nothing:

cfg <- cfg0
cfg$frm <- ref_case$frm
cfg$test_term <- ref_case$term
cfg$contrast <- ""
cfg$reference_levels <- cfg0$reference_levels[names(cfg0$reference_levels) %in%
  all.vars(ref_case$frm)]

progress("reference:", deparse(ref_case$frm), "/ test_term", ref_case$term)

su <- setup(cfg, "ref")
if(!is.null(su$err)) {
  cat("ERROR: the test_term reference run could not be set up:", su$err, "\n",
    file=stderr())
  quit(status=1)
}

pvals[["ref"]] <- list()
lfcs[["ref"]] <- list()

for(nm in methods) {
  n0 <- log_len()
  res <- run_test(nm, su$state, su$config)
  if(inherits(res, "try-error")) {
    notes <- c(notes, paste("ref", nm, why(res, n0), sep="\t"))
  } else {
    pvals[["ref"]][[nm]] <- standard_of(res, "pval")
    lfcs[["ref"]][[nm]] <- standard_of(res, "logfc")
  }
}

for(k in seq_along(cases)) {

  cs <- cases[[k]]
  cfg <- cfg0
  cfg$frm <- cs$frm
  cfg$test_term <- ""                    ## mutually exclusive with config$contrast
  cfg$contrast <- cs$con

  ## only the factor variables this formula uses: init_state() rejects a reference level
  ##   declared for a variable that is not in config$frm:

  cfg$reference_levels <- cfg0$reference_levels[names(cfg0$reference_levels) %in%
    all.vars(cs$frm)]

  progress("case", k, "of", length(cases), ":", deparse(cs$frm), "/", cs$con)

  su <- setup(cfg, k)
  if(!is.null(su$err)) {
    notes <- c(notes, su$err)
    next
  }

  key <- as.character(k)
  pvals[[key]] <- list()
  lfcs[[key]] <- list()

  ## the two references, on the gene-level aggregate, which is the matrix every row-wise
  ##   method is handed and the one deqms builds internally. The design is over the
  ##   observations, so it is the same matrix whether the features are peptides or genes:

  agg <- suppressMessages(combine_features(su$state, su$config))
  refs[[key]] <- list(
    eff=ind_effect(agg$state$expression, su$design$X, su$design$contrast),
    p=ind_nested_p(agg$state$expression, su$design$X, su$design$X_red)
  )

  ## test_prolfqua() with moderation and trend off is an exact nested F test, so it is the
  ##   one engine whose p-value can be held to the reference rather than only its effect
  ##   size; the moderated engines are checked on the effect size alone:

  cf_ex <- agg$config
  cf_ex$test_moderate <- FALSE
  cf_ex$test_trend <- FALSE
  n0 <- log_len()
  res_ex <- try(suppressMessages(test_h0(agg$state, cf_ex, method="prolfqua",
    is_log_transformed=TRUE)), silent=TRUE)

  if(inherits(res_ex, "try-error")) {
    notes <- c(notes, paste(k, "prolfqua (unmoderated reference)", why(res_ex, n0),
      sep="\t"))
  } else {
    unmod[[key]] <- standard_of(res_ex, "pval")
  }

  row <- c(
    case=key,
    frm=gsub("[[:space:]]", "", paste(deparse(cs$frm), collapse="")),
    con=cs$con,
    what=cs$lab,
    wts=as.character(length(su$design$cols_test)),
    df=as.character(su$design$df_intend),
    marg=if(su$marg) "y" else "n",
    up=if(su$up_ok) "ok" else "ERR"
  )

  for(nm in methods) {
    n0 <- log_len()
    res <- run_test(nm, su$state, su$config)
    if(inherits(res, "try-error")) {
      row[nm] <- "no"
      notes <- c(notes, paste(k, nm, why(res, n0), sep="\t"))
    } else {
      row[nm] <- describe(nm, res)
      pvals[[key]][[nm]] <- standard_of(res, "pval")
      lfcs[[key]][[nm]] <- standard_of(res, "logfc")
    }
  }

  ## the two reference columns. ref is the worst deviation of a reported effect size from
  ##   the least squares weighted sum, over the methods that must reproduce it exactly;
  ##   cor is the weakest correlation with it over all nine, which is what catches an
  ##   engine assembling a different contrast from the rest even where its coefficients
  ##   are legitimately its own:

  eff <- refs[[key]]$eff
  d_ex <- vapply(intersect(exact_ref, names(lfcs[[key]])),
    function(nm) abs_diff(lfcs[[key]][[nm]], eff), numeric(1))
  cors <- vapply(names(lfcs[[key]]),
    function(nm) cor_with(lfcs[[key]][[nm]], eff), numeric(1))

  row["ref"] <- if(!length(d_ex)) "NA" else
    paste0(if(max(d_ex, na.rm=TRUE) > 1e-6) "*" else "",
      formatC(max(d_ex, na.rm=TRUE), format="e", digits=0))
  row["cor"] <- if(!length(cors)) "NA" else
    formatC(min(cors, na.rm=TRUE), format="f", digits=3)

  tab <- rbind(tab, row)
}

tab <- as.data.frame(tab, stringsAsFactors=FALSE)
rownames(tab) <- NULL

###############################################################################
## identities. Each is exact arithmetic on the same fit, not a modelling assumption, so a
##   failure is a bug rather than a difference between engines. The last is exact only in
##   the same sense that a re-parameterization is: the two designs span the same space, so
##   an engine fitting the same model to the same values must reach the same p-value, but
##   it arrives there through a different optimization and gets a looser tolerance:

##   $expect names the methods whose deviation is a documented property of the engine's
##   route rather than a defect, reported separately and not marked in the table:

proda_lr <- paste("test_proda() sends any config$contrast to proDA's likelihood ratio",
  "branch on the constrained design, deliberately, so that proDA never re-parses the",
  "contrast; a one-column config$test_term run takes proDA's own Wald contrast instead.",
  "So the two are LR and Wald on the same fit. See test_proda.R, test_proda()")

ids <- list(
  list(a="1",  b="ref", fc=function(x) x,    tol=1e-9,
    lab="a contrast naming one coefficient is the test_term test of that term",
    expect=c(proda=proda_lr)),
  list(a="5",  b="2",   fc=function(x) 2 * x, tol=1e-9,
    lab="scaling a contrast leaves the p-value and scales the effect (case 5 vs 2)",
    expect=character(0)),
  list(a="6",  b="3",   fc=function(x) -x,   tol=1e-9,
    lab="reversing a contrast leaves the p-value and negates the effect (case 6 vs 3)",
    expect=character(0)),

  ## the two designs span the same space, so this is a re-parameterization: an engine
  ##   fitting the same model to the same values must reach the same p-value. The three
  ##   identities above are exact arithmetic on one fit and hold to floating point; this
  ##   one is two fits reached by two optimizations, so it is given a tolerance at the
  ##   level those converge to rather than at floating point:

  list(a="8",  b="ref", fc=function(x) x,    tol=1e-4,
    lab="~0 + sex states the same hypothesis as test_term sex on ~sex (case 8 vs ref)",
    expect=c(proda=proda_lr))
)

id_fail <- list()             ## id_fail[[case key]] = methods whose identity failed

## the measured deviation, whatever the verdict: an identity that misses by 5e-06 and one
##   that misses by 3e-02 are different findings, and reporting only pass or fail against a
##   threshold hides which is which:

dev_id <- function(id, nm) {
  pa <- pvals[[id$a]][[nm]]
  pb <- pvals[[id$b]][[nm]]
  la <- lfcs[[id$a]][[nm]]
  lb <- lfcs[[id$b]][[nm]]

  return(list(
    p=if(is.null(pa) || is.null(pb)) NA_real_ else rel_diff(pa, pb),
    l=if(is.null(la) || is.null(lb)) NA_real_ else abs_diff(la, id$fc(lb)),
    run=!is.null(pa) && !is.null(pb)
  ))
}

check_id <- function(id, nm) {
  d <- dev_id(id, nm)
  if(!d$run) return("not run")

  bad <- character(0)
  if(is.na(d$p)) bad <- c(bad, "p not comparable") else
    if(d$p > id$tol) bad <- c(bad, paste0("p rel ", signif(d$p, 3)))
  if(is.na(d$l)) bad <- c(bad, "effect not reported") else
    if(d$l > id$tol) bad <- c(bad, paste0("effect ", signif(d$l, 3)))

  if(!length(bad)) return("ok")
  return(paste(bad, collapse=", "))
}

id_lines <- character(0)

for(id in ids) {
  if(is.null(pvals[[id$a]])) {
    id_lines <- c(id_lines, paste0("  case ", id$a, " did not run; ", id$lab), "")
    next
  }
  verdict <- vapply(methods, function(nm) check_id(id, nm), character(1))
  skip <- names(verdict)[verdict %in% "not run"]
  bad <- names(verdict)[!verdict %in% c("ok", "not run")]
  known <- intersect(bad, names(id$expect))
  bad <- setdiff(bad, known)

  id_lines <- c(id_lines, paste0("  ", id$lab, ":"),
    paste0("    holds to ", format(id$tol, scientific=TRUE), " for ",
      sum(verdict %in% "ok"), " of ", length(methods), " methods",
      if(length(skip)) paste0("; not run: ", paste(skip, collapse=", ")) else "", "."))

  ## the deviation of every method that ran both, largest over the features, relative for
  ##   the p-value and scaled by the larger effect for the effect size:

  devs <- vapply(setdiff(methods, skip), function(nm) {
    d <- dev_id(id, nm)
    return(paste0(nm, " ", signif(d$p, 2), "/", signif(d$l, 2)))
  }, character(1))

  id_lines <- c(id_lines,
    strwrap(paste("p/effect deviation:", paste(devs, collapse="; ")),
      width=98, indent=4, exdent=6))

  for(nm in bad) {
    id_lines <- c(id_lines, paste0("    FAIL ", nm, ": ", verdict[[nm]]))
    id_fail[[id$a]] <- unique(c(id_fail[[id$a]], nm))
  }

  for(nm in known) {
    id_lines <- c(id_lines,
      strwrap(paste0("expected, not a defect: ", nm, ": ", verdict[[nm]], ". ",
        id$expect[[nm]]), width=98, indent=4, exdent=6))
  }

  id_lines <- c(id_lines, "")
}

###############################################################################
## the external anchors, summarized over the cases: what each method's reported effect
##   size did against the least squares weighted sum, and what the one unmoderated engine's
##   p-value did against the exact nested F. Worst case over the thirteen, since a
##   per-case-per-method listing would be 117 lines saying "0" a hundred times:

worst <- function(nm, what) {
  keys <- names(refs)
  vals <- vapply(keys, function(k) {
    got <- lfcs[[k]][[nm]]
    if(is.null(got)) return(NA_real_)
    if(what %in% "dev") return(abs_diff(got, refs[[k]]$eff))
    return(cor_with(got, refs[[k]]$eff))
  }, numeric(1))
  if(all(is.na(vals))) return(list(v=NA_real_, k=NA_character_))
  i <- if(what %in% "dev") which.max(vals) else which.min(vals)
  return(list(v=vals[[i]], k=keys[[i]]))
}

ref_lines <- c(
  "Effect size against the least squares weighted sum L'beta, computed here from",
  "  stats::lm.fit() on the same aggregate, one value per gene per case. Worst of the 13",
  "  cases per method; dev is the largest absolute difference over the genes, scaled by",
  "  the larger effect and 1, and cor is the correlation over the genes:")

for(nm in methods) {
  d <- worst(nm, "dev")
  c1 <- worst(nm, "cor")
  ref_lines <- c(ref_lines, paste0("  ", pad(nm, 15),
    if(nm %in% exact_ref) "exact:    " else "reported: ",
    "dev ", if(is.na(d$v)) "NA" else formatC(d$v, format="e", digits=1),
    " (case ", d$k, ");  cor ",
    if(is.na(c1$v)) "NA" else formatC(c1$v, format="f", digits=4),
    " (case ", c1$k, ")"))
}

ref_bad <- vapply(intersect(exact_ref, methods),
  function(nm) isTRUE(worst(nm, "dev")$v > 1e-6), logical(1))

ref_lines <- c(ref_lines,
  if(any(ref_bad)) paste0("  FAIL, beyond 1e-06: ",
    paste(names(ref_bad)[ref_bad], collapse=", ")) else
    "  All four methods that must reproduce L'beta exactly do so, in every case.",
  "",
  "P-value against an exact nested F test of the same two designs, also computed here.",
  "  test_prolfqua() with config$test_moderate and config$test_trend off is that test, so",
  "  it is the one engine whose p-value can be anchored rather than only its effect size;",
  "  every other engine moderates, and its p-value is meant to differ:")

p_dev <- vapply(names(refs), function(k) {
  if(is.null(unmod[[k]])) return(NA_real_)
  return(rel_diff(unmod[[k]], refs[[k]]$p))
}, numeric(1))

ref_lines <- c(ref_lines,
  if(all(is.na(p_dev))) "  the unmoderated reference run did not complete for any case." else
    strwrap(paste0("prolfqua unmoderated: largest relative deviation ",
      formatC(max(p_dev, na.rm=TRUE), format="e", digits=1), " over the ",
      sum(!is.na(p_dev)), " cases it ran, worst at case ",
      names(p_dev)[which.max(p_dev)], "; ",
      if(max(p_dev, na.rm=TRUE) > 1e-8) "beyond 1e-08, so LISTED BELOW." else
        "so every case is the exact test of the contrast."),
      width=98, indent=2, exdent=4),
  if(isTRUE(max(p_dev, na.rm=TRUE) > 1e-8)) paste0("  case ", names(p_dev),
    ": ", formatC(p_dev, format="e", digits=1))[!is.na(p_dev) & p_dev > 1e-8] else NULL,
  "")

## a cell whose identity failed is marked in the table, so the table and the checks below
##   cannot be read apart from one another:

for(key in names(id_fail)) {
  i <- which(tab$case %in% key)
  if(length(i)) {
    for(nm in id_fail[[key]]) tab[[nm]][i] <- paste0(tab[[nm]][i], "*")
  }
}

###############################################################################
## render:

## every config$contrast run is one degree of freedom by construction, so the df column is
##   dropped and stated once in the prose; the same for the upstream column when every
##   case passed the filters. Both reappear as soon as a case differs:

df_all_1 <- all(tab$df %in% "1")
up_all_ok <- all(tab$up %in% "ok")
marg_any <- any(tab$marg %in% "y")

## the table itself is allowed past 100 characters: nine method columns do not fit
##   otherwise, and abbreviating the header costs more in legibility than the width saves.
##   The surrounding prose is still wrapped to 100.

mw <- pmax(7, nchar(methods) + 1)          ## keep a gap ahead of the widest name

hdr <- paste0(pad("case", 5), pad("formula", 12), pad("config$contrast", 28),
  pad("wts", 5, TRUE), if(df_all_1) "" else pad("df", 4, TRUE),
  if(marg_any) pad("marg", 6, TRUE) else "",
  if(up_all_ok) "" else pad("up", 5, TRUE),
  pad("ref", 8, TRUE), pad("cor", 7, TRUE), " ",
  paste(vapply(seq_along(methods), function(j) pad(methods[j], mw[j], TRUE),
    character(1)), collapse=""))

lines <- c(hdr, strrep("-", nchar(hdr)))

for(i in seq_len(nrow(tab))) {
  lines <- c(lines,
    paste0(pad(tab$case[i], 5), pad(tab$frm[i], 12), pad(tab$con[i], 28),
      pad(tab$wts[i], 5, TRUE), if(df_all_1) "" else pad(tab$df[i], 4, TRUE),
      if(marg_any) pad(tab$marg[i], 6, TRUE) else "",
      if(up_all_ok) "" else pad(tab$up[i], 5, TRUE),
      pad(tab$ref[i], 8, TRUE), pad(tab$cor[i], 7, TRUE), " ",
      paste(vapply(seq_along(methods), function(j) pad(tab[[methods[j]]][i], mw[j], TRUE),
        character(1)), collapse="")))
}

## what each case is for, since the contrast expression alone does not say:

what_lines <- unlist(lapply(seq_len(nrow(tab)), function(i)
  paste0("  ", pad(tab$case[i], 4), pad(tab$con[i], 28), tab$what[i])))

out <- c(
  "h0testr: what each test_*() method reports for a config$contrast run, and the",
  "  identities that hold between such runs. Same fixture as the config$test_term table:",
  "  24 observations, fully crossed, grp (3 levels), sex and batch (2 levels each), age",
  "  continuous; 40 genes, 1 to 6 peptides each, no missing values, already log scale.",
  "",
  "wts  = coefficients config$contrast gives a non-zero weight. It is not the degrees of",
  if(df_all_1) c(
  "         freedom of the test: a contrast is a statement that one weighted sum of",
  "         coefficients is zero, so every case below is a 1 df test, however many",
  "         coefficients it weights and whatever config$frm looks like.") else
  "         freedom of the test, which the df column gives.",
  if(marg_any) c(
  "marg = whether f.warn_contrast_marginality() warned that the contrast weights a",
  "         coefficient of a term that a higher-order term in config$frm also contains,",
  "         so the comparison is at the reference level of the other variable rather",
  "         than averaged over it.") else NULL,
  if(up_all_ok) c(
  "Every case below runs init_state(), filter_features_by_formula() and",
  "  filter_features_by_estimability() without error, so all of them are contrasts",
  "  h0testr accepts upstream of the test.") else
  "up   = whether init_state() and the two formula-aware filters ran (ok) or not.",
  "",
  "Cells: t   = single-coefficient test (moderated t or Wald), the contrast having been",
  "               made the one coefficient of the fit;",
  "       F1  = nested-model F on one degree of freedom, the reduced model being",
  "               config$frm constrained so that the contrast is zero; this is what",
  "               test_lm() reports for any hypothesis;",
  "       no  = refused with an informative error (quoted below);",
  "       *   = an identity that must hold for this cell does not; see the checks below.",
  "ref  = the largest deviation of a reported effect size from the least squares weighted",
  "         sum L'beta, computed outside h0testr, over the four methods that must",
  "         reproduce it exactly; * marks a deviation beyond 1e-06. See below for the",
  "         other five, which fit something else and are reported rather than checked.",
  "cor  = the weakest correlation with that same L'beta over all nine methods, reported and",
  "         not checked. It was added to catch a contrast assembled over the wrong",
  "         coefficients in an engine whose coefficients are legitimately its own, which the",
  "         identities below cannot catch, since scaling a misaligned weight vector still",
  "         scales the effect and negating it still negates it. On this fixture it does not",
  "         do that job: the values are simulated noise with no planted effect, so the five",
  "         engines that reweight, fit robustly, model dropout or fit the peptides differ",
  "         from least squares by as much as the signal, and no threshold separates that",
  "         from a real defect. 1/test_contrast.R makes the same comparison work by planting",
  "         an effect of 3, which dominates the spread; planting one here would cost the",
  "         shared fixture that makes this table and the config$test_term one comparable.",
  "         The four exact methods are 1.0000 in every case, which is what the ref column",
  "         already says more sharply.",
  "The standardized logfc column holds the signed value of the contrast itself, the",
  "  weighted sum of the coefficients it names, in every cell that runs; unlike a joint",
  "  test of several coefficients, a contrast has a single signed effect size. See",
  "  f.logfc_effect().",
  "config$contrast names design matrix coefficients, so each name has to be written the",
  "  way R names it, and backquoted where that is not a syntactic R name: `(Intercept)`",
  "  and, for every interaction coefficient, `grpb:sexM`. Written bare, grpb:sexM parses",
  "  as the ':' operator rather than as a name and is refused, ':' not being one of the",
  "  operators a weighted sum may contain.",
  "",
  lines,
  "",
  "What each case is for:",
  what_lines,
  "",
  strwrap(paste0("Cases run, of ", nrow(tab), " accepted upstream: ",
    paste(vapply(methods, function(m) paste0(m, " ",
      sum(!tab[[m]] %in% "no")), character(1)), collapse="; "), "."),
    width=98, exdent=2),
  "",
  "Identities checked, over the features of every method that ran both runs. The",
  "  deviation reported for each method is the largest over the features: relative for the",
  "  p-value, and for the effect size the absolute difference scaled by the larger of the",
  "  two effects and 1. NA means the method reported no effect size for one of the runs.",
  id_lines,
  ref_lines,
  ""
)

if(length(notes)) {

  ## the same refusal recurs across cases, differing only in the term and coefficient
  ##   names it quotes, so group by the shape of the message:

  parts <- do.call(rbind, strsplit(notes, "\t", fixed=TRUE))
  shape <- gsub("[0-9]+", "N", gsub("\\([^)]*\\)", "(...)",
    gsub("'[^']*'", "'...'", parts[, 3])))
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
} else {
  out <- c(out, "No method refused any case above.", "")
}

cat(paste(out, collapse="\n"), "\n")

if(!is.null(outfile)) {
  writeLines(out, outfile)
  progress("wrote", outfile)
}

progress("done; expected-error messages in", log_file)

quit(status=0)
