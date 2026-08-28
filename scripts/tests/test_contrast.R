## Tests config$contrast: a weighted sum of the coefficients of config$frm, tested
##   within the full model, as an alternative to config$test_term. A contrast is one
##   degree of freedom however many coefficients it weights, which is why it reaches
##   engines a joint test does not, DEqMS above all. Two things are checked hardest:
##   that the weights the parser builds are the ones the expression says, and that the
##   hypothesis the engines run is the one those weights name, verified against an
##   independent nested F test and an independent computation of the weighted sum.

usage <- function(msg=NULL) {
  if(!is.null(msg)) cat("ERROR:", msg, "\n\n", file=stderr())
  cat(
    "Test config$contrast in h0testr: the parser that turns the expression into",
    "coefficient weights and its errors, the constrained design derived from those",
    "weights, the estimability check, the marginality warning, the refusal to set",
    "both config$contrast and config$test_term, per-feature estimability screening,",
    "and all seven test_method engines, including test_deqms(), which a contrast",
    "reaches as a single coefficient. Checks the reported effect size against an",
    "independent computation of the weighted sum of the coefficients, and the",
    "unmoderated p-values against an independent nested F test.",
    "",
    "Usage: Rscript test_contrast.R <r_dir>",
    "",
    "Required positional arguments:",
    "  <r_dir>  Path to the h0testr package R/ source directory; all .R files",
    "             there are sourced (the installed package is not used).",
    "",
    "Output: one PASS/FAIL line per assertion to stdout, then a count of",
    "  passes and failures. Log output, including messages from expected",
    "  errors, is written to a temporary file whose path is reported at the end.",
    "",
    "Exit codes: 0 all assertions passed; 1 one or more failed; 2 usage error;",
    "  3 a required test package is not installed.",
    "",
    "Examples:",
    "  Rscript test_contrast.R C:/path/to/h0test/h0testr/R",
    "  Rscript test_contrast.R ../../h0test/h0testr/R",
    "  Rscript test_contrast.R C:/path/to/h0testr/R > test_contrast.out 2>&1",
    sep="\n", file=stderr()
  )
  quit(status=2)
}

args <- commandArgs(trailingOnly=TRUE)
if(length(args) != 1) usage("wrong number of arguments")

r_dir <- args[1]
if(!dir.exists(r_dir)) usage(paste("r_dir not a directory:", r_dir))

for(pkg in c("limma", "lmtest", "DEqMS", "msqrob2", "proDA", "prolfqua")) {
  if(!requireNamespace(pkg, quietly=TRUE)) {
    cat("ERROR:", pkg, "not installed\n", file=stderr())
    quit(status=3)
  }
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

log_len <- function() {
  if(!file.exists(log_file)) return(0L)
  return(length(readLines(log_file, warn=FALSE)))
}

log_has <- function(n0, ...) {
  txt <- readLines(log_file, warn=FALSE)
  if(length(txt) <= n0) return(FALSE)
  txt <- paste(txt[-seq_len(n0)], collapse=" ")
  for(pat in c(...)) if(!grepl(pat, txt, fixed=TRUE)) return(FALSE)
  return(TRUE)
}

## an expected error: TRUE iff the call fails and the log picks up every pattern:

errs_with <- function(expr, ...) {
  n0 <- log_len()
  out <- try(suppressMessages(expr), silent=TRUE)
  if(!inherits(out, "try-error")) return(FALSE)
  if(length(list(...)) %in% 0) return(TRUE)
  return(log_has(n0, ...))
}

close_enough <- function(a, b, tol=1e-8) {
  if(length(a) != length(b)) return(FALSE)
  i <- !is.na(a) & !is.na(b)
  if(!any(i)) return(FALSE)
  if(any(is.na(a) != is.na(b))) return(FALSE)
  return(max(abs(a[i] - b[i])) < tol)
}

###############################################################################
## Twenty-four observations, fully crossed: grp has three levels, so it carries the
##   contrasts of interest, sex and batch two each, age continuous. Sixty genes with
##   two peptides each; the first ten are up in grp c only, so a contrast involving
##   level c should find them and one between a and b should not.

set.seed(101)

nobs <- 24
ngene <- 60

samps <- expand.grid(rep=1:2, batch=c("b1", "b2"), sex=c("F", "M"),
  grp=c("a", "b", "c"), stringsAsFactors=FALSE)
samps$obs <- paste0("o", sprintf("%02d", 1:nobs))
samps$age <- round(rnorm(nobs, 50, 10), 1)
samps <- samps[, c("obs", "grp", "sex", "batch", "age")]

gene_eff <- matrix(0, nrow=ngene, ncol=nobs)
gene_eff[1:10, samps$grp %in% "c"] <- 3

## peptide counts vary by gene, since DEqMS fits its variance against that count and
##   test_deqms() refuses a run in which every gene has the same number of them, which
##   would leave the engine this suite most wants to reach untested. Variances vary by
##   gene too, so that the moderating engines have a prior to estimate:

npep <- rep(c(2, 3, 2, 1), length.out=ngene)
gene <- rep(paste0("g", sprintf("%02d", 1:ngene)), times=npep)
nfeat <- length(gene)

sds <- rep(sqrt(1 / rgamma(ngene, shape=3, rate=3)), times=npep)
exprs <- matrix(rnorm(nfeat * nobs, 20, rep(sds, times=nobs)), nrow=nfeat)
exprs <- exprs + gene_eff[rep(1:ngene, times=npep), , drop=FALSE]
rownames(exprs) <- paste0(gene, ".p", unlist(lapply(npep, seq_len)))
colnames(exprs) <- samps$obs

feats <- data.frame(pep=rownames(exprs), gene=gene)

cfg0 <- new_config()
cfg0$obs_id_col <- cfg0$obs_col <- "obs"
cfg0$sample_id_col <- "obs"
cfg0$feat_id_col <- cfg0$feat_col <- "pep"
cfg0$gene_id_col <- "gene"
cfg0$frm <- ~grp + sex + batch + age
cfg0$test_term <- "grp"
cfg0$reference_levels <- c(grp="a", sex="F", batch="b1")
cfg0$estimability <- "test"
cfg0$df_resid_min <- 2
cfg0$save_state <- FALSE
cfg0$permute_var <- ""
cfg0$is_log_transformed <- TRUE
cfg0$impute_method <- "none"
cfg0$n_features_min <- 5
cfg0$test_prior_df <- 5
cfg0$log_file <- log_file

out <- initialize(list(expression=exprs, features=feats, samples=samps), cfg0,
  minimal=TRUE)
state0 <- out$state
config0 <- out$config

## a config testing a contrast: config$test_term must be "", the two being mutually
##   exclusive, which is asserted in its own section below:

cfg_con <- function(frm, contrast, method=NULL) {
  cfg <- config0
  cfg$frm <- frm
  cfg$test_term <- ""
  cfg$contrast <- contrast
  cfg$reference_levels <-
    config0$reference_levels[names(config0$reference_levels) %in% all.vars(frm)]
  if(!is.null(method)) cfg$test_method <- method
  return(cfg)
}

cfg_trm <- function(frm, test_term, method=NULL) {
  cfg <- config0
  cfg$frm <- frm
  cfg$test_term <- test_term
  cfg$contrast <- ""
  cfg$reference_levels <-
    config0$reference_levels[names(config0$reference_levels) %in% all.vars(frm)]
  if(!is.null(method)) cfg$test_method <- method
  return(cfg)
}

## gene 1 to 10 are the known effects, whether the ids are peptides or genes:

is_effect <- function(ids) {
  as.integer(sub("^g0*([0-9]+).*$", "\\1", as.character(ids))) <= 10
}

## the weighted sum of the least squares coefficients, computed here rather than by
##   any h0testr code, one value per row of exprs:

ind_effect <- function(mat, X, L) {
  out <- apply(mat, 1, function(y) {
    b <- stats::lm.fit(X, y)$coefficients
    b[is.na(b)] <- 0
    return(sum(L * b))
  })
  return(out)
}

## the exact F test of the full design against the constrained one, again computed
##   here; this is what an unmoderated engine should report:

ind_nested_p <- function(mat, X, X_red) {
  df1 <- qr(X)$rank - qr(X_red)$rank
  df2 <- ncol(mat) - qr(X)$rank
  out <- apply(mat, 1, function(y) {
    rss1 <- sum(stats::lm.fit(X, y)$residuals^2)
    rss0 <- sum(stats::lm.fit(X_red, y)$residuals^2)
    stat <- ((rss0 - rss1) / df1) / (rss1 / df2)
    return(stats::pf(stat, df1, df2, lower.tail=FALSE))
  })
  return(out)
}

###############################################################################
section("the parser: weights are what the expression says")

X0 <- stats::model.matrix(~grp + sex + batch + age, data=state0$samples)
cols0 <- colnames(X0)

wts <- function(txt, cols=cols0) {
  cfg <- config0
  cfg$test_term <- ""
  cfg$contrast <- txt
  return(f.contrast_vector(cfg, cols))
}

L1 <- wts("grpb - grpc")
report(is.numeric(L1) && identical(names(L1), cols0) &&
  close_enough(as.numeric(L1[c("grpb", "grpc")]), c(1, -1)) &&
  all(L1[setdiff(cols0, c("grpb", "grpc"))] %in% 0),
  "parser: 'grpb - grpc' weights those two coefficients 1 and -1 and no others")

report(close_enough(as.numeric(wts("(grpb + grpc)/2")[c("grpb", "grpc")]),
  c(0.5, 0.5)), "parser: '(grpb + grpc)/2' averages the two levels")

report(close_enough(as.numeric(wts("2*grpb - 0.5*grpc")[c("grpb", "grpc")]),
  c(2, -0.5)), "parser: constants scale a coefficient on either side of '*'")

report(close_enough(as.numeric(wts("-grpb")["grpb"]), -1),
  "parser: unary minus negates a weight")

report(close_enough(as.numeric(wts("grpb - (grpc - grpb)")["grpb"]), 2),
  "parser: parentheses group, and a repeated coefficient accumulates")

Xi <- stats::model.matrix(~grp * sex, data=state0$samples)
report(close_enough(as.numeric(wts("`grpb:sexM` - `grpc:sexM`",
  colnames(Xi))[c("grpb:sexM", "grpc:sexM")]), c(1, -1)),
  "parser: a backquoted interaction coefficient name is accepted")

report(close_enough(as.numeric(wts("`(Intercept)` + grpb")["(Intercept)"]), 1),
  "parser: the intercept can be weighted, under its backquoted name")

###############################################################################
section("the parser: errors name what was wrong")

report(errs_with(wts("grpz - grpb"), "is not a coefficient of the design matrix",
  "grpz", "coefficients available"),
  "parser error: an unknown coefficient name is named, with what was available")

report(errs_with(wts("log(grpb)"), "not supported in config$contrast"),
  "parser error: a function call is refused")

report(errs_with(wts("grpb:grpc"), "not supported in config$contrast"),
  "parser error: ':' is refused, an interaction coefficient being one name")

report(errs_with(wts("grpb * grpc"), "multiplies two", "must be a constant"),
  "parser error: multiplying two coefficient names is refused")

report(errs_with(wts("grpb / grpc"), "divisor of '/'", "must be", "constant"),
  "parser error: dividing by a coefficient name is refused")

report(errs_with(wts("grpb - grpb"), "every coefficient weight", "is zero"),
  "parser error: weights that cancel to zero leave no hypothesis")

report(errs_with(wts("grpb -"), "does not parse as an R expression"),
  "parser error: an unparseable expression says so")

###############################################################################
section("the constrained design")

cfg_bc <- cfg_con(~grp + sex + batch + age, "grpc - grpb")
d_bc <- f.design_test_cols(state0, cfg_bc)

report(!is.null(d_bc$contrast) && d_bc$df_intend %in% 1,
  "design: a contrast is one degree of freedom")

report(ncol(d_bc$X_red) == ncol(d_bc$X) - 1 &&
  f.design_rank(d_bc$X_red) == f.design_rank(d_bc$X) - 1,
  "design: the constrained design is one column and one rank below the full one")

report(f.design_rank(cbind(d_bc$X, d_bc$X_red)) == f.design_rank(d_bc$X),
  "design: the constrained design is nested in the full one")

report(identical(colnames(d_bc$X)[d_bc$cols_test], c("grpb", "grpc")),
  "design: cols_test are the coefficients carrying a non-zero weight")

d_trm <- f.design_test_cols(state0, cfg_trm(~grp + sex + batch + age, "grp"))
report(is.null(d_trm$contrast) && d_trm$df_intend %in% 2 &&
  identical(d_trm$X_red, d_trm$X[, -d_trm$cols_test, drop=FALSE]),
  "design: a term test is unchanged, and its X_red is a subset of the columns")

###############################################################################
section("estimability of the contrast itself")

## dup is grp under another set of level names, so its columns duplicate grp's and
##   the design is rank deficient. A contrast along the aliased direction is not
##   estimable; the sum of the two duplicated coefficients is:

samps2 <- state0$samples
samps2$dup <- c(a="A", b="B", c="C")[as.character(samps2$grp)]
st2 <- state0
st2$samples <- samps2

cfg_al <- cfg_con(~grp + dup, "grpb - dupB")
cfg_al$reference_levels <- c(grp="a", dup="A")
cfg_al <- initialize(st2, cfg_al, minimal=TRUE)$config

report(errs_with(f.design_test_cols(st2, cfg_al), "is not", "estimable",
  "row space of the design matrix"),
  "estimability: a contrast along an aliased direction is a hard error")

cfg_al2 <- cfg_al
cfg_al2$contrast <- "grpb + dupB"
d_al2 <- try(suppressMessages(f.design_test_cols(st2, cfg_al2)), silent=TRUE)
report(!inherits(d_al2, "try-error") && d_al2$df_intend %in% 1,
  "estimability: a contrast that is estimable despite the aliasing runs")

###############################################################################
section("config$contrast and config$test_term are mutually exclusive")

cfg_both <- cfg_con(~grp + sex + batch + age, "grpc - grpb")
cfg_both$test_term <- "grp"

report(errs_with(check_config(cfg_both), "both name", "one run tests one hypothesis",
  "set config$test_term"),
  "config: setting both test_term and contrast is refused, saying how to fix it")

report(!inherits(try(suppressMessages(check_config(cfg_bc)), silent=TRUE),
  "try-error"), "config: a contrast with test_term \"\" passes check_config()")

report(!inherits(try(suppressMessages(report_config(cfg_bc)), silent=TRUE),
  "try-error"),
  "config: report_config() does not look for a test_term in a contrast run")

st_init <- try(suppressMessages(initialize(
  list(expression=exprs, features=feats, samples=samps), cfg_bc, minimal=TRUE)),
  silent=TRUE)
report(!inherits(st_init, "try-error"),
  "config: initialize() accepts a config testing a contrast")

###############################################################################
section("the marginality warning")

n0 <- log_len()
cfg_marg <- cfg_con(~grp * sex, "grpc - grpb")
d_marg <- suppressMessages(f.design_test_cols(state0, cfg_marg))

report(log_has(n0, "WARNING", "weights coefficients of the term", "grp:sex",
  "sex = F"),
  "marginality: a contrast inside an interaction warns, naming the held level")

report(d_marg$df_intend %in% 1,
  "marginality: the contrast is still a one degree of freedom test")

n0 <- log_len()
d_nomarg <- suppressMessages(f.design_test_cols(state0, cfg_bc))
report(!log_has(n0, "weights coefficients of the term"),
  "marginality: no warning when no higher-order term contains the tested term")

###############################################################################
section("per-feature estimability screening")

## one peptide seen only in grp a, so neither grpb nor grpc, and therefore no
##   contrast between them, is estimable for it:

st3 <- state0
st3$expression[1, !(st3$samples$grp %in% "a")] <- NA

cfg_est <- cfg_bc
cfg_est$estimability <- "test"
st3b <- suppressMessages(filter_features_by_estimability(st3, cfg_est))

report(!(rownames(state0$expression)[1] %in% rownames(st3b$expression)),
  "estimability: a feature whose data cannot separate the contrast is dropped")

report(nrow(st3b$expression) == nrow(state0$expression) - 1,
  "estimability: the features that can separate it are all kept")

report(all(st3b$features[[config0$df_test_col]] %in% 1),
  "estimability: df_test is 1 per feature for a contrast")

###############################################################################
section("the engines: all seven run the contrast")

methods <- c("lm", "trend", "voom", "deqms", "msqrob", "proda", "prolfqua")
res <- list()

for(m in methods) {
  cfg <- cfg_con(~grp + sex + batch + age, "grpc - grpb", method=m)
  out1 <- try(suppressMessages(test(state0, cfg)), silent=TRUE)
  res[[m]] <- out1
  ok <- !inherits(out1, "try-error") && is.data.frame(out1$standard) &&
    identical(names(out1$standard), c("feature", "expr", "logfc", "stat", "lod",
      "pval", "adj_pval")) && any(!is.na(out1$standard$pval))
  report(ok, paste0("engine ", m, ": runs the contrast and returns the standard table"))
}

## the known effect is in grp c, so this contrast should rank those genes first:

for(m in methods) {
  out1 <- res[[m]]
  if(inherits(out1, "try-error")) {
    report(FALSE, paste0("engine ", m, ": ranks the known effects first"))
    next
  }
  tbl <- out1$standard
  tbl <- tbl[order(tbl$pval), ]
  top <- utils::head(tbl$feature, if(m %in% c("deqms", "msqrob")) 10 else 20)
  report(sum(is_effect(top)) >= if(m %in% c("deqms", "msqrob")) 8 else 16,
    paste0("engine ", m, ": ranks the known effects first"))
}

###############################################################################
section("the effect size is the weighted sum of the coefficients")

L_bc <- d_bc$contrast
eff_ind <- ind_effect(state0$expression, d_bc$X, L_bc)

for(m in c("lm", "trend", "prolfqua")) {
  tbl <- res[[m]]$standard
  got <- tbl$logfc[match(names(eff_ind), tbl$feature)]
  report(close_enough(got, as.numeric(eff_ind), tol=1e-6),
    paste0("engine ", m, ": logfc is the weighted sum of the least squares",
      " coefficients"))
}

## voom weights, proDA models dropout and msqrob2 fits robustly, so their
##   coefficients are their own; the effect size should still agree closely in
##   direction and size with the least squares one:

## proDA is the loosest of the three, its coefficients coming from a dropout model
##   fitted with its own regularization rather than from a reweighting of the same
##   least squares problem, so it gets a looser bound than voom rather than none:

for(m in c("voom", "proda")) {
  tbl <- res[[m]]$standard
  got <- tbl$logfc[match(names(eff_ind), tbl$feature)]
  i <- !is.na(got)
  lim <- if(m %in% "proda") 0.8 else 0.99
  report(sum(i) > 100 && stats::cor(got[i], as.numeric(eff_ind)[i]) > lim,
    paste0("engine ", m, ": logfc tracks the least squares weighted sum"))
}

report(all(!is.na(res[["msqrob"]]$standard$logfc)) &&
  stats::cor(res[["msqrob"]]$standard$logfc,
    tapply(eff_ind, feats$gene, mean)[res[["msqrob"]]$standard$feature]) > 0.99,
  "engine msqrob: logfc is the contrast, from msqrob2::hypothesisTest()")

## the simulated difference between grp c and grp b is 3, and every gene carrying it
##   is estimated further from zero than every gene that does not, without a threshold
##   having to be picked:

i_eff <- is_effect(names(eff_ind))
report(all(abs(eff_ind[i_eff] - 3) < 1.5) &&
  min(eff_ind[i_eff]) > max(abs(eff_ind[!i_eff])),
  "effect size: the contrast recovers the simulated size of 3, separating the two")

###############################################################################
section("the p-values are the test of the contrast")

## test_prolfqua() with moderation off is an exact nested F test, so it can be
##   checked against one computed here from the same two designs:

cfg_ex <- cfg_con(~grp + sex + batch + age, "grpc - grpb", method="prolfqua")
cfg_ex$test_moderate <- FALSE
cfg_ex$test_trend <- FALSE
res_ex <- suppressMessages(test(state0, cfg_ex))
p_ind <- ind_nested_p(state0$expression, d_bc$X, d_bc$X_red)
got <- res_ex$standard$pval[match(names(p_ind), res_ex$standard$feature)]

report(close_enough(got, as.numeric(p_ind), tol=1e-8),
  "prolfqua unmoderated: p-values match an independent nested F test exactly")

## the same statistic from limma, whose moderated t on the contrast squares to an F
##   on one denominator degree of freedom more than the unmoderated test has:

tbl_tr <- res[["trend"]]$standard
i_tr <- is_effect(tbl_tr$feature)
report(all(tbl_tr$stat[i_tr] > 0) && sum(tbl_tr$pval[i_tr] < 1e-3) >= 15 &&
  all(tbl_tr$pval[i_tr] < 0.05),
  "trend: the contrast is signed, and positive for genes up in grp c")

## a contrast between the two levels with no simulated difference finds nothing,
##   which the joint test of grp cannot distinguish from the effect above:

cfg_ab <- cfg_con(~grp + sex + batch + age, "grpb", method="trend")
res_ab <- suppressMessages(test(state0, cfg_ab))
report(sum(res_ab$standard$adj_pval < 0.05, na.rm=TRUE) <= 2,
  "trend: a contrast on the level with no simulated effect finds almost nothing")

###############################################################################
section("a contrast is one coefficient, whatever the engine could do without it")

## this section tested what the contrast bought: deqms was capped at one coefficient,
##   and a contrast was the only way it could reach a factor with more than two levels.
##   No engine is capped any more, so what is left to check is that a contrast still
##   arrives at the limma-family engines as a single coefficient, which is what makes it
##   one degree of freedom there and a different hypothesis from the term it lives in:

report(!exists("f.test_max_cols") && !exists("f.design_test_cols_max"),
  "caps: the cap helpers are gone, deqms having been the last capped engine")

res_dj <- try(suppressMessages(test_deqms(state0,
  cfg_trm(~grp + sex + batch + age, "grp", method="deqms"))), silent=TRUE)
report(!inherits(res_dj, "try-error") &&
  all(res_dj$hits$sca.df.num %in% 2) && all(is.na(res_dj$hits$sca.t)),
  "caps: test_deqms() runs the joint test of a three-level factor, on 2 df")

report(!inherits(res[["deqms"]], "try-error") &&
  any(!is.na(res[["deqms"]]$standard$pval)),
  "caps: the same three-level factor is also reachable by contrast under deqms")

d_con <- suppressMessages(f.design_test_cols(state0, cfg_bc))
report(!is.null(d_con$contrast) && length(d_con$cols_test) %in% 2,
  "caps: f.design_test_cols() returns the contrast over both of its coefficients")

hits_dc <- suppressMessages(test_deqms(state0, cfg_bc))$hits
report(all(hits_dc$sca.df.num %in% 1) && !any(is.na(hits_dc$sca.t)),
  "caps: that contrast still reaches deqms as one coefficient, on 1 df")

## msqrob takes the moderated t path, not the joint Wald path of f.msqrob_wald():

hits_msq <- suppressMessages(test_msqrob(state0,
  cfg_con(~grp + sex + batch + age, "grpc - grpb", method="msqrob")))$hits
report(all(c("logFC", "se", "df", "t", "pval", "adjPval") %in% names(hits_msq)) &&
  !("f_statistic" %in% names(hits_msq)),
  "caps: test_msqrob() answers a contrast with its own moderated t, not a joint F")

## and the contrast is a different hypothesis from the term it lives in:

cfg_j <- cfg_trm(~grp + sex + batch + age, "grp", method="trend")
res_j <- suppressMessages(test(state0, cfg_j))
p_j <- res_j$standard$pval[match(tbl_tr$feature, res_j$standard$feature)]
report(!close_enough(p_j, tbl_tr$pval, tol=1e-6),
  "the contrast and the joint test of the same factor are different hypotheses")

###############################################################################

cat("\n## log file:", log_file, "\n")
cat("## elapsed:", round(as.numeric(difftime(Sys.time(), t0, units="secs"))),
  "s\n")
cat("\n## passes:", n_pass, "; failures:", n_fail, "; elapsed:",
  round(as.numeric(difftime(Sys.time(), t0, units="secs"))), "s\n")

quit(save="no", status=if(n_fail > 0) 1 else 0)
