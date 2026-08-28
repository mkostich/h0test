## Tests for the trended variance prior reaching every engine that can honor it.
##   h0testr::test() gained a trend argument, config$test_trend answers when it is not
##   given, and the resolved value is passed to test_method "deqms", where it sets
##   limma::eBayes(trend=) beneath DEqMS's own count-based prior, and to "prolfqua", where
##   it fits the prior against mean feature intensity. Before this, config$test_trend was
##   read by "prolfqua" alone and the trended DEqMS fit could not be reached through
##   test() at all. A request that the resolved method cannot honor is a warning rather
##   than a refusal, since h0testr::tune() sets one configuration and varies the method.
##   The three methods are therefore in three positions, and this file checks each: only
##   "prolfqua" trends in a way that moves what test() reports, which is what
##   f.trend_methods() now names; "deqms" honors the argument in limma::eBayes() but
##   refits its prior from the spectra counts afterward, so its reported p-value does not
##   move and f.note_trend() says so rather than letting the run be silent; and the rest
##   cannot trend at all.
##   Covers f.is_trend(), f.trend_methods() and f.note_trend() directly, then confirms the
##   pass-through with real fits.

usage <- function(msg=NULL) {
  if(!is.null(msg)) cat("ERROR:", msg, "\n\n", file=stderr())
  cat(
    "Test that the trended variance prior reaches the engines that can fit one: that",
    "h0testr::test(trend=) overrides config$test_trend and config$test_trend answers when",
    "the argument is absent, that both reach test_method 'deqms' and 'prolfqua', that",
    "'trend' always trends and only an explicit FALSE draws a remark there, that a",
    "trended 'deqms' run is warned about because DEqMS refits the prior past it, and",
    "that a method which cannot trend warns rather than refusing.",
    "",
    "Usage: Rscript test_trend_flag.R <r_dir>",
    "",
    "Required positional arguments:",
    "  <r_dir>  Path to the h0testr package R/ source directory; all .R files",
    "             there are sourced (the installed package is not used).",
    "",
    "Optional named arguments: none.",
    "",
    "Needs the DEqMS and prolfqua packages, whose priors are the ones under test, and",
    "  limma beneath both.",
    "",
    "Output: one PASS/FAIL line per assertion to stdout, then a count of passes and",
    "  failures and the elapsed time. Engine messages go to a temporary log file, whose",
    "  path is reported at the end.",
    "",
    "Exit codes: 0 all assertions passed; 1 one or more failed; 2 usage error.",
    "",
    "Examples:",
    "  Rscript test_trend_flag.R C:/path/to/h0test/h0testr/R",
    "  Rscript test_trend_flag.R ../../h0test/h0testr/R",
    "  Rscript test_trend_flag.R C:/path/to/h0testr/R > trend.out 2>&1",
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

threw <- function(expr) inherits(try(expr, silent=TRUE), "try-error")

mark <- function() if(file.exists(log_file)) length(readLines(log_file)) else 0

log_since <- function(since=0) {
  if(!file.exists(log_file)) return(character(0))
  txt <- readLines(log_file)
  if(since >= length(txt)) return(character(0))
  return(txt[(since + 1):length(txt)])
}

logged <- function(pat, since=0) any(grepl(pat, log_since(since), fixed=TRUE))
log_count <- function(pat, since=0) sum(grepl(pat, log_since(since), fixed=TRUE))

## save_state FALSE and dir_out a temporary directory: new_config() ships save_state TRUE and
##   dir_out ".", and the fixtures here reach test(), which writes its results table into the
##   working directory named after length(config$run_order) + 3 -- an 8.results.*.tsv pair easy
##   to mistake for the output of a real run. Same guard as 1/test_check_config.R:

base_cfg <- function() {
  cfg <- new_config()
  cfg$log_file <- log_file
  cfg$save_state <- FALSE
  cfg$dir_out <- tempdir()
  return(cfg)
}

###############################################################################
section("f.is_trend(), which settles the value one place for all of them")

cfg <- base_cfg()
report(f.is_trend(NULL, cfg) %in% FALSE,
  "the h0testr::new_config() default is a flat prior")

cfg$test_trend <- NULL
report(f.is_trend(NULL, cfg) %in% FALSE,
  "and so is a config that never mentions it, absence not being an error here")

cfg <- base_cfg()
cfg$test_trend <- TRUE
report(f.is_trend(NULL, cfg) %in% TRUE, "config$test_trend answers when nothing is passed")
report(f.is_trend("", cfg) %in% TRUE, "an empty string counts as nothing passed")

## unlike f.is_log_transformed(), which errors when the two disagree: the scale of the
##   data is a fact and only one of them can be right, whereas this is a preference about
##   how to fit a prior, so a caller passing one is overriding the configuration on purpose:

report(f.is_trend(FALSE, cfg) %in% FALSE,
  "an argument that disagrees with config$test_trend simply wins")

cfg <- base_cfg()
report(f.is_trend(TRUE, cfg) %in% TRUE, "and wins in the other direction too")

m0 <- mark()
report(threw(f.is_trend(NA, cfg)), "an NA is refused, absence being how it is left unset")
report(logged("trend is not TRUE or FALSE", m0), "and the message names the offense")
report(threw(f.is_trend("yes", cfg)), "so is a value that is not logical")
report(threw(f.is_trend(c(TRUE, FALSE), cfg)), "so is a vector of them")

cfg <- base_cfg()
cfg$test_trend <- "TRUE"
report(threw(f.is_trend(NULL, cfg)),
  "a config$test_trend that is not TRUE or FALSE is refused as well")

###############################################################################
section("f.trend_methods() and f.note_trend(), what is said to the rest")

## the list names the methods where trending moves what test() reports, which is prolfqua
##   alone. deqms was here until it was observed that its reported p-value comes from a
##   prior refitted past the flag, so a trended deqms run was silent while being no
##   different from a flat one; it is now its own case below:

report(setequal(f.trend_methods(), "prolfqua"),
  "one method takes the flag and reports a different answer for it")
report(all(f.trend_methods() %in% test_methods()),
  "and it is a name h0testr::test_methods() offers")
report(!("deqms" %in% f.trend_methods()),
  "deqms is not among them, its own prior being refitted past the flag")
report(!("trend" %in% f.trend_methods()),
  "test_method 'trend' is not among them, taking no argument because it always trends")

cfg <- base_cfg()

for(nm in f.trend_methods()) {
  report(f.note_trend(nm, TRUE, TRUE, cfg) %in% FALSE,
    paste("nothing is said to method", nm, "which honors the request"))
}

report(f.note_trend("trend", TRUE, TRUE, cfg) %in% FALSE,
  "nor to method trend, which is already the fit that was asked for")

m0 <- mark()
report(f.note_trend("trend", FALSE, TRUE, cfg) %in% TRUE,
  "but an explicit FALSE there is a request it cannot honor")
report(logged("regardless", m0), "and says the prior is trended anyway")
report(logged("'deqms'", m0), "pointing at the method whose limma prior can be flat")

report(f.note_trend("trend", FALSE, FALSE, cfg) %in% FALSE,
  "while the default FALSE of a config is silent, or every run would carry the remark")

cannot <- c("lm", "voom", "msqrob", "msqrob_agg", "proda", "prolfqua_lmer")

for(nm in cannot) {
  m0 <- mark()
  report(f.note_trend(nm, TRUE, FALSE, cfg) %in% TRUE,
    paste("method", nm, "cannot trend, so a TRUE is remarked on"))
  report(logged("config$test_trend", m0), paste("naming the configuration key for", nm))
  report(!logged("Stopping", m0), paste("as a warning and not a refusal for", nm))
}

m0 <- mark()
report(f.note_trend("voom", TRUE, TRUE, cfg) %in% TRUE, "an explicit TRUE too")
report(logged("the trend argument", m0),
  "with the argument named instead, that being what there is to change")
report(logged("the methods that do are prolfqua", m0),
  "and the methods that would honor it listed")

for(nm in c(cannot, "trend", "deqms", f.trend_methods())) {
  report(f.note_trend(nm, FALSE, FALSE, cfg) %in% FALSE,
    paste("a flat prior is unremarkable on method", nm))
}

## deqms is neither of the two cases above: the argument is honored, in limma::eBayes(),
##   but DEqMS refits the prior from the spectra counts and the reported p-value does not
##   move. So a TRUE draws a remark, as it does for a method that cannot trend, but not the
##   same remark: the generic one says the method does not fit its prior against mean
##   feature intensity, which of deqms would be false:

m0 <- mark()
report(f.note_trend("deqms", TRUE, FALSE, cfg) %in% TRUE,
  "a trended deqms run is remarked on rather than being silent")
report(logged("config$test_trend", m0), "naming the configuration key it came from")
report(logged("test_deqms() passes it to limma::eBayes()", m0),
  "and saying the argument is honored, which it is")
report(logged("DEqMS refits the prior from the spectra counts", m0),
  "and why that does not change the answer")
report(logged("P.Value, t, B, s2.prior and s2.post", m0),
  "naming what does move, which is limma's own columns of the returned table")
report(logged("test_method 'trend'", m0),
  "and pointing at the method whose trended prior does move the reported answer")
report(!logged("does not fit", m0),
  "without the generic clause, which would be false of deqms")
report(!logged("Stopping", m0), "as a warning and not a refusal")

m0 <- mark()
report(f.note_trend("deqms", TRUE, TRUE, cfg) %in% TRUE,
  "an explicit TRUE draws it too")
report(logged("the trend argument", m0),
  "with the argument named instead, that being what there is to change")

###############################################################################
section("test_method deqms, where the flag could not be reached before")

## DEqMS fits its variance prior against the number of features behind each gene, so the
##   fixture needs that number to vary; limma's prior underneath it is the one this flag
##   sets. Raw scale, then log2 after initialize(), so that the zeros became NA first, as
##   normalize() would do in a full workflow:

set.seed(202)
n_obs <- 12
peps_per <- rep(c(2, 3, 4, 5), length.out=40)
gene <- rep(paste0("g", sprintf("%02d", seq_along(peps_per))), peps_per)
pep <- paste0(gene, "_p", unlist(lapply(peps_per, seq_len)))

mu <- rnorm(length(pep), mean=12, sd=1.5)
lg <- matrix(rnorm(length(pep) * n_obs, mean=mu, sd=0.6), nrow=length(pep))
grp <- c(rep("ctl", n_obs / 2), rep("trt", n_obs / 2))
lg[gene %in% paste0("g", sprintf("%02d", 1:10)), grp %in% "trt"] <- 1 +
  lg[gene %in% paste0("g", sprintf("%02d", 1:10)), grp %in% "trt"]

rownames(lg) <- pep
colnames(lg) <- paste0("s", sprintf("%02d", 1:n_obs))

state_a <- list(
  expression=2 ^ lg,
  features=data.frame(pep=pep, gene=gene, stringsAsFactors=FALSE),
  samples=data.frame(obs=colnames(lg), grp=grp, stringsAsFactors=FALSE)
)

cfg_a <- base_cfg()
cfg_a$feat_id_col <- "pep"
cfg_a$gene_id_col <- "gene"
cfg_a$obs_id_col <- cfg_a$sample_id_col <- "obs"
cfg_a$frm <- ~grp
cfg_a$test_term <- "grp"
cfg_a$reference_levels <- c(grp="ctl")
cfg_a$test_method <- "deqms"

out_a <- try(suppressMessages(initialize(state_a, cfg_a, minimal=TRUE)), silent=TRUE)
report(!inherits(out_a, "try-error"), "the peptide level fixture initializes")

out_a$state$expression <- log2(out_a$state$expression + 1)
out_a$config$is_log_transformed <- TRUE

report(length(unique(f.gene_counts(out_a$state, out_a$config))) > 1,
  "and the number of peptides per gene varies, which is what DEqMS needs")

m0 <- mark()
res_flat <- try(suppressMessages(test(out_a$state, out_a$config)), silent=TRUE)
report(!inherits(res_flat, "try-error"), "the untrended fit runs")
report(logged("(trend=FALSE)", m0),
  "and the log records which limma prior was fitted, the hit table not carrying it")
report(!logged("WARNING: test:", m0),
  "with nothing remarked on, a flat prior being what deqms reports from either way")

cfg_t <- out_a$config
cfg_t$test_trend <- TRUE
m0 <- mark()
res_cfg <- try(suppressMessages(test(out_a$state, cfg_t)), silent=TRUE)
report(!inherits(res_cfg, "try-error"), "so does the trended one")
report(logged("(trend=TRUE)", m0),
  "and config$test_trend reaches deqms, which it did not before")
report(logged("WARNING: test:", m0),
  "but the run is remarked on now, the reported answer not moving for it")
report(logged("DEqMS refits the prior from the spectra counts", m0),
  "with the reason, rather than the run being silent as it was")

m0 <- mark()
res_arg <- try(suppressMessages(test(out_a$state, out_a$config, trend=TRUE)),
  silent=TRUE)
report(!inherits(res_arg, "try-error") && logged("(trend=TRUE)", m0),
  "the argument reaches it from a config that says otherwise")
report(logged("the trend argument", m0),
  "and the remark names the argument, that being what there is to change")

m0 <- mark()
res_off <- try(suppressMessages(test(out_a$state, cfg_t, trend=FALSE)), silent=TRUE)
report(!inherits(res_off, "try-error") && logged("(trend=FALSE)", m0),
  "and turns it back off against a config that asks for it")
report(!logged("WARNING: test:", m0),
  "with no remark, an honored FALSE being nothing to say anything about")

report(identical(res_arg$standard$pval, res_cfg$standard$pval),
  "the two routes to the trended fit give the same p-values")
report(identical(res_off$standard$pval, res_flat$standard$pval),
  "as do the two routes to the flat one")

## and what the flag does and does not move here, which is the whole of why deqms is
##   documented as taking it without it changing the answer: DEqMS::spectraCounteBayes()
##   fits its prior from fit$sigma, fit$df.residual and fit$count, none of which
##   limma::eBayes() alters, and f.deqms_moderated_f() reports from that prior:

report(identical(res_flat$standard$pval, res_cfg$standard$pval),
  "the reported statistic is the same either way, coming from DEqMS's own prior")
report(!isTRUE(all.equal(res_flat$original$P.Value, res_cfg$original$P.Value)),
  "while limma's own p-value, carried alongside in the original table, does move")
report(!isTRUE(all.equal(res_flat$original$B, res_cfg$original$B)),
  "as does limma's log odds, the two priors being different quantities")
report(logged("which is unaffected", 0),
  "and the log says which of the two priors the flag reached")

###############################################################################
section("test_method prolfqua, which read config$test_trend already")

## a feature level fixture: the least squares path fits one model per row, and the
##   trended prior needs more rows than the spline has degrees of freedom:

set.seed(303)
exprs <- sim1(n_obs=12, n_feats=40)$mat
exprs <- log2(exprs + 1)

state_b <- list(
  expression=exprs,
  features=data.frame(feature_id=rownames(exprs), stringsAsFactors=FALSE),
  samples=data.frame(observation_id=colnames(exprs),
    grp=c(rep("ctl", 6), rep("trt", 6)), stringsAsFactors=FALSE)
)

cfg_b <- base_cfg()
cfg_b$feat_id_col <- cfg_b$gene_id_col <- "feature_id"
cfg_b$obs_id_col <- cfg_b$sample_id_col <- "observation_id"
cfg_b$frm <- ~grp
cfg_b$test_term <- "grp"
cfg_b$reference_levels <- c(grp="ctl")
cfg_b$test_method <- "prolfqua"
cfg_b$is_log_transformed <- TRUE

out_b <- try(suppressMessages(initialize(state_b, cfg_b, minimal=TRUE)), silent=TRUE)
report(!inherits(out_b, "try-error"), "the feature level fixture initializes")

m0 <- mark()
pro_flat <- try(suppressMessages(test(out_b$state, out_b$config)), silent=TRUE)
report(!inherits(pro_flat, "try-error") && !any(pro_flat$original$trend),
  "the flat prior is what a default config asks for")
report(logged("against a flat prior", m0), "and the log says so")

m0 <- mark()
pro_arg <- try(suppressMessages(test(out_b$state, out_b$config, trend=TRUE)),
  silent=TRUE)
report(!inherits(pro_arg, "try-error") && all(pro_arg$original$trend),
  "the argument reaches the prior this method was already able to fit")
report(logged("against a trended prior", m0), "and the log says that instead")
report(!logged("WARNING: test:", m0), "with nothing remarked on here either")

cfg_b2 <- out_b$config
cfg_b2$test_trend <- TRUE
pro_cfg <- try(suppressMessages(test(out_b$state, cfg_b2)), silent=TRUE)
report(!inherits(pro_cfg, "try-error") && all(pro_cfg$original$trend),
  "config$test_trend still reaches it, the argument only overriding")
report(identical(pro_arg$standard$pval, pro_cfg$standard$pval),
  "by the same route, so the two agree")

pro_off <- try(suppressMessages(test(out_b$state, cfg_b2, trend=FALSE)), silent=TRUE)
report(!inherits(pro_off, "try-error") && !any(pro_off$original$trend),
  "and an argument of FALSE overrides a config that asks for the trend")
report(!isTRUE(all.equal(pro_flat$standard$pval, pro_cfg$standard$pval)),
  "the two priors giving different p-values")

###############################################################################
section("the methods that cannot trend, and the one that always does")

## a request that cannot be honored costs a line in the log and nothing else: tune() sets
##   one config and varies the method, so refusing would kill a sweep on the first method
##   that does not trend:

cfg_lm <- out_b$config
cfg_lm$test_trend <- TRUE
m0 <- mark()
res_lm <- try(suppressMessages(test(out_b$state, cfg_lm, method="lm")), silent=TRUE)
report(!inherits(res_lm, "try-error"), "test_method lm runs with config$test_trend TRUE")
report(logged("does not fit", m0), "with a warning that the run is unaffected")
report(logged("config$test_trend", m0), "naming the key to change")

m0 <- mark()
res_tr <- try(suppressMessages(test(out_b$state, cfg_lm, method="trend")), silent=TRUE)
report(!inherits(res_tr, "try-error"), "test_method trend runs with it TRUE as well")
report(!logged("does not fit", m0),
  "and says nothing, that method already fitting the prior that was asked for")

m0 <- mark()
res_tr2 <- try(suppressMessages(test(out_b$state, out_b$config, method="trend",
  trend=FALSE)), silent=TRUE)
report(!inherits(res_tr2, "try-error"), "an explicit FALSE there still runs")
report(logged("regardless", m0), "but is remarked on, the method having no flat prior")
report(identical(res_tr$standard$pval, res_tr2$standard$pval),
  "and the fit is the same one either way")

report(log_count("test: method: trend", 0) %in% 2,
  "the resolved value is reported once per run in the method line")

###############################################################################
cat("\n## passes:", n_pass, "; failures:", n_fail, "; elapsed:",
  round(as.numeric(difftime(Sys.time(), t0, units="secs"))), "s\n")
cat("## log file:", log_file, "\n")

quit(save="no", status=if(n_fail > 0) 1 else 0)
