#' Get vector of \code{test(method=)} method options
#' @description
#'   Get a vector with acceptable values of \code{method} parameter for \code{h0testr::test()}.
#' @return
#'   Character vector with names of acceptable values for \code{h0testr::normalize(method=)}.
#' @examples
#' test_methods <- h0testr::test_methods()
#' cat("Available test methods:\n")
#' for(method in test_methods) {
#'   cat("method:", method, "\n")
#' }

test_methods <- function() {
  return(
    c("lm", "trend", "deqms", "msqrob", "msqrob_agg", "proda", "prolfqua",
      "prolfqua_lmer", "voom")
  )
}

## helper for test(): the settings only some engines read, so that one set on a run whose
##   method does not consult it can be said to have had no effect. Reported rather than
##   refused, and only when the value differs from the one new_config() ships: tune() sets
##   one config and varies the method, so most combinations of a sweep carry a setting
##   meant for one of the others and flagging those would be noise. Said from here rather
##   than from check_config(), which every step of the workflow calls and which would
##   repeat it once per step; this is the same place, and once per run, as the notes
##   test_prolfqua() and test_msqrob() make about the settings they do read. Keyed on the
##   method test() resolved rather than on config$test_method, test(method=) overriding
##   it. Returns the settings noted, so that a caller can check without reading the log:

f.note_ignored_settings <- function(method, config) {

  ## config$test_moderate is read by the least squares prolfqua path and by nothing else,
  ##   so every other method silently ignored a FALSE, while the same mistake with
  ##   config$test_ridge drew a note. Complementary to the note f.prolfqua_mixed() makes
  ##   rather than a duplicate of it: this one speaks only when the value differs from the
  ##   one new_config() ships, so only for FALSE, and that one speaks only when it is TRUE,
  ##   so prolfqua_lmer draws one message or the other and never both:

  knobs <- list(
    test_random_obs=list(default=TRUE, methods=c("prolfqua_lmer", "msqrob_agg")),
    test_ridge=list(default=FALSE, methods="msqrob_agg"),
    test_moderate=list(default=TRUE, methods="prolfqua")
  )

  ignored <- character(0)

  for(nom in names(knobs)) {
    if(!(nom %in% names(config)) || !(length(config[[nom]]) %in% 1)) next
    if(is.na(config[[nom]]) || config[[nom]] %in% knobs[[nom]]$default) next
    if(method %in% knobs[[nom]]$methods) next
    ignored <- c(ignored, paste0("config$", nom, "=", config[[nom]],
      " (read by test_method ", paste(knobs[[nom]]$methods, collapse=", "), ")"))
  }

  if(length(ignored)) {
    f.msg("NOTE: test: method", method, "does not consult",
      paste(ignored, collapse="; "), "\n",
      " so", if(length(ignored) > 1) "those settings have" else "that setting has",
      "no effect on this run", config=config)
  }

  return(ignored)
}

## The test methods where fitting the variance prior against mean feature intensity, which
##   is what limma::eBayes(trend=TRUE) does and what config$test_trend, or test()'s trend
##   argument, asks for, changes the answer test() reports. Only prolfqua: its moderation
##   is the limma one, so the trended prior is limma::eBayes(trend=TRUE)'s and it moves the
##   reported p-value. test_method "trend" is not here although it always trends: these are
##   the methods a TRUE can be handed to, and that one takes no such argument because
##   trending is its definition. deqms is not here either and used to be, which is what
##   made a trended deqms run silent: it does hand the argument to limma::eBayes(), and
##   that moves the limma columns of the table it returns, but DEqMS then refits the prior
##   from the spectra counts and every sca. column, which is what test() reports, comes out
##   unchanged; so a caller who asked for a trended answer did not get one and nothing said
##   so. See f.note_trend(), which now says it. The rest cannot trend at all: lm and msqrob
##   do not moderate against a covariate, proda fits its own prior, voom puts the
##   mean-variance relationship into precision weights so a trended prior on top of it
##   would count the same thing twice, and the mixed paths have a Satterthwaite denominator
##   rather than one residual variance to shrink:

f.trend_methods <- function() {
  return("prolfqua")
}

## helper for test(): a request to trend that the resolved method cannot honor. A warning
##   rather than an error, per the same reasoning as f.note_ignored_settings(): tune() sets
##   one config and varies the method, so a sweep would otherwise die on the first method
##   that does not trend. Keyed on the resolved value rather than on config$test_trend, so
##   that an explicit test(trend=TRUE) is caught too, and naming the source, since which of
##   the two to change differs. test_method "trend" gets the opposite message: it always
##   trends, so TRUE is already what it does and only an explicit FALSE is a request it
##   cannot honor. Returns TRUE when something was said, so that a test need not read the
##   log:

f.note_trend <- function(method, trend, given, config) {

  src <- if(given) "the trend argument" else "config$test_trend"

  if(method %in% "trend") {
    if(given && !trend) {
      f.msg("WARNING: test: trend is FALSE but test_method 'trend' is",
        "limma::eBayes(trend=TRUE), so the prior is fitted against mean gene",
        "intensity regardless;", "\n",
        " for a flat prior use test_method 'deqms', whose limma prior this argument",
        "does set, or 'lm', which does not moderate at all", config=config)
      return(TRUE)
    }
    return(FALSE)
  }

  ## deqms is the third case, and the reason it is no longer in f.trend_methods(). The
  ##   argument is honored: test_deqms() hands it to limma::eBayes(), which moves P.Value,
  ##   t, B, s2.prior and s2.post of the table it returns. But
  ##   DEqMS::spectraCounteBayes() then refits the prior from fit$sigma, fit$df.residual
  ##   and fit$count, and the reported statistic comes from that prior, so no sca. column
  ##   moves and neither does the p-value test() puts in the standardized table. Said here
  ##   rather than left to the generic message below, whose first clause -- that the method
  ##   does not fit its prior against mean feature intensity -- would be false of deqms:

  if(method %in% "deqms") {

    if(!trend) return(FALSE)

    f.msg("WARNING: test:", src, "is TRUE, and test_deqms() does hand it to",
      "limma::eBayes(), but DEqMS refits the prior from the spectra counts, so the",
      "p-value test() reports is the one a FALSE would have given;", "\n",
      " what moves is P.Value, t, B, s2.prior and s2.post of the returned table, not",
      "any sca. column;", "\n",
      " for a trended prior that changes the reported answer use test_method 'trend'",
      config=config)

    return(TRUE)
  }

  if(!trend || method %in% f.trend_methods()) return(FALSE)

  f.msg("WARNING: test:", src, "is TRUE, but test_method", method, "does not fit",
    "its variance prior against mean feature intensity, so the run is unaffected;",
    "\n", " the methods that do are", paste(f.trend_methods(), collapse=", "),
    "and 'trend', which always does", config=config)

  return(TRUE)
}

#' Hypothesis testing
#' @description
#'   Wrapper for various hypothesis testing methods.
#' @details
#'   Tests for differential expression using method specified in config. 
#'   See invididual \code{test_*} methods for more details. 
#'   The \code{method} setting meanings are: 
#'   \tabular{ll}{
#'     \code{lm}     \cr \tab Use \code{stats::lm()} on each feature. \cr
#'     \code{trend}  \cr \tab Use \code{limma::eBayes(trend=TRUE)}. \cr
#'     \code{deqms}  \cr \tab Use \code{DEqMS::spectraCounteBayes()}. \cr
#'     \code{msqrob} \cr \tab Use \code{msqrob2::msqrob()}. \cr
#'     \code{msqrob_agg} \cr \tab Use \code{msqrob2::msqrobAggregate()}: one mixed model per gene over the rows of its features, with the feature and the observation as random effects. \cr
#'     \code{proda}  \cr \tab Use \code{proDA::proDA()}. \cr
#'     \code{prolfqua} \cr \tab Use \code{prolfqua::strategy_lm()} on each feature. \cr
#'     \code{prolfqua_lmer} \cr \tab Use \code{prolfqua::strategy_lmer()}: one mixed model per gene over the rows of its features, with the feature and the observation as random effects. \cr
#'     \code{voom}   \cr \tab Use \code{limma::voom()}. \cr
#'   }
#'   Most methods are row-wise on \code{state$expression} and return one row per row of
#'     it, so their granularity follows the input. \code{"deqms"}, \code{"msqrob"},
#'     \code{"msqrob_agg"} and \code{"prolfqua_lmer"} return one row per gene whatever
#'     the input is: the first two aggregate internally, and the last two model the
#'     features of a gene instead of aggregating them, \code{"msqrob_agg"} aggregating
#'     only to carry the result. \code{h0testr::tune()} therefore does not call
#'     \code{combine_features()} before those four, and \code{h0testr::run()} should be
#'     given a \code{config$run_order} without it; \code{"prolfqua_lmer"} and
#'     \code{"msqrob_agg"} need feature level input and refuse a state whose features
#'     are already aggregated.
#'   The two feature level mixed model paths fit the same random structure by different
#'     engines. \code{"prolfqua_lmer"} takes Satterthwaite degrees of freedom for the
#'     tested contrast and is the better calibrated of the two; \code{"msqrob_agg"}
#'     reports \code{msqrob2}'s moderated t against \code{dfPosterior}, which is
#'     generous, and is mildly anti-conservative as a result. See
#'     \code{h0testr::test_msqrob()} for the simulated rejection rates.
#'   Several settings are read by some engines and not by others, so one set on a run
#'     whose method does not consult it does nothing at all:
#'     \code{config$test_random_obs} is read by \code{"prolfqua_lmer"} and
#'     \code{"msqrob_agg"}, and \code{config$test_ridge} by \code{"msqrob_agg"}. Such a
#'     setting is reported as a \code{NOTE} in the log rather than refused, and only when
#'     its value differs from the one \code{h0testr::new_config()} ships:
#'     \code{h0testr::tune()} sets one configuration and varies the method, so most
#'     combinations of a sweep carry a setting meant for one of the others. Said here,
#'     once per run, rather than in \code{h0testr::check_config()}, which every step of
#'     the workflow calls.
#'   The trended variance prior is the same kind of setting, but is resolved here for
#'     every method rather than read by each: the \code{trend} argument answers when
#'     given and \code{config$test_trend} otherwise, and the resolved value is passed to
#'     \code{"prolfqua"}, where it fits the prior against mean feature intensity, and to
#'     \code{"deqms"}, where it sets \code{limma::eBayes(trend=)} but does not change what
#'     is reported: \code{DEqMS} fits its own prior from quantities
#'     \code{limma::eBayes()} leaves alone, so only the \code{limma} columns of
#'     \code{original} move. \code{"trend"} always trends, that being its
#'     definition, so \code{TRUE} is silent there and only an explicit \code{FALSE} draws
#'     a remark. The remaining methods cannot trend: \code{"lm"} and \code{"msqrob"} do
#'     not moderate against a covariate, \code{"proda"} fits its own prior, \code{"voom"}
#'     carries the mean-variance relationship in its precision weights so a trended prior
#'     would count it twice, and the mixed paths have a Satterthwaite denominator rather
#'     than one residual variance to shrink. A \code{TRUE} reaching one of those is a
#'     \code{WARNING} in the log, not a refusal, for the same \code{h0testr::tune()}
#'     reason.
#'   See documentation for \code{h0testr::new_config()}
#'     for more detailed description of configuration parameters.
#' @param state List with elements formatted like the list returned by \code{read_data()}:
#'   \tabular{ll}{
#'     \code{expression} \cr \tab Numeric matrix with non-negative expression values. \cr
#'     \code{features}   \cr \tab A data.frame with feature meta-data for rows of expression. \cr
#'     \code{samples}    \cr \tab A data.frame with observation meta-data for columns of expression. \cr
#'   } 
#' @param config List with configuration values. Uses the following keys:
#'   \tabular{ll}{
#'     \code{feat_col}       \cr \tab Name of column in \code{state$fetaures} matching \code{rownames(state$expression)}. \cr
#'     \code{obs_col}        \cr \tab Name of column in \code{state$samples} matching \code{colnames(state$expression)}. \cr
#'     \code{gene_id_col}    \cr \tab Name of column in \code{state$fetaures} with gene/protein-group ids. \cr
#'     \code{frm}            \cr \tab Formula (formula) to be fit. \cr
#'     \code{test_term}      \cr \tab Term (character scalar) to be tested for non-zero coefficient. \cr
#'     \code{contrast}      \cr \tab Weighted sum (character scalar) of coefficients of \code{config$frm} to test instead of \code{config$test_term}; "" for none. \cr
#'     \code{reference_levels} \cr \tab Named character vector with the reference level of each factor variable in \code{config$frm} (see examples). \cr
#'     \code{test_method}    \cr \tab Character scalar in \code{h0testr::test_methods()}. \cr
#'   }
#' @param method Name of test method where
#'   \code{method \%in\% h0testr::test_methods()}. Defaults to
#'   \code{config$test_method}, and either being \code{""} counts as unset; both unset
#'   is an error. Also accepts \code{"none"}, which skips the test step and returns
#'   \code{NULL} without reading \code{state}; that is what
#'   \code{config$test_method="none"} is for in \code{h0testr::run()}, whose result then
#'   carries the processed state and \code{NULL} in place of a test. \code{"none"} is
#'   deliberately not one of \code{h0testr::test_methods()}, which names the engines;
#'   \code{h0testr::check_config()} lists what the key may hold.
#' @param is_log_transformed Logical scalar: whether \code{state$expression} has
#'   been log transformed. Only consulted for \code{method \%in\% c("proda",
#'   "prolfqua", "prolfqua_lmer")}. Defaults to \code{config$is_log_transformed}, which
#'   \code{h0testr::initialize()} and \code{h0testr::normalize()} maintain;
#'   passing both is an error unless they agree.
#' @param prior_df Prior degrees of freedom for method \code{proda};
#'   where \code{2 <= prior_df <= n_features}.
#' @param trend Logical scalar: whether the variance prior of the moderation is fitted
#'   against mean feature intensity rather than being flat. Changes the reported answer
#'   only for \code{method="prolfqua"}; \code{"trend"} always trends; \code{"deqms"}
#'   accepts it, and it moves the \code{limma} columns of \code{original}, but
#'   \code{DEqMS} refits the prior from the spectra counts afterward so the reported
#'   p-value is the one a \code{FALSE} gives; and the rest cannot trend at all. A request
#'   the resolved method cannot honor, or honors without the reported answer moving, is a
#'   \code{WARNING} rather than an error, so that an \code{h0testr::tune()} sweep over
#'   methods does not die on the first method that will not trend. Defaults to
#'   \code{config$test_trend}, and to \code{FALSE} when that is absent too; unlike
#'   \code{is_log_transformed}, an argument that disagrees with the configuration simply
#'   wins, this being a preference rather than a fact about the data. Unrelated to
#'   \code{method="trend"}, which names a different limma fit.
#' @return A list with the following elements: \cr
#'   \tabular{ll}{
#'     \code{original} \cr \tab A \code{data.frame} with results in native format returned by test. \cr
#'     \code{standard} \cr \tab A \code{data.frame} with results in a standardized format. \cr
#'     \code{fit}      \cr \tab Fitted model returned by the selected testing procedure. \cr
#'   } \cr
#'   Or \code{NULL} for \code{method="none"}, which skips the test step. \cr
#'   The \code{standard} \code{data.frame} has the following fields: \cr
#'   \tabular{ll}{
#'     \code{feature}   \cr \tab Name of feature tested. \cr
#'     \code{expr}      \cr \tab Average feature expression. \cr
#'     \code{logfc}     \cr \tab Estimated effect size; see below. \cr
#'     \code{stat}      \cr \tab Value of test statistic. \cr
#'     \code{lod}       \cr \tab Log-odds of differential expression. \cr
#'     \code{pval}      \cr \tab Raw p-value resulting from test. \cr
#'     \code{adj_pval}  \cr \tab Adjusted (for multiple testing) p-value. \cr
#'   }
#'   What \code{logfc} holds depends on how many design matrix columns
#'     \code{config$test_term} resolves to, which is a property of
#'     \code{config$frm} and \code{config$test_term} and so is the same for every
#'     row of one result:
#'   \itemize{
#'     \item One column: the coefficient of that column, signed. For a two level
#'       factor that is the difference between its levels on the scale of
#'       \code{state$expression}, so a log fold change when the input is log
#'       transformed. For a \strong{continuous} covariate it is instead the change
#'       \strong{per unit} of that covariate, which is not a fold change between
#'       groups and whose size depends on the units the covariate is recorded in:
#'       an effect per month is a twelfth of the same effect per year. Thresholds
#'       on \code{abs(logfc)} therefore have to be chosen with the covariate's
#'       units in mind, although rankings and p-values are unaffected.
#'     \item More than one column: a joint test of several coefficients has no
#'       single contrast to report, so what is reported is the total swing, the
#'       range over the observations of the fitted contribution of the terms under
#'       test. This is the largest difference those terms can account for between
#'       any two observations: for a factor of more than two levels, the largest
#'       difference between any two of its levels; for a continuous covariate, the
#'       slope times the range of the covariate, that is the total change across
#'       the observed range. It is unsigned, since several coefficients have no one
#'       direction, and it is on the same scale as the single coefficient case.
#'     \item \code{config$contrast}: the weighted sum of the coefficients it names,
#'       signed, which is the quantity being tested and is on the same scale as the
#'       single coefficient case. A contrast is one degree of freedom however many
#'       coefficients it weights, so unlike a joint test it always has one number to
#'       report. See \code{h0testr::new_config()} for how to write one.
#'   }
#'   Which side supplies \code{logfc} depends on \code{method}, although what the column
#'     means does not. On a one column test \code{"trend"}, \code{"voom"},
#'     \code{"deqms"}, \code{"msqrob"}, \code{"msqrob_agg"} and \code{"proda"} each
#'     carry an effect size in their own result table, and it is reported as it stands;
#'     \code{"lm"}, \code{"prolfqua"} and \code{"prolfqua_lmer"} do not, and no method
#'     does on a joint test or a \code{config$contrast}, so for those the value is
#'     computed here from the fitted coefficients of the columns under test. Both routes
#'     give the quantity described above, on the scale of \code{state$expression}, so
#'     \code{logfc} is comparable between methods; what differs between them is the
#'     model that produced the coefficient, not what the column holds.
#'   \code{lod} is narrower than the rest: it is \code{limma}'s \code{B}, the log-odds
#'     that a feature is differentially expressed, and only \code{method} \code{"trend"}
#'     and \code{"voom"} report it, and only on a one column test. Every other method
#'     leaves it \code{NA}, as does any joint test: \code{B} is a quantity of
#'     \code{limma}'s posterior odds calculation, nothing else here computes one, and
#'     unlike \code{logfc} and \code{expr} it is not filled in afterward. So an empty
#'     \code{lod} column says which method ran, not anything about the data.
#'   \code{feature} holds the value of \code{config$feat_col} identifying the row
#'     of \code{state$expression} the result describes, except for the gene level
#'     methods, \code{method \%in\% c("deqms", "msqrob", "msqrob_agg",
#'     "prolfqua_lmer")}, which report one row per value of
#'     \code{config$gene_id_col} whatever level the input is at. So a
#'     method other than those reports precursors when handed precursors and
#'     genes when handed the output of \code{h0testr::combine_features()}, which
#'     sets \code{config$feat_col} to \code{config$gene_id_col}. The feature
#'     metadata reported in \code{original} is at the matching level: for the gene
#'     level methods it is the gene level form of \code{state$features},
#'     built exactly as \code{h0testr::combine_features()} builds it.
#'   \code{expr} is the mean of the values handed to the test, over the
#'     observations where the feature was seen; for the gene level methods, which
#'     take feature level input and report gene level results, it
#'     is also over the features of each gene.
#'   \code{stat} is the statistic the engine itself reports, with two exceptions, both
#'     for a joint test and both for the same reason, that the engine's API answers one
#'     coefficient at a time while its fit supports the joint test.
#'     \code{method \%in\% c("msqrob", "msqrob_agg")} reports an F computed by
#'     \code{h0testr} from the fitted \code{msqrob2} models, since
#'     \code{msqrob2::hypothesisTest()} answers one contrast at a time; see
#'     \code{h0testr::test_msqrob()}. \code{method="deqms"} reports an F computed from
#'     the per-gene variance \code{DEqMS::spectraCounteBayes()} fits, since
#'     \code{DEqMS} moderates one coefficient's t-statistic and has no F-analogue; see
#'     \code{h0testr::test_deqms()}. Both reduce to the engine's own statistic at one
#'     numerator degree of freedom.
#'   A \code{config$contrast} run has no such exception: every engine reports its
#'     own statistic for a contrast, the two \code{msqrob} methods with
#'     \code{msqrob2::hypothesisTest()} and \code{"deqms"} with the single coefficient
#'     \code{limma::contrasts.fit()} leaves it.
#' @examples
#' set.seed(101)
#' ## no missing values: mnar_c0=-Inf, mnar_c1=0, mcar_p=0
#' samps <- h0testr::sim_samples(factors=list(condition=c("placebo", "drug")),
#'   n_per_cell=6)
#' sim <- h0testr::sim_design(samps, frm=~condition, test_term="condition", n_genes=25,
#'   n_genes_signif=5, effects=2, mnar_c0=-Inf, mnar_c1=0, mcar_p=0)
#' state <- sim$state
#' state$expression <- log2(state$expression + 1)
#' 
#' config <- sim$config   ## frm, test_term, the id columns and reference_levels, set
#'                        ##   to match what was simulated; save_state is FALSE
#' config$is_log_transformed <- TRUE   ## normalize() would; logged just above
#' rm(samps, sim)
#' 
#' ## set up and check covariates and parameters:
#' out <- h0testr::initialize(state, config, minimal=TRUE)
#' 
#' out <- h0testr::test(out$state, out$config, method="trend")
#' head(out$original)
#' head(out$standard)
#' summary(out$fit)

test <- function(state, config, method=NULL,
    is_log_transformed=NULL, prior_df=NULL, trend=NULL) {

  if(is.null(method) || method %in% "") method <- config$test_method
  if(is.null(method) || method %in% "") {
    f.err("test: method and config$test_method both unset", config=config)
  }

  ## "none" skips the test step, and returns before anything is read from state. The
  ##   branch used to sit at the end of the dispatch chain below, by which point the
  ##   covariates had been relevelled and f.design_test_cols() had built and rank checked
  ##   a design, so test(method="none") could fail at the step it had been told to skip,
  ##   on a state whose test term is not estimable. What makes the value worth having is
  ##   run() with config$test_method="none", as a workflow that normalizes, filters and
  ##   imputes and then stops, and a state that cannot be tested is exactly a case for
  ##   that. Not one of test_methods(), which names the engines; check_config() lists what
  ##   the key may hold:

  if(method %in% "none") {
    f.msg("skipping testing: method %in% 'none'", config=config)
    return(NULL)
  }

  if(is.null(prior_df)) prior_df <- config$test_prior_df
  if(method %in% "proda" && is.null(prior_df)) {
    f.err("test: method %in% 'proda' && is.null(prior_df)", config=config)
  }
  
  ## only these two methods are told the scale; resolved here rather than in
  ##   them so that an unusable combination is caught before the fit:

  if(method %in% c("proda", "prolfqua", "prolfqua_lmer")) {
    is_log_transformed <- f.is_log_transformed(is_log_transformed, config,
      "test")
  }
  
  ## resolved for every method rather than only the two that take it, so that the one
  ##   place the value is settled is also the place a method that cannot honor it says
  ##   so; the argument overrides config$test_trend. See f.is_trend():

  trend_given <- !(is.null(trend) || (is.character(trend) && all(trend %in% "")))
  trend <- f.is_trend(trend, config, "test")

  f.msg("test: method:", method, "; is_log_transformed:", is_log_transformed,
    "; prior_df:", prior_df, "; trend:", trend, config=config)

  ## a setting the resolved method does not read, said once here rather than once per
  ##   step by check_config(); see f.note_ignored_settings() and f.note_trend():

  f.note_ignored_settings(method, config)
  f.note_trend(method, trend, trend_given, config)

  ## the factor covariates of config$frm are rebuilt with the level ordering
  ##   config declares, before any design is derived and before the state reaches
  ##   an engine, so that every method reports the declared reference level rather
  ##   than one re-derived by sorting. initialize() has already done this, and
  ##   redoing it is a no-op there; what it is for is a direct caller of test(),
  ##   for whom config$reference_levels used to be silently ignored. It also drops
  ##   the factor levels no remaining observation is at, which is not a no-op after
  ##   initialize(): a filtering step can empty a level that was there when
  ##   initialize() recorded config$factor_levels, and model.matrix() codes such a
  ##   level as a column of zeros, which f.design_check_rank() refuses. See
  ##   f.relevel_state_covariates():

  state <- f.relevel_state_covariates(state, config, caller="test")

  ## every method below tests config$test_term against a reduced model, whether
  ##   by dropping terms or by contrasting coefficients, so none of them has
  ##   anything to report when the two models span the same space; checked once
  ##   here rather than in each method:

  design <- f.design_test_cols(state, config)

  ## the column of state$features that identifies the rows the engine will return, and
  ##   so the ids the standardized table is keyed by. Decided once here rather than
  ##   inside each f.format_*(), which used to read a column of its own choosing:
  ##   f.format_lm() read config$gene_id_col although test_lm() is row-wise on
  ##   state$expression, so a state whose features had not been aggregated came back
  ##   labelled with gene ids and failed the check below. See f.test_id_col():

  test_col <- f.test_id_col(method, config)

  if(method %in% "lm") {
    result <- test_lm(state, config)
    tbl2 <- f.format_lm(result$hits, test_col, config)
  } else if(method %in% "trend") {
    result <- test_trend(state, config)
    tbl2 <- f.format_limma(result$hits, config)
  } else if(method %in% "deqms") {
    result <- test_deqms(state, config, trend=trend)
    tbl2 <- f.format_deqms(result$hits, config)
  } else if(method %in% "msqrob") {
    result <- test_msqrob(state, config)
    tbl2 <- f.format_msqrob(result$hits, test_col, config)
  } else if(method %in% "msqrob_agg") {
    result <- test_msqrob(state, config, aggregate=TRUE)
    tbl2 <- f.format_msqrob(result$hits, test_col, config)
  } else if(method %in% "proda") {
    result <- test_proda(state, config,
      is_log_transformed=is_log_transformed, prior_df=prior_df)
    tbl2 <- f.format_proda(result$hits, config)
  } else if(method %in% "prolfqua") {
    result <- test_prolfqua(state, config,
      is_log_transformed=is_log_transformed, trend=trend)
    tbl2 <- f.format_prolfqua(result$hits, test_col, config)
  } else if(method %in% "prolfqua_lmer") {
    result <- test_prolfqua(state, config,
      is_log_transformed=is_log_transformed, mixed=TRUE)
    tbl2 <- f.format_prolfqua(result$hits, test_col, config)
  } else if(method %in% "voom") {
    result <- test_voom(state, config)
    tbl2 <- f.format_limma(result$hits, config)
  } else f.err("test: unexpected method:", method, config=config)

  ## the effect size and the average expression, where the engine's own table does not
  ##   carry them; see f.logfc_effect() for what logfc holds for a joint test:

  tbl2 <- f.fill_standard(tbl2, result, state, design, method, config)

  ## the feature metadata to report alongside the result. For a gene level method that
  ##   is the gene level form of the table, built by the same function combine_features()
  ##   uses, so that a gene row carries the same metadata whichever route produced it:

  feats <- state$features
  if(f.gene_level_method(method)) {
    feats <- f.gene_features(feats, config, "test")$features
  }
  rownames(feats) <- feats[[test_col]]

  if(!all(tbl2$feature %in% rownames(feats))) {
    i <- !(tbl2$feature %in% rownames(feats))
    f.err("test: method", method, "returned", sum(i), "of", nrow(tbl2), "results",
      "whose feature id is not in state$features[[", test_col, "]], so they cannot",
      "be matched back to the feature metadata;", "\n",
      "  first few returned:", utils::head(tbl2$feature[i], 5), "\n",
      "  first few available:", utils::head(rownames(feats), 5), "\n",
      "  config$feat_col:", config$feat_col,
      "; config$gene_id_col:", config$gene_id_col, config=config)
  }

  ## where the feature id sits in result$hits, which every engine answers
  ##   differently: the limma family puts it in the rownames, proDA::test_diff() in a
  ##   'name' column, msqrob2 and prolfqua in the id column named by config. The
  ##   f.format_*() functions each know which, but they return only the standardized
  ##   table, so the key is recovered here by taking the first candidate that
  ##   reproduces every id in it. Guessing instead is not visible in the standardized
  ##   table, only in the original one below, which comes back as a full set of NA
  ##   columns; hence taking a candidate only on a complete match, and erroring
  ##   rather than emitting that table when no candidate gives one:

  o <- NULL
  for(key in list(rownames(result$hits), result$hits[["feature"]],
    result$hits[["name"]], result$hits[[test_col]])) {

    if(is.null(key)) next
    o2 <- match(tbl2$feature, as.character(key))
    if(!any(is.na(o2))) {
      o <- o2
      break
    }
  }

  if(is.null(o)) {
    f.err("test: cannot match the standardized results back to the results",
      method, "returned, so the original results cannot be reported;", "\n",
      "feature ids:", utils::head(tbl2$feature, 5), "\n",
      "columns of the returned results:", names(result$hits), config=config)
  }

  tbl <- cbind(feats[tbl2$feature, , drop=F], result$hits[o, , drop=F])
  rownames(tbl) <- NULL
  
  if((!is.null(config$save_state)) && config$save_state) {
    
    file_out <- paste0(config$dir_out, "/", length(config$run_order) + 3, 
      config$result_mid_out, ".reformat", config$suffix_out)
    f.log("writing reformatted results to", file_out, config=config)
    f.save_tsv(tbl2, file_out, config)

    file_out <- paste0(config$dir_out, "/", length(config$run_order) + 3,
      config$result_mid_out, ".original", config$suffix_out)
    f.log("writing original results to", file_out, config=config)
    f.save_tsv(tbl, file_out, config)
  }
  
  return(list(original=tbl, standard=tbl2, fit=result$fit))
}
