## features with enough distinct values for one variable in the formula;
##   character scalar term, character scalar type in c("factor", "numeric")
##   (from f.covariate_types()), numeric matrix state$expression, data.frame
##   state$samples:

f.filter_features_by_term <- function(term, state, config, type="factor",
    n_non_na_min=2, n_distinct_min=2,
    n_groups_non_na_min=2, n_groups_distinct_min=1) {

  if(!(term %in% names(state$samples))) {
    f.err("f.filter_features_by_term: term: ", term,
      " not in names(state$samples): ", names(state$samples), config=config)
  }

  if(!(length(type) %in% 1 && type %in% c("factor", "numeric"))) {
    f.err("f.filter_features_by_term: type not scalar in c('factor',",
      "'numeric'); term:", term, "; type:", type, config=config)
  }

  if(type %in% "numeric") {

    ## a continuous variable has no groups whose members can be counted, and
    ##   counting its distinct values as groups would drop nearly every feature;
    ##   so the cheap screen is enough non-NA values, enough distinct values
    ##   among them, and at least two distinct values of the variable among the
    ##   observations that were measured. Whether the corresponding term is
    ##   actually estimable is settled by filter_features_by_estimability():

    x <- state$samples[[term]]

    f_ok_num <- function(v) {
      i <- !is.na(v)
      if(sum(i) < n_non_na_min) return(FALSE)
      if(length(unique(v[i])) < n_distinct_min) return(FALSE)
      return(length(unique(x[i])) >= 2)
    }

    i_ok <- apply(state$expression, 1, f_ok_num)

  } else {

    f_ok_fac <- function(v, term, meta) {

      ## n_non_na_min non-NA values in each of n_groups_non_na_min groups:
      i <- !is.na(v)
      mat <- as.matrix(table(meta[[term]][i]))
      check1 <- sum(mat[, 1, drop=T] >= n_non_na_min, na.rm=T) >= n_groups_non_na_min

      ## n_distinct_min unique non-NA values in n_groups_distinct_min groups:
      i <- !(is.na(v) | duplicated(v))
      mat <- as.matrix(table(meta[[term]][i]))
      check2 <- sum(mat[, 1, drop=T] >= n_distinct_min, na.rm=T) >= n_groups_distinct_min

      return(check1 && check2)
    }

    i_ok <- apply(state$expression, 1, f_ok_fac, term, state$samples)
  }

  return(i_ok)
}

#' Filter features without enough distinct values for each variable in \code{config$frm}.
#' @description
#'   Filter features without enough distinct values for each variable in \code{config$frm}.
#' @details
#'   Every variable in \code{config$frm} is screened, including variables that
#'     appear only within an interaction term. Interaction terms themselves are
#'     not screened here.
#'   This is a cheap pre-screen, intended to remove obviously untestable
#'     features before more expensive steps; whether a term is actually
#'     estimable for a feature is a separate question.
#'   For a factor variable, a feature is kept if at least
#'     \code{n_groups_non_na_min} levels have at least \code{n_non_na_min}
#'     non-\code{NA} values, and at least \code{n_groups_distinct_min} levels
#'     have at least \code{n_distinct_min} distinct non-\code{NA} values.
#'   For a numeric (continuous) variable there are no levels to count, so a
#'     feature is kept if it has at least \code{n_non_na_min} non-\code{NA}
#'     values, at least \code{n_distinct_min} distinct values among them, and
#'     the variable takes at least two distinct values across the observations
#'     where the feature was measured.
#'   Variables are classified as factor or numeric as described for
#'     \code{h0testr::initialize()}, and their values are checked as described
#'     there: missing, blank, non-finite and constant covariate values are
#'     errors here too, since the counts below would otherwise be taken over the
#'     observations whose covariates happen to be known, while every step that
#'     builds a design matrix refuses the same \code{state}.
#' @param state A list with elements like that returned by \code{read_data()}:
#'   \tabular{ll}{
#'     \code{expression} \cr \tab Numeric matrix with non-negative expression values. \cr
#'     \code{features}   \cr \tab A data.frame with feature meta-data for rows of expression. \cr
#'     \code{samples}    \cr \tab A data.frame with observation meta-data for columns of expression. \cr
#'   } 
#' @param config List with configuration values. Uses the following keys:
#'   \tabular{ll}{
#'     \code{frm}             \cr \tab Formula object specifying model to be fitted. \cr
#'     \code{reference_levels} \cr \tab Named character vector with the reference level of each factor variable in \code{config$frm}. \cr
#'     \code{covariate_types} \cr \tab Optional; classification of variables in \code{config$frm}, as set by \code{initialize()}. \cr
#'   }
#' @param n_non_na_min Minimum number of non-NA values per feature. Non-negative integer.
#' @param n_distinct_min Minimum number of distinct non-NA values per feature. Non-negative integer.
#' @param n_groups_non_na_min Minimum number of factor levels meeting \code{n_non_na_min}. Non-negative integer. Ignored for numeric variables.
#' @param n_groups_distinct_min Minimum number of factor levels meeting \code{n_distinct_min}. Non-negative integer. Ignored for numeric variables.
#' @return An updated \code{state} list with the following elements:
#'   \tabular{ll}{
#'     \code{expression} \cr \tab Numeric matrix with non-negative expression values. \cr
#'     \code{features}   \cr \tab A data.frame with feature meta-data for rows of expression. \cr
#'     \code{samples}    \cr \tab A data.frame with observation meta-data for columns of expression. \cr
#'   } 
#' @examples
#' set.seed(101)
#' exprs <- h0testr::sim1(n_obs=6, n_feats=12, mcar_p=0.5)$mat
#' feats <- data.frame(feature_id=rownames(exprs))
#' samps <- data.frame(observation_id=colnames(exprs), age=c(rep(4, 3), rep(12, 3)))
#' state <- list(expression=exprs, features=feats, samples=samps)
#' config <- list(frm=~age)
#' state2 <- h0testr::filter_features_by_formula(state, config)
#' print(state)
#' print(state2)

filter_features_by_formula <- function(state, config, 
    n_non_na_min=2, n_distinct_min=2, 
    n_groups_non_na_min=2, n_groups_distinct_min=1) {
  
  if(is.null(config$frm)) {
    f.err("filter_features_by_formula: is.null(config$frm)", config=config)
  }
  ## screened one variable at a time, including variables that appear only
  ##   within an interaction; interaction terms themselves are left to
  ##   filter_features_by_estimability(), since they have no cheap screen. The
  ##   criteria depend on whether the variable is a factor or continuous:

  types <- f.covariate_types(state, config)

  ## the counts below are per level of a factor and per observation for a continuous
  ##   variable, neither of which needs a design matrix, so nothing here would notice
  ##   a covariate that no fit can use: base::table() drops missing values silently,
  ##   so a feature would be screened over the observations whose covariates happen
  ##   to be known and then reported as kept, while every downstream step that builds
  ##   a design refuses the same state. Checked here so that the count this reports
  ##   is over the observations that would actually be fit:

  f.check_covariate_values(state, config, types=types,
    caller="filter_features_by_formula", warn_distinct=FALSE)

  vars <- names(types)

  result <- matrix(TRUE, nrow=nrow(state$expression), ncol=length(vars),
    dimnames=list(NULL, vars))

  for(idx in seq_along(vars)) {
    result[, idx] <- f.filter_features_by_term(
      vars[idx], state, config,
      type=types[[vars[idx]]],
      n_non_na_min=n_non_na_min,
      n_distinct_min=n_distinct_min,
      n_groups_non_na_min=n_groups_non_na_min,
      n_groups_distinct_min=n_groups_distinct_min
    )
  }

  i <- apply(result, 1, all, na.rm=F)
  i[is.na(i)] <- F
  
  f.msg("filtering", sum(!i), "features by formula; keeping", sum(i), config=config)
  state$expression <- state$expression[i, , drop=F]
  state$features <- state$features[i, , drop=F]
  ## state$samples <- state$samples
  
  return(state)
}

#' Filter features whose model coefficients are not estimable.
#' @description
#'   Filter features for which the coefficients of \code{config$test_term}
#'     cannot be estimated from the values actually observed for that feature.
#' @details
#'   This is the last and most expensive of three tiers of filtering. The
#'     earlier tiers ask necessary but not sufficient questions:
#'     \code{prefilter()} counts non-\code{NA} values without reference to
#'     \code{config$frm}, and \code{filter_features_by_formula()} screens one
#'     variable at a time, marginally. Only this function looks at the design
#'     matrix as a whole, so only this function can answer whether the requested
#'     test is estimable for a feature.
#'   A value is missing if and only if it is \code{NA}; see
#'     \code{h0testr::initialize()}. For each feature, let \code{S} be the
#'     observations where it is not \code{NA}, \code{X} the design matrix built
#'     from \code{config$frm}, and \code{X_red} the same matrix with the columns
#'     of \code{config$test_term} removed. Then:
#'     \tabular{ll}{
#'       \code{df_test}   \cr \tab \code{rank(X[S, ]) - rank(X_red[S, ])}; estimable degrees of freedom for \code{config$test_term}. \cr
#'       \code{df_resid}  \cr \tab \code{length(S) - rank(X[S, ])}; residual degrees of freedom. \cr
#'       \code{df_intend} \cr \tab \code{df_test} recomputed over all observations; the test that was asked for. \cr
#'       \code{df_deficit} \cr \tab \code{ncol(X) - rank(X[S, ])}; coefficients of the requested model that are not estimable. \cr
#'     }
#'   \code{X} is built once over all observations and then subset by row, so the
#'     factor level ordering set by \code{initialize()} is preserved: a level
#'     absent from \code{S} leaves an all-zero column, which is exactly the rank
#'     deficiency to be detected. Transformations in \code{config$frm} are not
#'     supported (see \code{h0testr::initialize()}), so row subsetting and
#'     per-feature refitting agree.
#'   The reduced model is formed by dropping term labels, not by editing the
#'     formula text, and the labels come from the same helper the hypothesis
#'     tests use. A \code{config$test_term} naming a variable therefore drops
#'     every term containing that variable, keeping the reduced model
#'     hierarchical; a \code{config$test_term} naming an interaction drops that
#'     term alone.
#'   \code{config$estimability} selects one of three nested requirements, each
#'     strictly stronger than the one before it. \code{df_deficit} is the sum of
#'     the nuisance-side deficit and \code{df_intend - df_test}, so requiring
#'     \code{"full"} entails \code{"term"}, which entails \code{"test"}:
#'     \tabular{ll}{
#'       \code{"test"} \cr \tab \code{df_test >= 1}; the term is testable. The hypothesis tested, and the covariate adjustment applied, may differ between features. \cr
#'       \code{"term"} \cr \tab \code{df_test == df_intend}; every feature is tested against the same hypothesis, but the covariate adjustment may still differ. \cr
#'       \code{"full"} \cr \tab \code{df_deficit == 0}; every coefficient of the requested model is estimable for every feature. \cr
#'     }
#'   Independently of that, a feature is dropped when
#'     \code{df_resid < df_resid_min}, which is a question of residual precision
#'     rather than of estimability.
#'   Counts of dropped features, broken out by reason, are logged. Features whose
#'     covariate columns collapsed are counted in the log even when
#'     \code{config$estimability} is too permissive to drop them: such a feature
#'     is tested without the adjustment that was asked for, which
#'     \code{df_test} cannot reveal, since removing the test columns removes the
#'     same rank from both models.
#'   Ranks are computed once per distinct missingness pattern rather than once
#'     per feature, which on real data is usually a large saving.
#'   Following \code{add_filter_stats()}, \code{df_test} and \code{df_resid} are
#'     written into \code{state$features} for the features that survive.
#'   See documentation for \code{h0testr::new_config()}
#'     for more detailed description of configuration parameters.
#' @param state A list with elements like that returned by \code{read_data()}:
#'   \tabular{ll}{
#'     \code{expression} \cr \tab Numeric matrix with expression values. \cr
#'     \code{features}   \cr \tab A data.frame with feature meta-data for rows of expression. \cr
#'     \code{samples}    \cr \tab A data.frame with observation meta-data for columns of expression. \cr
#'   }
#' @param config List with configuration values. Uses the following keys:
#'   \tabular{ll}{
#'     \code{frm}           \cr \tab Formula object specifying model to be fitted. \cr
#'     \code{test_term}     \cr \tab Term (character) in \code{config$frm} to test for significance. \cr
#'     \code{contrast}      \cr \tab Weighted sum (character scalar) of coefficients of \code{config$frm} to test instead of \code{config$test_term}; "" for none. \cr
#'     \code{estimability}  \cr \tab Requirement placed on \code{config$test_term}; scalar character in \code{c("test", "term", "full")}. \cr
#'     \code{df_resid_min}  \cr \tab Minimum residual degrees of freedom (non-negative numeric) to keep feature. \cr
#'     \code{df_test_col}   \cr \tab Name (character) of new column in feature metadata to hold \code{df_test}. \cr
#'     \code{df_resid_col}  \cr \tab Name (character) of new column in feature metadata to hold \code{df_resid}. \cr
#'     \code{n_samples_min} \cr \tab Optional; only used to warn when it is inconsistent with \code{df_resid_min}. \cr
#'   }
#' @param estimability Requirement placed on \code{config$test_term}; scalar character in \code{c("test", "term", "full")}. Overrides \code{config$estimability}.
#' @param df_resid_min Minimum residual degrees of freedom to keep feature. Non-negative numeric. Overrides \code{config$df_resid_min}.
#' @return An updated \code{state} list with the following elements:
#'   \tabular{ll}{
#'     \code{expression} \cr \tab Numeric matrix with expression values. \cr
#'     \code{features}   \cr \tab A data.frame with feature meta-data for rows of expression, with \code{df_test} and \code{df_resid} added. \cr
#'     \code{samples}    \cr \tab A data.frame with observation meta-data for columns of expression. \cr
#'   }
#' @examples
#' set.seed(101)
#' exprs <- h0testr::sim1(n_obs=8, n_feats=12, mcar_p=0.4)$mat
#' feats <- data.frame(feature_id=rownames(exprs))
#' samps <- data.frame(observation_id=colnames(exprs), grp=rep(c("ctl", "trt"), 4))
#' state <- list(expression=exprs, features=feats, samples=samps)
#' config <- list(frm=~grp, test_term="grp", reference_levels=c(grp="ctl"),
#'   estimability="test", df_resid_min=2, df_test_col="df_test",
#'   df_resid_col="df_resid")
#' state2 <- h0testr::filter_features_by_estimability(state, config)
#' print(nrow(state$expression))
#' print(state2$features)

filter_features_by_estimability <- function(state, config, estimability=NULL,
    df_resid_min=NULL) {

  if(!is.matrix(state$expression)) {
    f.err("filter_features_by_estimability: !is.matrix(state$expression);",
      "class(state$expression):", class(state$expression), config=config)
  }

  if(is.null(config$frm)) {
    f.err("filter_features_by_estimability: is.null(config$frm)", config=config)
  }

  if(is.null(estimability)) estimability <- config$estimability
  if(is.null(estimability) || estimability %in% "") {
    estimability <- "test"
    f.msg("filter_features_by_estimability: estimability and",
      "config$estimability both unset; using default:", estimability,
      config=config)
  }

  allowed <- c("test", "term", "full")
  if(!(length(estimability) %in% 1 && estimability %in% allowed)) {
    f.err("filter_features_by_estimability: estimability not scalar in",
      paste0("c('", paste(allowed, collapse="', '"), "')"), "; value:",
      estimability, config=config)
  }

  if(is.null(df_resid_min)) df_resid_min <- config$df_resid_min
  if(is.null(df_resid_min)) {
    df_resid_min <- 2
    f.msg("filter_features_by_estimability: df_resid_min and",
      "config$df_resid_min both unset; using default:", df_resid_min,
      config=config)
  }

  if(!(is.numeric(df_resid_min) && length(df_resid_min) %in% 1 &&
      df_resid_min >= 0)) {
    f.err("filter_features_by_estimability: df_resid_min not a non-negative",
      "numeric scalar; value:", df_resid_min, config=config)
  }

  ## the design, the columns carrying the test, and the df of the test that was
  ##   asked for, measured over all observations; derived by the same helper the
  ##   hypothesis tests use, so the reduced model screened here is the one that
  ##   will actually be tested. Throws an informative error if config$test_term
  ##   does not fit config$frm, or if the test it names is vacuous. The covariates
  ##   are releveled first, for the same reason test() relevels them, and into a
  ##   state of their own: which features are estimable has to be decided on the
  ##   design that will be fitted, but this function returns a filtered state
  ##   rather than a refitted one, so the classes of the columns it was given are
  ##   left alone. See f.relevel_state_covariates():

  st_fit <- f.relevel_state_covariates(state, config,
    caller="filter_features_by_estimability")
  design <- f.design_test_cols(st_fit, config)
  X <- design$X
  X_red <- design$X_red
  cols_test <- design$cols_test
  rank_all <- design$rank_all
  df_intend <- design$df_intend

  f.msg("filter_features_by_estimability: estimability:", estimability,
    "; df_resid_min:", df_resid_min, ";",
    if(is.null(design$contrast)) paste("test_term:", config$test_term)
      else paste("contrast:", config$contrast),
    "; design columns:", ncol(X), "; test columns:", length(cols_test),
    "; df_intend:", df_intend, config=config)

  ## n_samples_min screens on the same axis as df_resid_min, so an
  ##   n_samples_min below what df_resid_min implies does no work ahead of this
  ##   filter; not an error, since the two are set independently:

  n_needed <- df_resid_min + rank_all
  if(length(config$n_samples_min) %in% 1 && config$n_samples_min < n_needed) {
    f.msg("WARNING: filter_features_by_estimability: config$n_samples_min",
      config$n_samples_min, "is below the", n_needed, "observations a feature",
      "needs for df_resid >=", df_resid_min, "with a full rank design;", "\n",
      "the n_samples_min screen in filter_features() therefore does no work",
      "ahead of this filter", config=config)
  }

  ## ranks depend only on which observations are missing, so compute them once
  ##   per distinct missingness pattern:

  na_mat <- is.na(state$expression)
  keys <- apply(na_mat, 1, function(v) paste0(as.integer(v), collapse=""))
  u_keys <- unique(keys)
  i_rep <- match(u_keys, keys)

  f.msg("filter_features_by_estimability:", nrow(state$expression),
    "features in", length(u_keys), "distinct missingness patterns",
    config=config)

  df_test_u <- integer(length(u_keys))
  df_resid_u <- integer(length(u_keys))
  df_deficit_u <- integer(length(u_keys))
  t_last <- Sys.time()

  for(idx in seq_along(u_keys)) {

    i_obs <- !na_mat[i_rep[idx], ]
    r_full <- f.design_rank(X[i_obs, , drop=F])

    ## the reduced model comes from f.design_test_cols() rather than by dropping
    ##   columns here, so that a config$contrast run screens features against the
    ##   constrained design the tests will compare against, which is a
    ##   re-parameterization of X and not a subset of its columns:

    df_test_u[idx] <- r_full - f.design_rank(X_red[i_obs, , drop=F])
    df_resid_u[idx] <- sum(i_obs) - r_full
    df_deficit_u[idx] <- ncol(X) - r_full

    if(as.numeric(difftime(Sys.time(), t_last, units="secs")) >= 15) {
      f.msg("filter_features_by_estimability: pattern", idx, "of",
        length(u_keys), config=config)
      t_last <- Sys.time()
    }
  }

  idxs <- match(keys, u_keys)
  df_test <- df_test_u[idxs]
  df_resid <- df_resid_u[idxs]
  df_deficit <- df_deficit_u[idxs]

  ok_test <- df_test >= 1
  ok_term <- df_test >= df_intend
  ok_full <- df_deficit %in% 0
  ok_resid <- df_resid >= df_resid_min

  i <- ok_test & ok_resid
  if(estimability %in% c("term", "full")) i <- i & ok_term
  if(estimability %in% "full") i <- i & ok_full

  ## every dropped feature attributed to the first applicable reason, so the
  ##   counts sum to the number dropped:

  why <- rep("", length(i))
  why[!i & !ok_test] <- "test_term_not_estimable"
  why[!i & why %in% "" & !ok_resid] <- "too_few_residual_df"
  why[!i & why %in% "" & !ok_term] <- "test_term_partly_estimable"
  why[!i & why %in% ""] <- "model_not_full_rank"

  f.msg("filtering", sum(!i), "features by estimability; keeping", sum(i),
    config=config)
  for(nom in c("test_term_not_estimable", "too_few_residual_df",
      "test_term_partly_estimable", "model_not_full_rank")) {
    f.msg("  dropped,", nom, ":", sum(why %in% nom), config=config)
  }

  ## reported whether or not they were dropped: a feature that keeps df_test but
  ##   loses covariate columns is tested without the adjustment that was asked
  ##   for, and df_test cannot show that:

  f.msg("  features with test_term only partly estimable:", sum(!ok_term),
    "; features whose requested model is not full rank:", sum(!ok_full),
    config=config)
  f.msg("  df_test:", paste(names(table(df_test)), table(df_test), sep=":"),
    config=config)
  f.quantile(df_resid, config, digits=0)

  df_test_col <- config$df_test_col
  if(is.null(df_test_col) || df_test_col %in% "") df_test_col <- "df_test"
  df_resid_col <- config$df_resid_col
  if(is.null(df_resid_col) || df_resid_col %in% "") df_resid_col <- "df_resid"

  state$features[[df_test_col]] <- df_test
  state$features[[df_resid_col]] <- df_resid

  state$expression <- state$expression[i, , drop=F]
  state$features <- state$features[i, , drop=F]
  ## state$samples <- state$samples

  return(state)
}

#' Filter features based on number of samples
#' @description
#'   Filter features based on number of samples in which feature was measured.
#' @details
#'   A value is missing if and only if it is \code{NA}. Raw zeros are converted
#'     to \code{NA} by \code{h0testr::initialize()}, so a feature is counted as
#'     measured in a sample whenever its value there is not \code{NA},
#'     regardless of sign.
#'   Feature constant if \code{length(unique(expression_values)) \%in\% 1}.
#'   See documentation for \code{h0testr::new_config()} 
#'     for more detailed description of configuration parameters. 
#' @param state A list with elements like that returned by \code{read_data()}:
#'   \tabular{ll}{
#'     \code{expression} \cr \tab Numeric matrix with non-negative expression values. \cr
#'     \code{features}   \cr \tab A data.frame with feature meta-data for rows of expression. \cr
#'     \code{samples}    \cr \tab A data.frame with observation meta-data for columns of expression. \cr
#'   } 
#' @param config List with configuration values. Uses the following keys:
#'   \tabular{ll}{
#'     \code{n_samples_min} \cr \tab Minimum number (non-negative numeric) of samples with a non-NA value for the feature to keep feature. \cr
#'   }
#' @param n_samples_min Minimum number of samples with a non-NA value for the feature. Non-negative numeric.
#' @param remove_constant Logical scalar: if constant features of \code{state$expression} should be removed.
#' @param filter_by_formula Logical scalar: if \code{filter_features_by_formula()} should be run after other filters.
#' @return An updated \code{state} list with the following elements:
#'   \tabular{ll}{
#'     \code{expression} \cr \tab Numeric matrix with non-negative expression values. \cr
#'     \code{features}   \cr \tab A data.frame with feature meta-data for rows of expression. \cr
#'     \code{samples}    \cr \tab A data.frame with observation meta-data for columns of expression. \cr
#'   } 
#' @examples
#' set.seed(101)
#' exprs <- h0testr::sim1(n_obs=6, n_feats=12, mcar_p=0.5)$mat
#' feats <- data.frame(feature_id=rownames(exprs))
#' samps <- data.frame(observation_id=colnames(exprs), age=c(rep("4m", 3), rep("12m", 3)))
#' state <- list(expression=exprs, features=feats, samples=samps)
#' ## age is character, so its reference level has to be declared:
#' config <- list(n_samples_min=3, frm=~age, reference_levels=c(age="4m"))
#' state2 <- h0testr::filter_features(state, config)
#' print(state)
#' print(state2)

filter_features <- function(state, config, 
    n_samples_min=NULL, remove_constant=TRUE, filter_by_formula=TRUE) {
  
  check_config(config)
  
  if(!is.matrix(state$expression)) {
    f.err("filter_features: !is.matrix(state$expression)", "\n",
      "class(state$expression):", class(state$expression), config=config)
  }
  
  if(is.null(n_samples_min)) n_samples_min <- config$n_samples_min
  if(is.null(n_samples_min)) {
    f.err("filter_features: n_samples_min and config$n_samples_min both unset", 
      config=config)
  }
  
  ## NA is the only indicator of a missing value; see f.zeros_to_na():

  f <- function(v) sum(!is.na(v)) >= n_samples_min
  i <- apply(state$expression, 1, f)
  
  f.msg("filtering", sum(!i), "features, keeping", sum(i), config=config)
  state$expression <- state$expression[i, , drop=F]
  state$features <- state$features[i, , drop=F]
  f.msg("features filtered: nrow(state$expression):", nrow(state$expression), 
    "; ncol(state$features):", ncol(state$features), config=config)
  
  if(remove_constant) {
    f <- function(v) {
      i_na <- is.na(v)
      if(any(!i_na)) {
        return(length(unique(v[!i_na])) >= 2)
      } else {
        return(FALSE)
      }
    }
    i <- apply(state$expression, 1, f)
    i[is.na(i)] <- F
    f.msg("filtering", sum(!i), "constant features; keeping", sum(i), config=config)
    state$expression <- state$expression[i, , drop=F]
    state$features <- state$features[i, , drop=F]
    f.msg("constant features removed:", sum(!i), config=config)
    f.msg("non-constant features left:", nrow(state$features), config=config)
  }
  
  if(filter_by_formula) {
    state <- filter_features_by_formula(state, config)
  }
  
  return(state)
}

#' Filter samples based on number of features
#' @description
#'   Filter samples based on number of features with a non-\code{NA} value.
#' @details
#'   A value is missing if and only if it is \code{NA}. Raw zeros are converted
#'     to \code{NA} by \code{h0testr::initialize()}, so a feature is counted as
#'     measured in a sample whenever its value there is not \code{NA},
#'     regardless of sign.
#'   Sample constant if \code{length(unique(expression_values)) \%in\% 1)}.
#'   See documentation for \code{h0testr::new_config()} 
#'     for more detailed description of configuration parameters. 
#' @param state A list with elements like that returned by \code{read_data()}:
#'   \tabular{ll}{
#'     \code{expression} \cr \tab Numeric matrix with non-negative expression values. \cr
#'     \code{features}   \cr \tab A data.frame with feature meta-data for rows of expression. \cr
#'     \code{samples}    \cr \tab A data.frame with observation meta-data for columns of expression. \cr
#'   } 
#' @param config List with configuration values. Uses the following keys:
#'   \tabular{ll}{
#'     \code{n_features_min} \cr \tab Minimum number (non-negative numeric) of features with a non-NA value in the observation to keep observation. \cr
#'   }
#' @param n_features_min Minimum number of features with a non-NA value per sample. Non-negative numeric.
#' @param remove_constant Logical scalar: if constant observations of \code{state$expression} should be removed.
#' @return An updated \code{state} list with the following elements:
#'   \tabular{ll}{
#'     \code{expression} \cr \tab Numeric matrix with non-negative expression values. \cr
#'     \code{features}   \cr \tab A data.frame with feature meta-data for rows of expression. \cr
#'     \code{samples}    \cr \tab A data.frame with observation meta-data for columns of expression. \cr
#'   } 
#' @examples
#' set.seed(101)
#' exprs <- h0testr::sim1(n_obs=6, n_feats=12, mcar_p=0.5)$mat
#' feats <- data.frame(feature_id=rownames(exprs))
#' samps <- data.frame(observation_id=colnames(exprs))
#' state <- list(expression=exprs, features=feats, samples=samps)
#' config <- list(n_features_min=4)
#' state2 <- h0testr::filter_observations(state, config)
#' print(state)
#' print(state2)

filter_observations <- function(state, config, 
    n_features_min=NULL, remove_constant=TRUE) {
  
  check_config(config)
  
  if(!is.matrix(state$expression)) {
    f.err("filter_observations: !is.matrix(state$expression)", "\n",
      "class(state$expression):", class(state$expression), config=config)
  }
  
  if(is.null(n_features_min)) n_features_min <- config$n_features_min
  if(is.null(n_features_min)) {
    f.err("filter_observations: n_features_min and config$n_features_min both unset", 
      config=config)
  }
  
  ## NA is the only indicator of a missing value; see f.zeros_to_na():

  f <- function(v) sum(!is.na(v)) >= n_features_min
  i <- apply(state$expression, 2, f)
  f.msg("filtering", sum(!i), "observations, keeping", sum(i), config=config)
  state$expression <- state$expression[, i, drop=F]
  state$samples <- state$samples[i, , drop=F]
  f.msg("n_features_min:", n_features_min, config=config)
  f.msg("observations filtered:", sum(!i), config=config)
  f.msg("observations left:", ncol(state$expression), config=config)
  
  if(remove_constant) {
    f <- function(v) {
      i_na <- is.na(v)
      if(any(!i_na)) {
        return(length(unique(v[!i_na])) >= 2)
      } else {
        return(FALSE)
      }
    }
    i <- apply(state$expression, 2, f)
    f.msg("filtering", sum(!i), "constant observations, keeping", sum(i), config=config)
    state$expression <- state$expression[, i, drop=F]
    state$samples <- state$samples[i, , drop=F]
    f.msg("constant observations removed:", sum(!i), config=config)
    f.msg("non-constant observations left:", ncol(state$expression), config=config)
  }
  
  return(state)
}

#' Number of samples in which each feature was measured
#' @description
#'   Calculates the number of samples in which each feature was measured
#' @details A feature counts as measured in a sample whenever its value there is
#'   not \code{NA}, regardless of sign; raw zeros are converted to \code{NA} by
#'   \code{h0testr::initialize()}.
#' @param state A list with elements like that returned by \code{read_data()}:
#'   \tabular{ll}{
#'     \code{expression} \cr \tab Numeric matrix with non-negative expression values. \cr
#'     \code{features}   \cr \tab A data.frame with feature meta-data for rows of expression. \cr
#'     \code{samples}    \cr \tab A data.frame with observation meta-data for columns of expression. \cr
#'   } 
#' @param config List with configuration values. Does not use any params, so can pass empty list.
#' @return A numeric vector of length \code{nrow(state$expression)} with non-negative sample 
#'   counts for each feature.
#' @examples
#' set.seed(101)
#' exprs <- h0testr::sim1(n_obs=6, n_feats=12, mcar_p=0.5)$mat
#' feats <- data.frame(feature_id=rownames(exprs))
#' samps <- data.frame(observation_id=colnames(exprs))
#' state <- list(expression=exprs, features=feats, samples=samps)
#' config <- list()
#' h0testr::samples_per_feature(state, config)

samples_per_feature <- function(state, config) {
  
  if(!is.matrix(state$expression)) {
    f.err("samples_per_feature: !is.matrix(state$expression)", "\n",
      "class(state$expression):", class(state$expression), config=config)
  }
  
  ## NA is the only indicator of a missing value; see f.zeros_to_na():

  n <- apply(state$expression, 1, function(v) sum(!is.na(v)))

  return(n)
}

#' Median expression of each feature in each expressing sample
#' @description
#'   Calculates the median expression of each feature in each expressing sample. 
#' @details 
#'   Median over the non-\code{NA} values of the feature. Values are not
#'     screened by sign, so a transformed value of zero or below counts toward
#'     the median; raw zeros are already \code{NA} by this point, having been
#'     converted by \code{h0testr::initialize()}. A feature measured in no
#'     sample has no median, and gets \code{NA}.
#' @param state A list with elements like that returned by \code{read_data()}:
#'   \tabular{ll}{
#'     \code{expression} \cr \tab Numeric matrix with non-negative expression values. \cr
#'     \code{features}   \cr \tab A data.frame with feature meta-data for rows of expression. \cr
#'     \code{samples}    \cr \tab A data.frame with observation meta-data for columns of expression. \cr
#'   } 
#' @param config List with configuration values. Does not use any params, so can pass empty list.
#' @return A numeric vector of length \code{nrow(state$expression)} with median
#'   expression over the expressing samples, and \code{NA} for any feature with
#'   no expressing samples.
#' @examples
#' set.seed(101)
#' exprs <- h0testr::sim1(n_obs=6, n_feats=12)$mat
#' feats <- data.frame(feature_id=rownames(exprs))
#' samps <- data.frame(observation_id=colnames(exprs))
#' state <- list(expression=exprs, features=feats, samples=samps)
#' config <- list()
#' h0testr::feature_median_expression(state, config)

feature_median_expression <- function(state, config) {
  
  if(!is.matrix(state$expression)) {
    f.err("feature_median_expression: !is.matrix(state$expression)", "\n",
      "class(state$expression):", class(state$expression), config=config)
  }
  
  ## a feature measured in no sample has no median, and NA says so. Standing a
  ##   number in for it would place a never-measured feature at a real position
  ##   on the scale, in a statistic reported to the user; on a log scale zero is
  ##   not even a neutral choice, but near the top of the range. Nothing in the
  ##   package reads this column back, and combine_features() carries feature
  ##   metadata forward by taking the first row of each group rather than by
  ##   arithmetic, so an NA rides through untouched:

  m <- apply(state$expression, 1, stats::median, na.rm=T)
  
  return(m)
}

#' Number of measured features per sample
#' @description
#'   Calculates the number of measured features in each sample.
#' @details 
#'   A feature counts as measured in a sample whenever its value there is not
#'     \code{NA}, regardless of sign; raw zeros are converted to \code{NA} by
#'     \code{h0testr::initialize()}.
#' @param state A list with elements like that returned by `read_data()`:
#'   \tabular{ll}{
#'     \code{expression} \cr \tab Numeric matrix with non-negative expression values. \cr
#'     \code{features}   \cr \tab A data.frame with feature meta-data for rows of expression. \cr
#'     \code{samples}    \cr \tab A data.frame with observation meta-data for columns of expression. \cr
#'   }
#' @param config List with configuration values. Does not use any params, 
#'   so can pass empty list.
#' @return A numeric vector of length \code{ncol(state$expression)} with 
#'   number of features measured in each sample.
#' @examples
#' set.seed(101)
#' exprs <- h0testr::sim1(n_obs=6, n_feats=12)$mat
#' feats <- data.frame(feature_id=rownames(exprs))
#' samps <- data.frame(observation_id=colnames(exprs))
#' state <- list(expression=exprs, features=feats, samples=samps)
#' config <- list()
#' h0testr::features_per_sample(state, config)

features_per_sample <- function(state, config) {
  
  if(!is.matrix(state$expression)) {
    f.err("features_per_sample: !is.matrix(state$expression)", "\n",
      "class(state$expression):", class(state$expression), config=config)
  }
  
  ## NA is the only indicator of a missing value; see f.zeros_to_na():

  n <- apply(state$expression, 2, function(v) sum(!is.na(v)))

  return(n)
}

## prefilter features; ensure possibility (not guarantee) of 2+ distinct non-NA values 
##   in one group and 2+ non-NA values in another group:

f.prefilter_features <- function(state, min1=3, min2=4) {
  
  n1 <- apply(state$expression, 1, function(v) sum(!(is.na(v) | duplicated(v)), na.rm=T))
  i1 <- n1 >= min1
  
  n2 <- apply(state$expression, 1, function(v) sum(!(is.na(v)), na.rm=T))
  i2 <- n2 >= min2
  
  i <- i1 & i2
  state$expression <- state$expression[i, , drop=F]
  state$features <- state$features[i, , drop=F]
  ## state$samples <- state$samples
  
  return(state)
}

#' Prefilter data
#' @description
#'   Remove features and observations too sparse to be worth carrying through
#'     normalization and aggregation.
#' @details
#'   A cheap first pass, meant to run before anything has been normalized. It removes
#'     features with fewer than 3 distinct non-\code{NA} values or fewer than 4
#'     non-\code{NA} values, which leaves the possibility, though not the guarantee, of
#'     two distinct values in one group and two values in another, and then removes
#'     observations with fewer than \code{n_features_min} non-\code{NA} values. Neither
#'     screen is a statement about \code{config$frm}: the filtering that knows the design
#'     is \code{h0testr::filter()}, which runs after normalization and aggregation.
#'   \code{config$n_samples_min} and \code{config$n_features_min} are deliberately not
#'     consulted here. They are the thresholds of \code{h0testr::filter()}, and they are
#'     meant to be applied to the aggregated data that step sees, not to the precursor
#'     level table this one gets: \code{config$n_features_min} in particular defaults to
#'     1000, which is a sensible count of protein groups per observation and would
#'     discard nearly every observation if applied here. The feature thresholds above are
#'     therefore fixed, and \code{n_features_min} is this function's own argument.
#'   A value is missing if and only if it is \code{NA}; raw zeros are converted to
#'     \code{NA} by \code{h0testr::initialize()}, which runs first.
#'   Reports the state before and after each of the two screens.
#'   See documentation for \code{h0testr::new_config()}
#'     for more detailed description of configuration parameters.
#' @param state List with elements formatted like the list returned by \code{read_data()}:
#'   \tabular{ll}{
#'     \code{expression} \cr \tab Numeric matrix with non-negative expression values. \cr
#'     \code{features}   \cr \tab A data.frame with feature meta-data for rows of expression. \cr
#'     \code{samples}    \cr \tab A data.frame with observation meta-data for columns of expression. \cr
#'   } 
#' @param config List with configuration values. Requires the following keys:
#'   \tabular{ll}{
#'     \code{feat_col}  \cr \tab Name of column (character) in \code{state$features} matching \code{rownames(state$expression)}. \cr
#'     \code{obs_col}   \cr \tab Name of column (character) in \code{state$samples} matching \code{colnames(state$expression)}. \cr
#'   }
#' @param n_features_min Minimum number of features with a non-NA value per observation;
#'   numeric >= 2. This function's own threshold; \code{config$n_features_min} is not
#'   consulted, being the threshold of \code{h0testr::filter()}. There is no
#'   corresponding argument for the feature screen, whose thresholds are fixed.
#' @return A list (the filtered state) with the following elements:
#'   \tabular{ll}{
#'     \code{expression} \cr \tab Numeric matrix with non-negative expression values. \cr
#'     \code{features}   \cr \tab A data.frame with feature meta-data for rows of expression. \cr
#'     \code{samples}    \cr \tab A data.frame with observation meta-data for columns of expression. \cr
#'   } 
#' @examples
#' set.seed(101)
#' exprs <- h0testr::sim1(n_obs=6, n_feats=12, mcar_p=0.75)$mat
#' feats <- data.frame(feature_id=rownames(exprs))
#' samps <- data.frame(observation_id=colnames(exprs))
#' state <- list(expression=exprs, features=feats, samples=samps)
#' config <- list(feat_col="feature_id", obs_col="observation_id")
#' state2 <- h0testr::prefilter(state, config)
#' print(state2)

prefilter <- function(state, config, n_features_min=2) {

  check_config(config)

  f.log_block("prefilter features and samples", config=config)
  f.msg("before filtering features:", config=config)
  f.report_state(state, config)

  ## f.prefilter_features() rather than filter_features(): the thresholds here are its
  ##   own fixed ones and not config$n_samples_min, which belongs to filter(), and there
  ##   is no design to screen against before normalization and aggregation:

  state <- f.prefilter_features(state)
  f.msg("after filtering features", config=config)
  f.report_state(state, config)
  
  state <- filter_observations(state, config, n_features_min=n_features_min)
  f.msg("after filtering samples", config=config)
  f.report_state(state, config)
  
  f.check_state(state, config)
  
  return(state)
}

## helper for add_filter_stats(): each statistic there is written into a metadata column
##   whose rows are identified by position, so the names the statistic carries have to
##   match the ids the metadata carries. `==` returns logical(0) when either side is NULL,
##   and all(logical(0)) is TRUE, so comparing the two directly passed vacuously in
##   exactly the cases where the alignment could not be verified: a config$feat_col or
##   config$obs_col that names no column, which includes the "" that new_config() ships,
##   and a state$expression with no dimnames. Refused rather than assumed:

f.check_stat_ids <- function(nms, ids, stat, key, config) {

  if(is.null(ids)) {
    f.err("add_filter_stats:", paste0("config$", key), "names no column of the",
      "metadata, so the alignment of", stat, "with it cannot be checked;", "\n",
      " ", paste0("config$", key, ":"), config[[key]], "\n",
      "  to fix, run h0testr::initialize(), which sets it, or set it yourself",
      config=config)
  }

  if(is.null(nms)) {
    f.err("add_filter_stats:", stat, "has no names, state$expression having no",
      "dimnames, so its alignment with", paste0("config$", key), "cannot be checked",
      config=config)
  }

  if(length(nms) != length(ids)) {
    f.err("add_filter_stats: length of", stat, ":", length(nms), "!= length of the",
      paste0("config$", key), "column:", length(ids), config=config)
  }

  if(!all(nms == ids)) {
    i <- which(nms != ids)
    f.err("add_filter_stats:", stat, "is not aligned with the column",
      paste0("config$", key), "names;", length(i), "of", length(nms), "differ,",
      "first at index", i[1], ":", nms[i[1]], "vs", ids[i[1]], config=config)
  }

  return(invisible(TRUE))
}

#' Add filter statistics
#' @description
#'   Adds filtering-related statistics to \code{state$features}, 
#'     and \code{state$samples}.
#' @details
#'   Wrapper for \code{samples_per_feature()}, \code{feature_median_expression()},
#'     \code{features_per_sample()}. Also reports quantiles of distributions.
#'   Each statistic is written into a metadata column by position, so
#'     \code{config$feat_col} and \code{config$obs_col} have to name the columns holding
#'     the ids that \code{rownames(state$expression)} and
#'     \code{colnames(state$expression)} carry, and the two have to agree. A key that
#'     names no column, which includes the \code{""} that \code{h0testr::new_config()}
#'     ships, a \code{state$expression} without dimnames, and a genuine mismatch are all
#'     errors. \code{h0testr::initialize()} sets both keys.
#'   See documentation for \code{h0testr::new_config()}
#'     for more detailed description of configuration parameters.
#' @param state List with elements formatted like the list returned by \code{read_data()}:
#'   \tabular{ll}{
#'     \code{expression} \cr \tab Numeric matrix with non-negative expression values. \cr
#'     \code{features}   \cr \tab A data.frame with feature meta-data for rows of expression. \cr
#'     \code{samples}    \cr \tab A data.frame with observation meta-data for columns of expression. \cr
#'   } 
#' @param config List with configuration values. Requires the following keys:
#'   \tabular{ll}{
#'     \code{feat_col}  \cr \tab Name of column (character) in \code{feature_file_in} that corresponds to rows of \code{data_file_in}. \cr
#'     \code{obs_col}   \cr \tab Name of column (character) in \code{sample_file_in} that corresponds to columns of \code{expression}. \cr
#'   }
#' @return A list (the processed state) with the following elements:
#'   \tabular{ll}{
#'     \code{expression} \cr \tab Numeric matrix with non-negative expression values. \cr
#'     \code{features}   \cr \tab A data.frame with feature meta-data for rows of expression. \cr
#'     \code{samples}    \cr \tab A data.frame with observation meta-data for columns of expression. \cr
#'   } 
#' @examples
#' exprs <- h0testr::sim1(n_obs=6, n_feats=8)$mat
#' feats <- data.frame(feature_id=rownames(exprs))
#' samps <- data.frame(observation_id=colnames(exprs))
#' state <- list(expression=exprs, features=feats, samples=samps)
#' config <- h0testr::new_config()      ## defaults
#' config$save_state <- FALSE             ## default is TRUE
#' config$feat_col <- config$feat_id_col
#' config$obs_col <- config$obs_id_col
#' state <- h0testr::add_filter_stats(state, config)
#' print(state)

add_filter_stats <- function(state, config) {
  
  check_config(config)
  
  n <- samples_per_feature(state, config)
  f.check_stat_ids(names(n), state$features[[config$feat_col]],
    "samples_per_feature()", "feat_col", config)
  state$features[[config$n_samples_expr_col]] <- n

  m <- feature_median_expression(state, config)
  f.check_stat_ids(names(m), state$features[[config$feat_col]],
    "feature_median_expression()", "feat_col", config)
  state$features[[config$median_raw_col]] <- m

  n <- features_per_sample(state, config)
  f.check_stat_ids(names(n), state$samples[[config$obs_col]],
    "features_per_sample()", "obs_col", config)
  state$samples[[config$n_features_expr_col]] <- n
  
  n <- apply(state$expression, 1, function(v) sum(!is.na(v)))
  f.msg("samples per feature:", config=config)
  f.quantile(n, config, digits=0)

  n <- apply(state$expression, 2, function(v) sum(!is.na(v)))
  f.msg("features per sample", config=config)
  f.quantile(n, config, digits=0)
  
  return(state)
}

#' Filter features and samples
#' @description
#'   Filter features and samples based on expression. 
#' @details 
#'   Filters out features measured in too few samples, and filters out samples
#'     with too few measured features. A value is missing if and only if it is
#'     \code{NA}; raw zeros are converted to \code{NA} by
#'     \code{h0testr::initialize()}.
#'   Features are filtered before observations, and then
#'     \code{filter_features_by_estimability()} runs last, since dropping an
#'     observation changes every feature's missingness pattern. This is a single
#'     pass: features dropped for lack of estimability are not fed back into
#'     \code{filter_observations()}.
#'   Features and/or samples considered constant if
#'     \code{length(unique(expression_values)) \%in\% 1}.
#'   Stops with an error if any feature is left with no measured value at all.
#'     Such a feature has no intensity for any downstream step to model, so an
#'     imputer would have to invent one; the filtering criteria above already
#'     remove them, and one surviving means a criterion was disabled or a
#'     threshold set too low.
#'   See documentation for \code{h0testr::new_config()}
#'     for more detailed description of configuration parameters.
#' @param state A list with elements like that returned by `read_data()`:
#'   \tabular{ll}{
#'     \code{expression} \cr \tab Numeric matrix with non-negative expression values. \cr
#'     \code{features}   \cr \tab A data.frame with feature meta-data for rows of expression. \cr
#'     \code{samples}    \cr \tab A data.frame with observation meta-data for columns of expression. \cr
#'   }
#' @param config List with configuration values. Uses the following keys:
#'   \tabular{ll}{
#'     \code{feat_col}            \cr \tab Name of column in \code{state$features} matching \code{rownames(state$expression)}.
#'     \code{obs_col}             \cr \tab Name of column in \code{state$samples} matching \code{colnames(state$expression)}.
#'     \code{n_features_min}      \cr \tab Minimum number (non-negative numeric) of features with a non-NA value in the observation to keep observation. \cr
#'     \code{n_samples_min}       \cr \tab Minimum number (non-negative numeric) of samples with a non-NA value for the feature to keep feature. \cr
#'     \code{median_raw_col}      \cr \tab Name (character) of new column in feature metadata to hold median expression in expressing samples. \cr
#'     \code{n_samples_expr_col}  \cr \tab Name (character) of new column in feature metadata to hold number of expressing samples. \cr
#'     \code{n_features_expr_col} \cr \tab Name (character) of new column in sample metadata to hold number of measured features. \cr
#'     \code{estimability}        \cr \tab Requirement placed on \code{config$test_term}; scalar character in \code{c("test", "term", "full")}. \cr
#'     \code{df_resid_min}        \cr \tab Minimum residual degrees of freedom (non-negative numeric) to keep feature. \cr
#'     \code{df_test_col}         \cr \tab Name (character) of new column in feature metadata to hold \code{df_test}. \cr
#'     \code{df_resid_col}        \cr \tab Name (character) of new column in feature metadata to hold \code{df_resid}. \cr
#'   }
#' @param remove_constant Logical scalar: if constant rows and columns of \code{state$expression} should be removed.
#' @param filter_by_formula Logical scalar: if \code{filter_features_by_formula()} should be run after other feature filters.
#' @param filter_by_estimability Logical scalar: if \code{filter_features_by_estimability()} should be run after features and observations have been filtered.
#' @return A list with elements like that returned by \code{read_data()}:
#'   \tabular{ll}{
#'     \code{expression} \cr \tab Numeric matrix with non-negative expression values. \cr
#'     \code{features}   \cr \tab A data.frame with feature meta-data for rows of expression. \cr
#'     \code{samples}    \cr \tab A data.frame with observation meta-data for columns of expression. \cr
#'   }
#' @examples
#' set.seed(101)
#' exprs <- h0testr::sim1(n_obs=8, n_feats=12, mcar_p=0.2)$mat
#' feats <- data.frame(feature_id=rownames(exprs))
#' samps <- data.frame(observation_id=colnames(exprs), age=c(rep("4m", 4), rep("12m", 4)))
#' state <- list(expression=exprs, features=feats, samples=samps)
#'
#' ## assume default median_raw_col, n_samples_expr_col, and n_features_expr_col are ok:
#' config <- h0testr::new_config()        ## defaults
#' config$save_state <- FALSE             ## default is TRUE
#' config$feat_col <- config$feat_id_col
#' config$obs_col <- config$obs_id_col
#' config$n_features_min <- 3
#' config$n_samples_min <- 2
#'
#' ## filter_features_by_estimability() needs test_term to fit frm:
#' config$frm <- ~age
#' config$test_term <- "age"
#' config$reference_levels <- c(age="4m")
#' out <- h0testr::filter(state, config)
#' print(out$state)
#' str(out$config)

filter <- function(state, config, remove_constant=TRUE, filter_by_formula=TRUE,
    filter_by_estimability=TRUE) {

  check_config(config)

  f.msg("filter: remove_constant:", remove_constant,
    "; filter_by_formula:", filter_by_formula,
    "; filter_by_estimability:", filter_by_estimability, config=config)

  state <- filter_features(state, config,
    remove_constant=remove_constant, filter_by_formula=filter_by_formula)

  state <- filter_observations(state, config,
    remove_constant=remove_constant)

  ## last, and after observations have been dropped: dropping an observation
  ##   changes every feature's missingness pattern, and so its ranks. Single
  ##   pass, so dropping features here does not re-qualify any observation:

  if(filter_by_estimability) {
    state <- filter_features_by_estimability(state, config)
  }

  state <- add_filter_stats(state, config)

  ## nothing downstream can work with a feature measured in no sample: it has no
  ##   intensity to model, so an imputer would have to invent one outright.
  ##   filter_features() and filter_features_by_estimability() both remove such
  ##   features, so one surviving to here means a criterion was disabled or a
  ##   threshold set too low. Reported plainly rather than left for an imputer to
  ##   trip over several steps later:

  n_obs <- rowSums(!is.na(state$expression))   ## a count, not a value

  if(any(n_obs %in% 0)) {

    i_none <- n_obs %in% 0
    nom <- state$features[[config$feat_col]][i_none]
    if(is.null(nom)) nom <- rownames(state$expression)[i_none]

    f.err("filter:", sum(i_none), "of", length(n_obs), "features have no",
      "measured values after filtering, so nothing downstream can model or",
      "impute them;", "\n",
      "  check the filtering criteria (remove_constant, filter_by_formula,",
      "filter_by_estimability) and their thresholds;", "\n",
      "  features:", paste(utils::head(nom, 5), collapse=", "), config=config)
  }

  f.check_state(state, config)
  f.report_state(state, config)
  
  prfx <- "filtered"
  if(!is.null(config$run_order)) {
    i <- config$run_order %in% "filter"
    if(any(i)) {
      prfx <- paste0(which(i)[1] + 2, ".filtered")
    }
  }
  f.save_state(state, config, prefix=prfx)
  
  return(list(state=state, config=config))
}
