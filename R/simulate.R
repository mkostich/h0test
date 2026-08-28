## helper; argument checking for simulators. Bounds are inclusive unless open=T; ninf_ok
##   admits -Inf, which mnar_c0 uses to turn MNAR off:

f.chk_num <- function(x, nom, lo=-Inf, hi=Inf, int=F, ninf_ok=F, open=F) {

  if(!is.numeric(x) || length(x) != 1 || is.na(x)) {
    stop(nom, " must be a scalar non-missing numeric; got: ",
      paste0(utils::head(as.character(x), 3), collapse=", "), " (length ", length(x), ")")
  }

  if(is.infinite(x)) {
    if(ninf_ok && x < 0) return(invisible(NULL))
    stop(nom, " must be finite; got: ", x)
  }

  if(int && x != round(x)) stop(nom, " must be a whole number; got: ", x)

  if(open) {
    if(x <= lo || x >= hi) stop(nom, " must lie in (", lo, ", ", hi, "); got: ", x)
  } else {
    if(x < lo || x > hi) stop(nom, " must lie in [", lo, ", ", hi, "]; got: ", x)
  }

  return(invisible(NULL))
}

## helper; log_m_mean and log_m_sd may each be 0, but not both: log(feature mean) redrawn
##   until strictly positive, and a mean of 0 with an sd of 0 has no positive value to reach:

f.chk_log_m <- function(log_m_mean, log_m_sd, nom) {

  if(log_m_mean <= 0 && log_m_sd <= 0) {
    stop(nom, ": log_m_mean and log_m_sd cannot both be 0, log(feature means) being drawn until ",
      "it is strictly positive, so a mean of 0 with an sd of 0 leaves no draw to make; ",
      "log_m_mean: ", log_m_mean, "; log_m_sd: ", log_m_sd)
  }

  return(invisible(NULL))
}

## helper; n normally distributed positive values with means m and sds s:

f.sim_rnorm_pos <- function(n, m, s) {

  if(length(m) != length(s)) stop("f.sim_rnorm_pos: length(m) != length(s)")

  ## retry below indexes m and s with length n logical; indexing shorter vector
  ##   past end gives NA:

  m <- rep_len(m, n)
  s <- rep_len(s, n)

  ## the parameter checks, finiteness first so that an NA mean is caught here rather than
  ##   reaching if(any(m < 0)) as an NA condition:

  f.bad <- function(i) {                 ## name the offenders, not the whole vector, since n
    i <- utils::head(which(i), 5)        ##   can be in the tens of thousands here:
    return(paste0(i, ":(m=", m[i], ", s=", s[i], ")", collapse=" "))
  }

  i_bad <- !is.finite(m) | !is.finite(s)
  if(any(i_bad)) {
    stop("f.sim_rnorm_pos: m and s must be finite; first offending positions: ", f.bad(i_bad))
  }

  i_bad <- m < 0
  if(any(i_bad)) {
    stop("f.sim_rnorm_pos: any(m < 0); first offending positions: ", f.bad(i_bad))
  }

  i_bad <- m <= 0 & s <= 0
  if(any(i_bad)) {
    stop("f.sim_rnorm_pos: no positive value is possible where m <= 0 and s <= 0; ",
      "first offending positions: ", f.bad(i_bad))
  }

  v <- stats::rnorm(n, mean=m, sd=s)
  i <- v <= 0
  i[is.na(i)] <- T

  while(any(i)) {
    v[i] <- stats::rnorm(sum(i), mean=m[i], sd=s[i])
    i <- v <= 0
    i[is.na(i)] <- T
  }

  return(v)
}

## helper; argument checking for per-term specification: either scalar, applied to every term
##   of frm, or vector named for terms it applies to:

f.chk_term_vec <- function(x, nom, terms, lo=-Inf, hi=Inf, int=F) {

  if(!is.numeric(x) || !length(x) || any(!is.finite(x))) {
    stop(nom, " must be non-missing finite numeric; got: ",
      paste0(utils::head(as.character(x), 3), collapse=", "), " (length ", length(x), ")")
  }

  if(is.null(names(x))) {

    if(length(x) != 1) {
      stop(nom, " must be a scalar or a vector named for the terms of frm; got an unnamed ",
        "vector of length ", length(x))
    }
    out <- rep(x, length(terms))
    names(out) <- terms

  } else {

    if(any(!nzchar(names(x)))) {
      stop(nom, " has an element with no name at position ", which(!nzchar(names(x)))[1])
    }
    if(any(duplicated(names(x)))) {
      stop(nom, " names a term more than once: ",
        paste0(unique(names(x)[duplicated(names(x))]), collapse=", "))
    }

    ## unmatched name is refused rather than ignored:

    bad <- setdiff(names(x), terms)
    if(length(bad)) {
      stop(nom, " names something that is not a term of frm: ", paste0(bad, collapse=", "),
        "; terms of frm: ", paste0(terms, collapse=", "))
    }

    out <- rep(0, length(terms))         ## term not named by caller gets no effect:
    names(out) <- terms
    out[names(x)] <- x
  }

  i_bad <- out < lo | out > hi
  if(any(i_bad)) {
    stop(nom, " must lie in [", lo, ", ", hi, "]; got: ",
      paste0(names(out)[i_bad], "=", out[i_bad], collapse=", "))
  }

  if(int) {
    i_bad <- out != round(out)
    if(any(i_bad)) {
      stop(nom, " must be whole numbers; got: ",
        paste0(names(out)[i_bad], "=", out[i_bad], collapse=", "))
    }
  }

  return(out)
}

## helper for sim1() and sim_design(); ensures values strictly positive (> 0):

f.sim0 <- function(n_obs, feat_means, feat_sds) {

  mat <- NULL

  for(i_obs in 1:n_obs) {
    if(is.matrix(feat_means)) {
      v <- f.sim_rnorm_pos(n=nrow(feat_means), m=feat_means[, i_obs], s=feat_sds[, i_obs])
    } else {
      v <- f.sim_rnorm_pos(n=length(feat_means), m=feat_means, s=feat_sds)
    }
    mat <- cbind(mat, v)
  }

  return(mat)
}

## technical replication; each replicate is drawn around its sample's value with CV cv_reps and
##   redrawn until strictly positive, like feature means are:

f.sim_tech_reps <- function(mat, reps_per_sample, cv_reps) {

  if(reps_per_sample < 2) return(mat)

  f <- function(v) {
    f0 <- function(val) f.sim_rnorm_pos(n=reps_per_sample, m=val, s=cv_reps * val)
    return(list(t(sapply(v, f0))))
  }
  tmp_list <- apply(mat, 2, f)

  for(nom in names(tmp_list)) {
    tmp_list[[nom]] <- tmp_list[[nom]][[1]]
    colnames(tmp_list[[nom]]) <- paste0(nom, "_rep", 1:reps_per_sample)
  }
  mat <- do.call(cbind, tmp_list)

  return(mat)
}

## helper for sim1() and sim_design(); estimate p(missing|log(intensity)) using
##   logit model. Drops cells in mat randomly based on p(missing|log(intensity)):

f.mnar <- function(mat, mnar_c0, mnar_c1, mnar_off=0.0001) {

  f <- function(v) {

    resp <- mnar_c0 + mnar_c1 * log(v + mnar_off)   ## logit(p_mnar) ~ c0 + c1 * log(intensity)

    ## stats::plogis() rather than exp(resp) / (1 + exp(resp)): 

    p_mnar <- stats::plogis(resp)                   ## inverse logit

    i_mnar <- as.logical(stats::rbinom(length(v), 1, p_mnar))
    v[i_mnar] <- NA

    return(v)
  }

  out <- apply(mat, 2, f)

  ## apply() returns a vector, not a one row matrix, when mat has a single row;  
  ##   restore shape and labels:

  if(!is.matrix(out)) {
    out <- matrix(out, nrow=nrow(mat), ncol=ncol(mat), dimnames=dimnames(mat))
  }

  return(out)
}

## mcar:

f.mcar <- function(mat, mcar_p) {

  i_mcar <- as.logical(stats::rbinom(length(c(mat)), 1, mcar_p))
  i_mcar <- matrix(i_mcar, nrow=nrow(mat), ncol=ncol(mat))
  mat[i_mcar] <- NA

  return(mat)
}

## Heterogenous number of peptides per gene if p_drop > 0. 
##   draw is per peptide; floor is per gene: 

f.pep_drop <- function(mat, peps_per_gene, p_drop, genes=NULL) {

  if(!(p_drop > 0 && peps_per_gene >= 2)) return(mat)

  if(is.null(genes)) genes <- sub("_pep[0-9]+$", "", rownames(mat))
  genes <- as.character(genes)

  if(length(genes) != nrow(mat)) {
    stop("f.pep_drop: length(genes) != nrow(mat); length(genes): ", length(genes),
      "; nrow(mat): ", nrow(mat))
  }

  i_drop <- as.logical(stats::rbinom(nrow(mat), 1, p_drop))

  ## only genes that lost something can lose everything, so loop is over those:

  for(gene in unique(genes[i_drop])) {
    i_gene <- which(genes %in% gene)
    if(all(i_drop[i_gene])) {
      i_drop[i_gene[sample(length(i_gene), 1)]] <- FALSE
    }
  }

  ## drop=F: a single surviving row would become a vector and lose its rownames,
  ##   which is where caller reads gene labels from:

  return(mat[!i_drop, , drop=F])
}

#' Simulate a one-condition dataset
#' @description
#'   Simulates non-negative data matrix with \code{n_feats} rows and \code{n_obs} columns.
#' @details
#'   Feature means: \code{log(feature_mean) ~ rnorm(mean=log_m_mean, sd=log_m_sd)}.
#'   That draw is redrawn until it is strictly positive, so really a normal truncated
#'     below 0 and \code{feat_mean} is never below 1. With the default
#'     \code{log_m_mean=11, log_m_sd=2.7} the truncation is 4 SDs out and so takes almost no
#'     draws, but a smaller \code{log_m_mean} or a larger \code{log_m_sd} will pull the realized
#'     mean of \code{log(feat_mean)} above \code{log_m_mean}.
#'   Dispersion of feature CVs: \code{log(feature_cv) ~ rnorm(mean=log_cv_mean, sd=log_cv_sd)}.
#'   MNAR: \code{logit(p(mnar|log(m))) ~ mnar_c0 + mnar_c1 * log(m + mnar_off)}.
#'   For no MNAR, set \code{mnar_c0=-Inf, mnar_c1=0}.
#'   For no MCAR, set \code{mcar_p=0}.
#' @param n_obs Scalar whole number of observations to simulate, with \code{n_obs >= 1}.
#' @param n_feats Scalar whole number of features to simulate, with \code{n_feats >= 1}.
#' @param log_m_mean Mean of \code{log(feature means)}. Scalar numeric,
#'   with \code{log_m_mean >= 0}.
#' @param log_m_sd Standard deviation of \code{log(feature means)} around
#'   \code{log_m_mean}. Scalar numeric, with \code{log_m_sd >= 0}. Either of
#'   \code{log_m_mean} and \code{log_m_sd} may be 0, but not both: the draw is redrawn until it
#'   is strictly positive, which a mean of 0 with an SD of 0 never is.
#' @param log_cv_mean Mean of \code{log(feature CVs)}, the CVs being of features around
#'   their respective feature means. Scalar numeric.
#' @param log_cv_sd Standard deviation of \code{log(feature CVs)} around
#'   \code{log_cv_mean}. Scalar numeric, with \code{log_cv_sd >= 0}.
#' @param mnar_c0 Intercept of logistic fit of \code{p(mnar|log(intensity))}. Scalar
#'   numeric; \code{-Inf} turns MNAR off.
#' @param mnar_c1 Slope of logistic fit of \code{p(mnar|log(intensity))}. Scalar numeric.
#' @param mnar_off Offset for taking logit of \code{p(mnar|log(intensity))}.
#'   Scalar numeric, with \code{0 < mnar_off < 1}.
#' @param mcar_p Probability of missing completely at random.
#'   Scalar numeric, with \code{0 <= mcar_p <= 1}.
#' @return List with elements:
#'   \tabular{ll}{
#'     \code{mat}       \cr \tab Numeric matrix with \code{n_feats} rows and \code{n_obs} columns,
#'       of positive whole numbers, with \code{NA} wherever MNAR or MCAR made a value
#'       missing. \cr
#'     \code{feat_mean} \cr \tab Numeric vector of strictly positive parameter means for each
#'       feature in \code{mat}. \cr
#'     \code{feat_cv}   \cr \tab Numeric vector of strictly positive parameter CVs for each
#'       feature in \code{mat}. \cr
#'   }
#' @seealso \code{\link{sim_design}}, for a matrix with a design imposed on it, of any number of
#'   groups, together with the annotation and configuration needed to run it.
#' @examples
#' ## default missing value settings:
#' rslt <- h0testr::sim1(n_obs=6, n_feats=8)
#' print(rslt$mat)
#' print(rslt$feat_mean)
#' print(rslt$feat_cv)
#'
#' ## no missing values:
#' rslt <- h0testr::sim1(n_obs=6, n_feats=8, mnar_c0=-Inf, mnar_c1=0, mcar_p=0)
#' print(rslt$mat)
#' print(rslt$feat_mean)
#' print(rslt$feat_cv)

sim1 <- function(n_obs, n_feats, log_m_mean=11, log_m_sd=2.7,
    log_cv_mean=-0.75, log_cv_sd=0.5, mnar_c0=4.65, mnar_c1=-0.5,
    mnar_off=0.0001, mcar_p=0.002) {

  f.chk_num(n_obs, "sim1: n_obs", lo=1, int=T)
  f.chk_num(n_feats, "sim1: n_feats", lo=1, int=T)
  f.chk_num(log_m_mean, "sim1: log_m_mean", lo=0)
  f.chk_num(log_m_sd, "sim1: log_m_sd", lo=0)
  f.chk_num(log_cv_mean, "sim1: log_cv_mean")
  f.chk_num(log_cv_sd, "sim1: log_cv_sd", lo=0)
  f.chk_num(mnar_c0, "sim1: mnar_c0", ninf_ok=T)
  f.chk_num(mnar_c1, "sim1: mnar_c1")
  f.chk_num(mnar_off, "sim1: mnar_off", lo=0, hi=1, open=T)
  f.chk_num(mcar_p, "sim1: mcar_p", lo=0, hi=1)
  f.chk_log_m(log_m_mean, log_m_sd, "sim1")

  ## feature mean, cv, and sd:
  m <- exp(f.sim_rnorm_pos(n=n_feats, m=log_m_mean, s=log_m_sd))
  cv <- exp(stats::rnorm(n_feats, mean=log_cv_mean, sd=log_cv_sd))
  s <- m * cv

  mat <- f.sim0(n_obs=n_obs, feat_means=m, feat_sds=s)

  rownames(mat) <- names(m) <- names(cv) <- paste0("feat_", 1:n_feats)
  colnames(mat) <- paste0("obs_", 1:n_obs)

  ## mnar:
  mat <- f.mnar(mat, mnar_c0=mnar_c0, mnar_c1=mnar_c1, mnar_off=mnar_off)

  ## mcar; f.mcar() is what these three lines used to spell out in place:
  mat <- f.mcar(mat, mcar_p=mcar_p)

  ## ceiling(), as in sim_design(): round() can take a strictly positive draw below 1 to 0, which
  ##   undoes what f.sim_rnorm_pos() went to the trouble of guaranteeing:

  return(list(mat=ceiling(mat), feat_mean=m, feat_cv=cv))
}

#' Build a samples table for sim_design()
#' @description
#'   Builds sample-level annotation \code{data.frame} for \code{\link{sim_design}} from a
#'   crossing of factor levels and a set of continuous covariates.
#' @details
#'   Rows are one biological sample each; \code{\link{sim_design}} expands technical replicates
#'   itself. Cells are generated by \code{expand.grid()}, which varies the first factor fastest,
#'   and that is the order \code{n_per_cell} is read in. The result is an ordinary
#'   \code{data.frame}, so anything not offered here can be done to it directly:
#'   \tabular{ll}{
#'     unbalance \cr \tab \code{samps <- samps[-c(3, 7), ]} \cr
#'     shuffle   \cr \tab \code{samps <- samps[sample(nrow(samps)), ]} \cr
#'     confound  \cr \tab \code{samps$age <- samps$age + 10 * (samps$sex == "M")} \cr
#'   }
#'   The default covariate draw is standard normal rather than uniform because a uniform
#'   covariate has no tails and so never presents a high-leverage observation.
#' @param factors Named list of level vectors, one element per factor variable; the first level
#'   of each is the reference level. \code{NULL} for no factors, in which case \code{n} is
#'   required.
#' @param covariates Named list of continuous covariate specifications, one element per
#'   covariate; each being either a length 2 numeric range, drawn uniformly, or a function of
#'   \code{n} returning \code{n} finite numeric values, or \code{NULL} for the default draw from
#'   \code{stats::rnorm(n, mean=0, sd=1)}. A character vector of names is shorthand for all
#'   defaults. \code{NULL} for no covariates.
#' @param n_per_cell Number of samples per cell of the factor crossing. Scalar whole number, or
#'   one whole number per cell in \code{expand.grid()} order, each \code{>= 0}; a cell given 0 is
#'   left out entirely, which is one way to reach an empty cell. Ignored when \code{factors} is
#'   \code{NULL}.
#' @param n Scalar whole number of samples, with \code{n >= 1}. Required when \code{factors} is
#'   \code{NULL}, and refused otherwise, \code{n_per_cell} determining the count in that case.
#' @return A \code{data.frame} with one row per sample and columns:
#'   \tabular{ll}{
#'     \code{sample_id} \cr \tab Character sample identifier, \code{samp1} onwards. \cr
#'     one per factor \cr \tab Factor, with levels in the order given, reference level first. \cr
#'     one per covariate \cr \tab Numeric, on whatever scale it was drawn. \cr
#'   }
#' @seealso \code{\link{sim_design}}, which consumes this.
#' @examples
#' ## balanced crossing, three samples per cell:
#' samps <- h0testr::sim_samples(factors=list(sex=c("F", "M"), geno=c("WT", "KO", "HET")))
#' print(samps)
#'
#' ## unbalanced, with a covariate drawn uniformly on a range:
#' samps <- h0testr::sim_samples(factors=list(sex=c("F", "M")),
#'   covariates=list(age=c(20, 60)), n_per_cell=c(5, 2))
#' print(samps)
#'
#' ## one cell left empty, and two covariates taking the default standard normal draw:
#' samps <- h0testr::sim_samples(factors=list(sex=c("F", "M"), geno=c("WT", "KO")),
#'   covariates=c("age", "intake"), n_per_cell=c(3, 3, 3, 0))
#' print(samps)
#'
#' ## no factors, so n is required; the covariate is drawn by a function of n:
#' samps <- h0testr::sim_samples(covariates=list(age=function(n) stats::rnorm(n, 50, 10)), n=8)
#' print(samps)

sim_samples <- function(factors=NULL, covariates=NULL, n_per_cell=3, n=NULL) {

  if(is.null(factors) && is.null(covariates)) {
    stop("sim_samples: give factors, covariates, or both; both were NULL")
  }

  if(!is.null(factors)) {

    if(!is.null(n)) {
      stop("sim_samples: give n_per_cell with factors, or n without them, not both; got n: ", n)
    }
    if(!is.list(factors) || !length(factors) || is.null(names(factors)) ||
        any(!nzchar(names(factors))) || any(duplicated(names(factors)))) {
      stop("sim_samples: factors must be a non-empty list with a distinct name per factor; ",
        "got names: ", paste0(names(factors), collapse=", "))
    }

    for(nom in names(factors)) {
      lev <- factors[[nom]]
      if(!length(lev) || any(is.na(lev)) || any(duplicated(lev))) {
        stop("sim_samples: factors$", nom, " must be distinct non-missing levels; got: ",
          paste0(utils::head(as.character(lev), 5), collapse=", "))
      }
    }

    ## expand.grid() varies first factor fastest, which is the order n_per_cell is read in:

    cells <- expand.grid(factors, stringsAsFactors=F, KEEP.OUT.ATTRS=F)

    if(!is.numeric(n_per_cell) || !length(n_per_cell) || any(!is.finite(n_per_cell)) ||
        any(n_per_cell < 0) || any(n_per_cell != round(n_per_cell))) {
      stop("sim_samples: n_per_cell must be whole numbers >= 0; got: ",
        paste0(utils::head(as.character(n_per_cell), 5), collapse=", "))
    }
    if(!(length(n_per_cell) %in% c(1, nrow(cells)))) {
      stop("sim_samples: n_per_cell must be length 1 or one per cell; got length ",
        length(n_per_cell), " for ", nrow(cells), " cells")
    }
    n_per_cell <- rep_len(n_per_cell, nrow(cells))
    if(sum(n_per_cell) < 1) {
      stop("sim_samples: sum(n_per_cell) is 0, so there would be no samples")
    }

    ## cell given 0 is left out entirely; one way to achieve an empty cell:

    samps <- cells[rep(1:nrow(cells), times=n_per_cell), , drop=F]

    for(nom in names(factors)) {
      samps[[nom]] <- factor(samps[[nom]], levels=factors[[nom]])
    }

  } else {

    if(is.null(n)) stop("sim_samples: n is required when no factors are given")
    f.chk_num(n, "sim_samples: n", lo=1, int=T)
    samps <- data.frame(row.names=1:n)
  }

  n_rows <- nrow(samps)

  if(!is.null(covariates)) {

    if(is.character(covariates)) {        ## bare names: every one takes the default draw
      covariates <- stats::setNames(vector("list", length(covariates)), covariates)
    }
    if(!is.list(covariates) || !length(covariates) || is.null(names(covariates)) ||
        any(!nzchar(names(covariates))) || any(duplicated(names(covariates)))) {
      stop("sim_samples: covariates must be a list with a distinct name per covariate, or a ",
        "character vector of names; got names: ", paste0(names(covariates), collapse=", "))
    }

    bad <- intersect(names(covariates), c("sample_id", names(samps)))
    if(length(bad)) {
      stop("sim_samples: covariate name is already in use: ", paste0(bad, collapse=", "))
    }

    for(nom in names(covariates)) {

      spec <- covariates[[nom]]

      if(is.null(spec)) {
        v <- stats::rnorm(n_rows, mean=0, sd=1)
      } else if(is.function(spec)) {
        v <- spec(n_rows)
        if(!is.numeric(v) || length(v) != n_rows || any(!is.finite(v))) {
          stop("sim_samples: covariates$", nom, " must return ", n_rows, " finite numeric ",
            "values; got ", length(v), " of class: ", paste0(class(v), collapse=", "))
        }
      } else if(is.numeric(spec) && length(spec) == 2 && all(is.finite(spec)) &&
          spec[1] < spec[2]) {
        v <- stats::runif(n_rows, min=spec[1], max=spec[2])
      } else {
        stop("sim_samples: covariates$", nom, " must be NULL, a function of n, or an ",
          "increasing length 2 numeric range; got: ",
          paste0(utils::head(as.character(spec), 3), collapse=", "))
      }

      samps[[nom]] <- v
    }
  }

  samps <- cbind(data.frame(sample_id=paste0("samp", 1:n_rows), stringsAsFactors=F), samps)
  rownames(samps) <- NULL

  return(samps)
}

#' Simulate a dataset for an arbitrary design
#' @description
#'   Simulates a non-negative data matrix, with the annotation and configuration needed to run
#'   it, for the design given by \code{frm} and \code{samps}: any number of groups, continuous
#'   covariates, adjustment covariates, and interactions.
#' @details
#'   The mean structure is
#'   \code{log2(mean[feat, obs]) = log2(feat_mean[feat]) + x[obs, ] \%*\% beta[feat, ]}, with
#'   \code{x} the \code{stats::model.matrix()} of \code{frm}, so any number of groups, a
#'   continuous covariate, an adjustment covariate and an interaction are all the same here.
#'   Around the mean, \code{log(feature_mean) ~ rnorm(mean=log_m_mean, sd=log_m_sd)} and
#'   \code{log(feature_cv) ~ rnorm(mean=log_cv_mean, sd=log_cv_sd)}; values are drawn on the raw
#'   scale with an SD proportional to the mean, and then come technical replication, MNAR
#'   (\code{logit(p(mnar|log(m))) ~ mnar_c0 + mnar_c1 * log(m + mnar_off)}, turned off by
#'   \code{mnar_c0=-Inf, mnar_c1=0}), MCAR (turned off by \code{mcar_p=0}) and feature dropout.
#'   The defaults of those layers were fit to a large real dataset. See \code{\link{sim1}} for
#'   the same intensity, noise and missingness model with no design at all.
#'
#'   An effect size is a log2 fold change, but of different things for the two kinds of term. A
#'   factor dummy stays 0/1, so its coefficient is the fold change between the level the effect
#'   was planted on and the reference level. A continuous variable is centered and scaled to its
#'   realized mean and SD for the model matrix only, the returned samples table keeping the raw
#'   values, so its coefficient is the fold change per SD of the variable. One number therefore
#'   means a comparable amount of signal either way, without meaning the identical thing. The
#'   centering matters in itself too: without it a covariate with a nonzero mean would shift
#'   overall intensity off \code{log_m_mean} and silently change the global missingness rate.
#'
#'   Each significant gene takes its effect on ONE randomly chosen column of its term, with an
#'   equal chance of an increase or a decrease, so a factor of more than two levels is changed in
#'   one level and left alone in the others. Only columns that vary across samples are eligible,
#'   since a constant column has no contrast to express an effect with: the dummy of a factor
#'   level that is declared but has no samples is all zero, and planting there would give a
#'   \code{truth} the data does not carry. A term whose columns all turn out constant gets a
#'   warning and no effect. All features of a gene share their gene's coefficients. Per-term counts
#'   of significant genes are useful: a realistic confounder moves many genes while the term
#'   of interest moves few, as in \code{n_genes_signif=c(sex=2, age=40)}.
#'
#'   The feature-mean draw is redrawn until strictly
#'   positive, so \code{log(feat_mean)} is a normal truncated below 0 and \code{feat_mean} is
#'   never below 1. And \code{p_drop=0.75} is aggressive: with \code{peps_per_gene=3} a gene has
#'   all three of its peptides dropped about 42\% of the time, and then keeps one of them, chosen
#'   at random, since the floor is per gene. So every gene of \code{truth} still has a feature in
#'   \code{expression}, and the peptide counts pile up at 1 rather than spreading; lower
#'   \code{p_drop} or raise \code{peps_per_gene} for counts spread more evenly.
#'
#'   A design that cannot be fitted is not refused. Only \code{X \%*\% beta} is needed here,
#'   which is defined even when the design is rank deficient or has empty cells.
#' @param samps Sample-level annotation \code{data.frame}, one row per biological sample, holding
#'   every variable named in \code{frm}; see \code{\link{sim_samples}}, which builds one. A
#'   \code{sample_id} column is used if present and generated otherwise. Factor variables should
#'   be factors, with first level being the reference; character and logical variables
#'   are levelled in sorted order. Treatment contrasts are what make a factor coefficient a fold
#'   change against that reference level, so an ordered factor, a variable carrying a
#'   \code{contrasts} attribute, and an \code{options("contrasts")} whose unordered entry is not
#'   \code{contr.treatment} are all refused rather than silently changing what \code{truth} and
#'   \code{config$reference_levels} mean.
#' @param frm Formula giving the design, e.g. \code{~sex + age} or \code{~geno * treatment}. Its
#'   variables must all be columns of \code{samps}, and it must have at least one term.
#' @param test_term Term (scalar character) of \code{frm} to be tested, e.g. \code{"sex"}, passed
#'   through to \code{config$test_term}. Must be one of \code{frm}'s terms.
#' @param n_genes Scalar whole number of gene groups to simulate, with \code{n_genes >= 2}.
#' @param n_genes_signif Number of genes given a nonzero effect. Either a scalar whole number,
#'   applied to every term of \code{frm}, or a vector of whole numbers named for the terms it
#'   applies to, a term not named getting 0; each with
#'   \code{0 <= n_genes_signif <= n_genes}. Genes are chosen independently per term, so one gene
#'   can carry effects of several terms.
#' @param effects Effect size, as a log2 fold change: between levels for a factor term, and per
#'   SD of the variable for a continuous one; see Details. Either a scalar
#'   numeric, applied to every term of \code{frm}, or a numeric vector named for the terms it
#'   applies to, a term not named getting 0. A term that \code{n_genes_signif} gives genes to and
#'   this gives 0 to gets no effect and a warning, since naming one term here zeroes the others.
#' @param peps_per_gene Scalar whole number of features (peptides or precursors) per gene before
#'   dropout, with \code{peps_per_gene >= 1}.
#' @param reps_per_sample Scalar whole number of technical replicate observations per sample,
#'   with \code{reps_per_sample >= 1}.
#' @param cv_reps CV of technical replicates around their sample value. Scalar numeric, with
#'   \code{cv_reps >= 0}. Replicate draws are redrawn until strictly positive, as the feature
#'   means are, so a \code{cv_reps} near or above 1 pulls a sample's replicates above that
#'   sample's own value.
#' @param log_m_mean Mean of \code{log(feature means)}. Scalar numeric, with
#'   \code{log_m_mean >= 0}.
#' @param log_m_sd Standard deviation of \code{log(feature means)} around \code{log_m_mean}.
#'   Scalar numeric, with \code{log_m_sd >= 0}. Either of \code{log_m_mean} and \code{log_m_sd}
#'   may be 0, but not both.
#' @param log_cv_mean Mean of \code{log(feature CVs)}, the CVs being of features around their
#'   respective feature means. Scalar numeric.
#' @param log_cv_sd Standard deviation of \code{log(feature CVs)} around \code{log_cv_mean}.
#'   Scalar numeric, with \code{log_cv_sd >= 0}.
#' @param p_drop Probability of dropping each feature, giving a heterogeneous number of features
#'   per gene. Scalar numeric, with \code{0 <= p_drop <= 1}. Has no effect unless
#'   \code{peps_per_gene >= 2}. Every gene keeps at least one feature, so \code{p_drop=1} leaves
#'   one feature per gene and \code{truth} keeps every gene whatever \code{p_drop} says.
#' @param mnar_c0 Intercept of logistic fit of \code{p(mnar|log(intensity))}. Scalar numeric;
#'   \code{-Inf} turns MNAR off.
#' @param mnar_c1 Slope of logistic fit of \code{p(mnar|log(intensity))}. Scalar numeric.
#' @param mnar_off Offset for taking logit of \code{p(mnar|log(intensity))}. Scalar numeric, with
#'   \code{0 < mnar_off < 1}.
#' @param mcar_p Probability of missing completely at random. Scalar numeric, with
#'   \code{0 <= mcar_p <= 1}.
#' @return List with elements:
#'   \tabular{ll}{
#'     \code{state} \cr \tab List of \code{expression} (numeric matrix of positive whole numbers,
#'       features by observations, with \code{NA} wherever MNAR or MCAR made a value missing),
#'       \code{features} and \code{samples}, ready for \code{\link{init_state}}. \cr
#'     \code{config} \cr \tab \code{\link{new_config}} with \code{frm}, \code{test_term}, the
#'       four id columns and \code{reference_levels} set to match; \code{save_state} set
#'       \code{FALSE}, so that a simulation writes no files; and \code{n_features_min} set to 1,
#'       its default of 1000 being meant for a real dataset and otherwise filtering away every
#'       observation of a simulation of fewer than 1000 genes. \code{covariate_types} and
#'       \code{factor_levels} are left for \code{\link{init_state}} to derive. \cr
#'     \code{truth} \cr \tab Numeric matrix of true coefficients on the log2 scale, one row per
#'       surviving gene and one column per column of the model matrix, the intercept excluded. \cr
#'     \code{feat_gene} \cr \tab Character vector naming the gene of each feature in
#'       \code{expression}, named by feature. \cr
#'     \code{feat_mean} \cr \tab Numeric vector of reference-level parameter means for each
#'       feature. \cr
#'     \code{feat_cv} \cr \tab Numeric vector of within-group parameter CVs for each feature. \cr
#'     \code{x} \cr \tab The model matrix \code{truth} multiplies: one row per sample and the
#'       same columns as \code{truth}, with continuous variables centered and scaled. \cr
#'   }
#' @seealso \code{\link{sim_samples}}, which builds \code{samps}; \code{\link{sim1}} for a
#'   design-free matrix.
#' @examples
#' ## one factor and one covariate, with the covariate moving more genes, as a confounder would:
#' samps <- h0testr::sim_samples(factors=list(sex=c("F", "M")), covariates=list(age=c(20, 60)),
#'   n_per_cell=3)
#' sim <- h0testr::sim_design(samps, frm=~sex + age, test_term="sex", n_genes=8,
#'   n_genes_signif=c(sex=2, age=4), mnar_c0=-Inf, mnar_c1=0, mcar_p=0)
#' print(sim$state$expression)
#' print(sim$truth)
#' print(sim$state$samples)
#'
#' ## three groups, peptides and technical replicates:
#' samps <- h0testr::sim_samples(factors=list(geno=c("WT", "KO", "HET")), n_per_cell=2)
#' sim <- h0testr::sim_design(samps, frm=~geno, test_term="geno", n_genes=6, n_genes_signif=3,
#'   peps_per_gene=3, reps_per_sample=2)
#' print(sim$truth)
#' print(sim$feat_gene)
#'
#' ## no factors at all, so a continuous term is the one tested:
#' samps <- h0testr::sim_samples(covariates=list(age=c(20, 60)), n=10)
#' sim <- h0testr::sim_design(samps, frm=~age, test_term="age", n_genes=6, n_genes_signif=2)
#' print(sim$truth)

sim_design <- function(samps, frm, test_term, n_genes, n_genes_signif=0, effects=0.5,
    peps_per_gene=1, reps_per_sample=1, cv_reps=0.1,
    log_m_mean=11, log_m_sd=2.7, log_cv_mean=-0.75, log_cv_sd=0.5,
    p_drop=0.75, mnar_c0=4.65, mnar_c1=-0.5, mnar_off=0.0001, mcar_p=0.002) {

  ## the design:

  if(!is.data.frame(samps)) {
    stop("sim_design: samps must be a data.frame; got class: ",
      paste0(class(samps), collapse=", "))
  }
  if(nrow(samps) < 1 || ncol(samps) < 1) {
    stop("sim_design: samps must have at least one row and one column; got ", nrow(samps),
      " rows and ", ncol(samps), " columns")
  }
  if(!inherits(frm, "formula")) {
    stop("sim_design: frm must be a formula; got class: ", paste0(class(frm), collapse=", "))
  }
  if(!is.character(test_term) || length(test_term) != 1 || is.na(test_term) ||
      !nzchar(test_term)) {
    stop("sim_design: test_term must be a scalar non-empty character; got: ",
      paste0(utils::head(as.character(test_term), 3), collapse=", "))
  }

  term_labs <- attr(stats::terms(frm), "term.labels")
  if(!length(term_labs)) {
    stop("sim_design: frm must have at least one term; got: ", paste0(deparse(frm), collapse=""))
  }
  if(!(test_term %in% term_labs)) {
    stop("sim_design: test_term is not a term of frm; test_term: ", test_term,
      "; terms of frm: ", paste0(term_labs, collapse=", "))
  }

  vars <- all.vars(frm)
  bad <- setdiff(vars, names(samps))
  if(length(bad)) {
    stop("sim_design: frm names variables missing from samps: ", paste0(bad, collapse=", "),
      "; samps has: ", paste0(names(samps), collapse=", "))
  }

  f.chk_num(n_genes, "sim_design: n_genes", lo=2, int=T)
  n_genes_signif <- f.chk_term_vec(n_genes_signif, "sim_design: n_genes_signif", term_labs,
    lo=0, hi=n_genes, int=T)
  effects <- f.chk_term_vec(effects, "sim_design: effects", term_labs)
  f.chk_num(peps_per_gene, "sim_design: peps_per_gene", lo=1, int=T)
  f.chk_num(reps_per_sample, "sim_design: reps_per_sample", lo=1, int=T)
  f.chk_num(cv_reps, "sim_design: cv_reps", lo=0)
  f.chk_num(log_m_mean, "sim_design: log_m_mean", lo=0)
  f.chk_num(log_m_sd, "sim_design: log_m_sd", lo=0)
  f.chk_num(log_cv_mean, "sim_design: log_cv_mean")
  f.chk_num(log_cv_sd, "sim_design: log_cv_sd", lo=0)
  f.chk_num(p_drop, "sim_design: p_drop", lo=0, hi=1)
  f.chk_num(mnar_c0, "sim_design: mnar_c0", ninf_ok=T)
  f.chk_num(mnar_c1, "sim_design: mnar_c1")
  f.chk_num(mnar_off, "sim_design: mnar_off", lo=0, hi=1, open=T)
  f.chk_num(mcar_p, "sim_design: mcar_p", lo=0, hi=1)
  f.chk_log_m(log_m_mean, log_m_sd, "sim_design")

  ## variable types, reference levels, and the scaling of continuous variables:

  ref_levels <- character(0)
  samps_z <- samps

  for(nom in vars) {

    v <- samps[[nom]]

    if(is.factor(v) || is.character(v) || is.logical(v)) {

      if(any(is.na(v))) {
        stop("sim_design: variable ", nom, " must be non-missing; got ", sum(is.na(v)),
          " missing of ", length(v))
      }
      if(!is.factor(v)) v <- factor(v)    ## an existing factor keeps the level order it has
      if(nlevels(v) < 2) {
        stop("sim_design: factor variable ", nom, " needs at least 2 levels; got: ",
          paste0(levels(v), collapse=", "))
      }

      ## factor coefficient is log2 fold change vs. reference level only under treatment contrasts:

      if(is.ordered(v)) {
        stop("sim_design: variable ", nom, " is an ordered factor, which stats::model.matrix() ",
          "gives polynomial contrasts, so its coefficients would be polynomial trends rather ",
          "than log2 fold changes against the reference level named in ",
          "config$reference_levels; declare it unordered, with factor(levels=c(",
          paste0("\"", levels(v), "\"", collapse=", "), "))")
      }
      if(!is.null(attr(v, "contrasts"))) {
        stop("sim_design: variable ", nom, " carries a contrasts attribute, which ",
          "stats::model.matrix() honors in place of treatment contrasts, so its coefficients ",
          "would not be log2 fold changes against the reference level named in ",
          "config$reference_levels; drop it with attr(samps$", nom, ", \"contrasts\") <- NULL")
      }

      samps[[nom]] <- samps_z[[nom]] <- v
      ref_levels[nom] <- levels(v)[1]

    } else if(is.numeric(v)) {

      if(any(!is.finite(v))) {
        stop("sim_design: numeric variable ", nom, " must be finite; got ",
          sum(!is.finite(v)), " non-finite of ", length(v))
      }
      s <- stats::sd(v)
      if(!(s > 0)) {
        stop("sim_design: numeric variable ", nom, " has no variance, so it cannot carry a ",
          "per-SD effect; every value is: ", v[1])
      }
      samps_z[[nom]] <- (v - mean(v)) / s

    } else {

      stop("sim_design: variable ", nom, " must be numeric, logical, character or factor; ",
        "got class: ", paste0(class(v), collapse=", "))
    }
  }

  ## the other way a factor gets contrasts that are not treatment contrasts:

  if(length(ref_levels)) {

    contr <- getOption("contrasts")
    contr_1 <- if(is.character(contr) && length(contr)) contr[[1]] else
      paste0("a ", paste0(class(contr), collapse="/"), " of length ", length(contr))

    if(!identical(contr_1, "contr.treatment")) {
      stop("sim_design: options(\"contrasts\") gives ", contr_1, " for unordered factors rather ",
        "than contr.treatment, so a factor coefficient would not be a log2 fold change against ",
        "the reference level named in config$reference_levels; restore the default with ",
        "options(contrasts=c(\"contr.treatment\", \"contr.poly\"))")
    }
  }

  ## sample ids, made up if the caller's samps does not carry them:

  if(!("sample_id" %in% names(samps))) {
    samps[["sample_id"]] <- paste0("samp", 1:nrow(samps))
  }
  samps[["sample_id"]] <- as.character(samps[["sample_id"]])
  if(any(is.na(samps[["sample_id"]]))) {
    stop("sim_design: samps$sample_id must be non-missing; got ",
      sum(is.na(samps[["sample_id"]])), " missing of ", nrow(samps))
  }
  if(any(duplicated(samps[["sample_id"]]))) {
    stop("sim_design: samps$sample_id must be unique; first repeated: ",
      utils::head(samps[["sample_id"]][duplicated(samps[["sample_id"]])], 1))
  }

  ## model matrix, and columns an effect can be planted on, which intercept is not:

  x_all <- stats::model.matrix(frm, data=samps_z)
  i_int <- colnames(x_all) == "(Intercept)"
  x <- x_all[, !i_int, drop=F]
  assign <- attr(x_all, "assign")[!i_int]

  if(!ncol(x)) {
    stop("sim_design: frm gives no model matrix columns besides the intercept; got: ",
      paste0(deparse(frm), collapse=""))
  }
  if(nrow(x) != nrow(samps)) {
    stop("sim_design: model.matrix() returned ", nrow(x), " rows for ", nrow(samps),
      " samples, so frm's variables are not all usable")
  }
  rownames(x) <- samps[["sample_id"]]

  ## column that takes same value in every sample cannot result in effect: 

  i_const <- !apply(x, 2, function(v) any(v != v[1]))

  ## the true coefficients, log2 scale, one row per gene and one column per model matrix column:

  genes <- paste0("gene", 1:n_genes)
  beta <- matrix(0, nrow=n_genes, ncol=ncol(x), dimnames=list(genes, colnames(x)))

  for(i_term in seq_along(term_labs)) {

    nom <- term_labs[i_term]
    i_cols <- which(assign == i_term & !i_const)

    ## term given genes but no effect size gets nothing:

    if(n_genes_signif[nom] && !effects[nom]) {
      warning("sim_design: no effect planted on term ", nom, ", because effects is 0 there ",
        "while n_genes_signif asks for ", n_genes_signif[nom], " gene(s); a term effects does ",
        "not name gets 0. effects: ", paste0(names(effects), "=", effects, collapse=", "))
    }

    if(!n_genes_signif[nom] || !effects[nom]) next

    if(!length(i_cols)) {
      warning("sim_design: no effect planted on term ", nom,
        ", because none of its model matrix columns vary across samples: ",
        paste0(colnames(x)[assign == i_term], collapse=", "))
      next
    }

    for(i_gene in sample(n_genes, n_genes_signif[nom])) {
      i_col <- i_cols[sample(length(i_cols), 1)]     ## sample(i_cols, 1) misreads a length 1
      sgn <- ifelse(as.logical(stats::rbinom(1, 1, prob=0.5)), 1, -1)
      beta[i_gene, i_col] <- sgn * effects[nom]
    }
  }

  ## features, and the gene each belongs to; expand.grid() varies peptide fastest:

  if(peps_per_gene >= 2) {
    feats <- do.call(paste0,
      expand.grid(paste0("_pep", 1:peps_per_gene), genes)[, c(2, 1), drop=F])
    feat_gene <- rep(genes, each=peps_per_gene)
  } else {
    feats <- genes
    feat_gene <- genes
  }
  names(feat_gene) <- feats

  ## Feature means and CVs, then per-observation means design implies. 
  ##   all features of a gene share its coefficients:

  n_feats <- length(feats)
  cv <- exp(stats::rnorm(n_feats, mean=log_cv_mean, sd=log_cv_sd))
  m <- exp(f.sim_rnorm_pos(n=n_feats, m=log_m_mean, s=log_m_sd))
  names(m) <- names(cv) <- feats

  mat_m <- m * 2^(beta[feat_gene, , drop=F] %*% t(x))
  mat_s <- mat_m * cv
  dimnames(mat_m) <- dimnames(mat_s) <- list(feats, samps[["sample_id"]])

  mat <- f.sim0(n_obs=nrow(samps), feat_means=mat_m, feat_sds=mat_s)
  dimnames(mat) <- list(feats, samps[["sample_id"]])

  ## technical replication, which is also where the observation ids come from:

  mat <- f.sim_tech_reps(mat, reps_per_sample=reps_per_sample, cv_reps=cv_reps)
  if(any(is.na(c(mat)) | c(mat) <= 0, na.rm=T)) {

    ## min() without na.rm=T reports NA:

    n_na <- sum(is.na(c(mat)))
    stop("sim_design: any(is.na(c(mat)) | c(mat) <= 0); NA: ", n_na, " of ", length(c(mat)),
      "; min of the rest: ", if(n_na < length(c(mat))) min(c(mat), na.rm=T) else "none left")
  }

  obs <- samps[rep(1:nrow(samps), each=reps_per_sample), , drop=F]
  obs[["observation_id"]] <- colnames(mat)
  obs <- obs[, c("observation_id", "sample_id",
    setdiff(names(obs), c("observation_id", "sample_id"))), drop=F]
  rownames(obs) <- NULL

  ## mnar -> mcar -> heterogenous number of features per gene:

  mat <- f.mnar(mat, mnar_c0=mnar_c0, mnar_c1=mnar_c1, mnar_off=mnar_off)
  mat <- f.mcar(mat, mcar_p=mcar_p)
  mat <- f.pep_drop(mat, peps_per_gene=peps_per_gene, p_drop=p_drop,
    genes=feat_gene[rownames(mat)])

  ## trim the parameters to the features that survived dropout:

  feat_gene <- feat_gene[rownames(mat)]
  m <- m[rownames(mat)]
  cv <- cv[rownames(mat)]
  beta <- beta[rownames(beta) %in% feat_gene, , drop=F]

  feats_df <- data.frame(feature_id=rownames(mat), gene_id=as.character(feat_gene),
    stringsAsFactors=F)

  ## configuration that matches what was just made; only the fields a user would set, 
  ##   since covariate_types and factor_levels are init_state()'s to derive:

  config <- new_config()
  config$frm <- frm
  config$test_term <- test_term
  config$feat_id_col <- "feature_id"
  config$gene_id_col <- "gene_id"
  config$obs_id_col <- "observation_id"
  config$sample_id_col <- "sample_id"
  config$reference_levels <- ref_levels
  config$save_state <- FALSE
  config$n_features_min <- 1

  state <- list(expression=ceiling(mat), features=feats_df, samples=obs)

  return(list(state=state, config=config, truth=beta, feat_gene=feat_gene, feat_mean=m,
    feat_cv=cv, x=x))
}
