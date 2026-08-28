## The effect size reported in the logfc column of the standardized table, from the
##   coefficients of the design matrix columns:

f.logfc_effect <- function(coefs, design, config) {

  x_test <- design$X[, design$cols_test, drop=F]

  if(!is.matrix(coefs) || ncol(coefs) != ncol(x_test)) {
    f.err("f.logfc_effect: coefs does not match the columns under test;",
      "columns under test:", ncol(x_test), "; columns of coefs:",
      if(is.matrix(coefs)) ncol(coefs) else paste("not a matrix:", class(coefs)),
      config=config)
  }

  ## config$contrast has a signed effect size:

  if(!is.null(design$contrast)) {
    out <- drop(coefs %*% design$contrast[design$cols_test])
    names(out) <- rownames(coefs)
    return(out)
  }

  if(ncol(x_test) %in% 1) {
    out <- coefs[, 1, drop=T]
    names(out) <- rownames(coefs)
    return(out)
  }

  ## contribution of tested terms to fitted value of every observation; 
  ##   one column per feature:

  contrib <- x_test %*% t(coefs)
  out <- apply(contrib, 2, function(v) diff(range(v)))
  names(out) <- rownames(coefs)

  return(out)
}

## Coefficients of design matrix columns, one row per feature and
##   one column per column of design$cols_test:

f.test_coefs <- function(result, method, design, config) {

  cols <- colnames(design$X)[design$cols_test]

  ## limma moderates variance and not coefficients, so fit$coefficients holds
  ##   least squares estimates for design:

  if(method %in% c("trend", "voom", "deqms")) {
    coefs <- result$fit$coefficients
    if(!is.matrix(coefs) || !all(cols %in% colnames(coefs))) return(NULL)
    return(coefs[, cols, drop=F])
  }

  ## test_lm() already reports coefficients under test, one column of its hit
  ##   table per:

  if(method %in% "lm") {
    nom <- ifelse(cols %in% "(Intercept)", "Intercept", make.names(cols))
    if(anyDuplicated(nom) || !all(nom %in% names(result$hits))) return(NULL)
    coefs <- as.matrix(result$hits[, nom, drop=F])
    dimnames(coefs) <- list(as.character(result$hits$feature), cols)
    return(coefs)
  }

  ## msqrob2 fits by robust regression and moderates variance rather than
  ##   coefficients, so each gene's fitted StatModel carries estimates for 
  ##   parameters msqrob2::msqrob() built from config$frm. A model that cannot 
  ##   be fit gives NA:

  if(method %in% c("msqrob", "msqrob_agg")) {

    dat <- SummarizedExperiment::rowData(result$fit[["genes"]])
    if(!all(c("msqrobModels", config$gene_id_col) %in% names(dat))) return(NULL)

    ## mixed fit carries estimates of random effects alongside fixed ones:

    nom <- f.msqrob_parms(cols,
      method %in% "msqrob_agg" && isTRUE(config$test_ridge))

    models <- dat$msqrobModels
    coefs <- matrix(as.numeric(NA), nrow=length(models), ncol=length(cols),
      dimnames=list(as.character(dat[[config$gene_id_col]]), cols))

    for(idx in seq_along(models)) {
      beta <- try(msqrob2::getCoef(models[[idx]]), silent=T)
      if(inherits(beta, "try-error") || is.null(names(beta))) next
      coefs[idx, ] <- beta[match(nom, names(beta))]
    }

    return(coefs)
  }

  ## proDA::proDA() renames intercept column:

  if(method %in% "proda") {
    coefs <- stats::coefficients(result$fit)
    nom <- ifelse(cols %in% "(Intercept)", "Intercept", cols)
    if(!is.matrix(coefs) || !all(nom %in% colnames(coefs))) return(NULL)
    coefs <- coefs[, nom, drop=F]
    colnames(coefs) <- cols
    return(coefs)
  }

  ## prolfqua::strategy_lm() fits a plain stats::lm() of response on columns
  ##   of design$X, under the make.names() forms of their names that test_prolfqua()
  ##   built, one model per feature in fit$modelDF. A model that could not be fit, or
  ##   that dropped a column as non-estimable for that feature, yields NA:

  ## mixed path fits one model per gene; falls back to stats::lm() for a gene
  ##   with a single feature:

  if(method %in% "prolfqua_lmer") {
    coefs <- result$coefs
    if(!is.matrix(coefs) || !all(cols %in% colnames(coefs))) return(NULL)
    return(coefs[, cols, drop=F])
  }

  if(method %in% "prolfqua") {

    mdf <- as.data.frame(result$fit$modelDF)
    id_col <- intersect(c(config$feat_col, config$gene_id_col, config$feat_id_col),
      names(mdf))
    if(!length(id_col) || !("linear_model" %in% names(mdf))) return(NULL)

    nom <- make.names(colnames(design$X), unique=T)[design$cols_test]
    coefs <- matrix(as.numeric(NA), nrow=nrow(mdf), ncol=length(cols),
      dimnames=list(as.character(mdf[[id_col[1]]]), cols))

    for(idx in seq_len(nrow(mdf))) {
      beta <- try(stats::coef(mdf$linear_model[[idx]]), silent=T)
      if(inherits(beta, "try-error") || is.null(names(beta))) next
      coefs[idx, ] <- beta[match(nom, names(beta))]
    }

    return(coefs)
  }

  return(NULL)
}

## average feature expression for the expr column of standardized table: 

f.feature_means <- function(state, method, config) {

  out <- rowMeans(state$expression, na.rm=T)
  out[is.nan(out)] <- NA                  ## a feature seen in no observation

  if(f.gene_level_method(method)) {
    return(tapply(out, f.gene_ids(state$features, config, "f.feature_means"),
      mean, na.rm=T))
  }

  names(out) <- as.character(state$features[[f.test_id_col(method, config)]])

  return(out)
}

## Fill two columns of standardized table so an engine's result table need not:

f.fill_standard <- function(tbl2, result, state, design, method, config) {

  ## logfc, only where engine reported none:

  if(any(is.na(tbl2$logfc))) {

    coefs <- f.test_coefs(result, method, design, config)

    if(is.null(coefs)) {
      f.msg("WARNING: test: method", method, "does not expose the coefficients of",
        "the columns carrying the test, so the logfc column is left empty",
        config=config)
    } else {
      eff <- f.logfc_effect(coefs, design, config)
      i <- is.na(tbl2$logfc) & tbl2$feature %in% names(eff)
      tbl2$logfc[i] <- eff[tbl2$feature[i]]
      f.msg("test: logfc:", if(!is.null(design$contrast)) {
        paste("config$contrast", trimws(config$contrast), ", the weighted sum of",
          length(design$cols_test), "coefficients")
      } else if(length(design$cols_test) %in% 1) {
        paste("coefficient of", colnames(design$X)[design$cols_test])
      } else {
        paste("total swing over", length(design$cols_test), "coefficients (",
          paste(colnames(design$X)[design$cols_test], collapse=", "), ")")
      }, config=config)
    }
  }

  if(all(is.na(tbl2$expr))) {
    means <- f.feature_means(state, method, config)
    i <- tbl2$feature %in% names(means)
    tbl2$expr[i] <- means[tbl2$feature[i]]
  }

  return(tbl2)
}

## helper for test_h0():

f.format_lm <- function(tbl, id_col, config) {

  if(!is.data.frame(tbl)) {
    f.err("f.format_lm: !is.data.frame(tbl); class(tbl): ",
      class(tbl), config=config)
  }

  if(!(id_col %in% names(tbl))) {
    f.err("f.format_lm: id_col:", id_col, "not %in% names(tbl); names(tbl):",
      names(tbl), config=config)
  }

  nom <- c("pval", "p.adj", "stat")
  if(!all(nom %in% names(tbl))) {
    f.err("f.format_lm: expected names not %in% names(tbl); names(tbl):",
      names(tbl), "; expected names:", nom, config=config)
  }

  ## stat is the F statistic test_lm() reports for its model comparison: 
  
  tbl <- data.frame(feature=tbl[[id_col]], expr=as.numeric(NA),
    logfc=as.numeric(NA), stat=tbl$stat, lod=as.numeric(NA),
    pval=tbl$pval, adj_pval=tbl$p.adj)
  
  tbl <- tbl[order(tbl$pval, decreasing=F), , drop=F]
  rownames(tbl) <- NULL
  
  return(tbl)
}

## helper for test_h0():

f.format_msqrob <- function(tbl, id_col, config) {

  if(!is.data.frame(tbl)) {
    f.err("f.format_msqrob: !is.data.frame(tbl); class(tbl): ",
      class(tbl), config=config)
  }

  if(!(id_col %in% names(tbl))) {
    f.err("f.format_msqrob: id_col:", id_col, "not %in% names(tbl); names(tbl):",
      names(tbl), config=config)
  }

  ## two shapes, since msqrob2::hypothesisTest() returns a moderated t and a log fold
  ##   change for the single contrast test_msqrob() runs when one design column carries
  ##   the test, and f.msqrob_wald() returns an F statistic and no single coefficient
  ##   when several do:

  if(all(c("logFC", "t", "pval", "adjPval") %in% names(tbl))) {

    tbl <- data.frame(feature=tbl[[id_col]], expr=as.numeric(NA),
      logfc=tbl$logFC, stat=tbl$t, lod=as.numeric(NA),
      pval=tbl$pval, adj_pval=tbl$adjPval)

  } else if(all(c("f_statistic", "pval", "adjPval") %in% names(tbl))) {

    tbl <- data.frame(feature=tbl[[id_col]], expr=as.numeric(NA),
      logfc=as.numeric(NA), stat=tbl$f_statistic, lod=as.numeric(NA),
      pval=tbl$pval, adj_pval=tbl$adjPval)

  } else {
    f.err("f.format_msqrob: expected names not %in% names(tbl); names(tbl):",
      names(tbl), config=config)
  }

  tbl <- tbl[order(tbl$pval, decreasing=F), , drop=F]
  rownames(tbl) <- NULL
  
  return(tbl)
}

## helper for test_h0():

f.format_proda <- function(tbl, config) {
  
  if(!is.data.frame(tbl)) {
    f.err("f.format_proda: !is.data.frame(tbl); class(tbl): ", 
      class(tbl), config=config)
  }

  ## two shapes, since proDA::test_diff() returns a t statistic and a difference
  ##   for the single contrast test_proda() runs when one design column carries the
  ##   test, and an F statistic with no difference for the likelihood ratio test it
  ##   runs when several do:

  if(all(c("name", "avg_abundance", "diff", "t_statistic", "pval", "adj_pval") %in%
    names(tbl))) {

    tbl <- data.frame(feature=tbl$name, expr=tbl$avg_abundance,
      logfc=tbl$diff, stat=tbl$t_statistic, lod=as.numeric(NA),
      pval=tbl$pval, adj_pval=tbl$adj_pval)

  } else if(all(c("name", "avg_abundance", "f_statistic", "pval", "adj_pval") %in%
    names(tbl))) {

    tbl <- data.frame(feature=tbl$name, expr=tbl$avg_abundance,
      logfc=as.numeric(NA), stat=tbl$f_statistic, lod=as.numeric(NA),
      pval=tbl$pval, adj_pval=tbl$adj_pval)

  } else {
    f.err("f.format_proda: expected names not %in% names(tbl); names(tbl):",
      names(tbl), config=config)
  }

  tbl <- tbl[order(tbl$pval, decreasing=F), , drop=F]
  rownames(tbl) <- NULL
  
  return(tbl)
}

## helper for test_h0():

f.format_prolfqua <- function(tbl, id_col, config) {

  if(!is.data.frame(tbl)) {
    f.err("f.format_prolfqua: !is.data.frame(tbl); class(tbl): ",
      class(tbl), config=config)
  }

  if(!(id_col %in% names(tbl))) {
    f.err("f.format_prolfqua: id_col:", id_col,
      "not %in% names(tbl); names(tbl):", names(tbl), config=config)
  }

  nom <- c("F.value", "p.value", "FDR")
  if(!all(nom %in% names(tbl))) {
    f.err("f.format_prolfqua: expected names not %in% names(tbl); names(tbl):", 
      names(tbl), config=config)
  }
  
  tbl <- data.frame(feature=tbl[[id_col]], expr=as.numeric(NA),
    logfc=as.numeric(NA), stat=tbl$F.value, lod=as.numeric(NA),
    pval=tbl$p.value, adj_pval=tbl$FDR)
  
  tbl <- tbl[order(tbl$pval, decreasing=F), , drop=F]
  rownames(tbl) <- NULL
  
  return(tbl)
}

## helper for test_h0(). Separate from f.format_limma() because DEqMS's table carries two
##   sets of statistics: the limma columns of the moderated fit, and the sca.*
##   columns DEqMS's count-based prior produced:

f.format_deqms <- function(tbl, config) {

  if(!is.data.frame(tbl)) {
    f.err("f.format_deqms: !is.data.frame(tbl); class(tbl): ",
      class(tbl), config=config)
  }

  nom <- c("gene", "AveExpr", "sca.t", "sca.F", "sca.df.num", "sca.P.Value",
    "sca.adj.pval")

  if(!all(nom %in% names(tbl))) {
    f.err("f.format_deqms: expected names not %in% names(tbl); missing:",
      paste(setdiff(nom, names(tbl)), collapse=", "), "; names(tbl):",
      names(tbl), config=config)
  }

  df_num <- unique(tbl$sca.df.num)

  if(length(df_num) != 1) {
    f.err("f.format_deqms: the table reports", length(df_num), "different numerator",
      "degrees of freedom (", paste(utils::head(df_num, 5), collapse=", "),
      "), but the hypothesis is a property of config$frm and config$test_term and so",
      "is the same for every gene", config=config)
  }

  if(df_num %in% 1) {

    if(!("logFC" %in% names(tbl))) {
      f.err("f.format_deqms: one coefficient carries the test but the table has no",
        "logFC column; names(tbl):", names(tbl), config=config)
    }

    tbl <- data.frame(feature=tbl$gene, expr=tbl$AveExpr,
      logfc=tbl$logFC, stat=tbl$sca.t, lod=as.numeric(NA),
      pval=tbl$sca.P.Value, adj_pval=tbl$sca.adj.pval)

  } else {

    tbl <- data.frame(feature=tbl$gene, expr=tbl$AveExpr,
      logfc=as.numeric(NA), stat=tbl$sca.F, lod=as.numeric(NA),
      pval=tbl$sca.P.Value, adj_pval=tbl$sca.adj.pval)
  }

  tbl <- tbl[order(tbl$pval, decreasing=F), , drop=F]
  rownames(tbl) <- NULL

  return(tbl)
}

## helper for test_h0():

f.format_limma <- function(tbl, config) {

  if(!is.data.frame(tbl)) {
    f.err("f.format_limma: !is.data.frame(tbl); class(tbl): ",
      class(tbl), config=config)
  }
    
  if(all(c("AveExpr", "logFC", "t", "B", "P.Value", "adj.P.Val") %in% names(tbl))) {
  
    tbl <- data.frame(feature=rownames(tbl), expr=tbl$AveExpr, 
      logfc=tbl$logFC, stat=tbl$t, lod=tbl$B, 
      pval=tbl$P.Value, adj_pval=tbl$adj.P.Val)
    
  } else if(all(c("AveExpr", "F", "P.Value", "adj.P.Val") %in% names(tbl))) {
  
    tbl <- data.frame(feature=rownames(tbl), expr=tbl$AveExpr, 
      logfc=as.numeric(NA), stat=tbl$F, lod=as.numeric(NA), 
      pval=tbl$P.Value, adj_pval=tbl$adj.P.Val)
    
  } else {
    f.err("f.format_limma: expected names not %in% names(tbl); names(tbl):", 
      names(tbl), config=config)
  }
  
  tbl <- tbl[order(tbl$pval, decreasing=F), , drop=F]
  rownames(tbl) <- NULL
  
  return(tbl)
}

