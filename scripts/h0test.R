## get default configuration:
config <- h0testr::f.new_config()

## customize:
config$frm <- ~age+strain+gender+age:strain   ## formula with variable of interest and covariates
config$test_term <- "age:strain"              ## term in FRM on which test is to be performed
config$permute_var <- ""                      ## variable to permute; "" for no permutation (normal execution)
config$reference_levels <- c(                  ## reference level of each factor variable
  age="4",                                    ## numeric treated as continuous unless named here
  strain="B6",                                ## remaining levels sorted: "A_J", "Balbc_J"
  gender="Male"
)
config$feat_id_col <- "gene_id"               ## features[, feat_id_col] == rownames(exprs)
config$obs_col <- "assay_id"                  ## unique id for observations; samples[, obs_col] == colnames(exprs)
config$sample_col <- "sample_id"              ## id for samples in samps; not unique if tech reps;

## run basic testing workflow:
out <- h0testr::f.run(config)


