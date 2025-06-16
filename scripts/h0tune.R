## configure:
config1 <- h0testr::f.new_config()
config1$frm <- ~age+strain+gender+age:strain   ## formula with variable of interest and covariates
config1$test_term <- "age:strain"              ## term in FRM on which test is to be performed
config1$permute_var <- ""                      ## variable to permute; "" for no permutation (normal execution)
config1$sample_factors <- list(                ## set levels of factor variables in sample metadata
  age=c("4", "12", "24"),                      ## numeric treated as numeric unless levels set here, then as factor
  strain=c("B6", "A_J", "Balbc_J"),
  gender=c("Male", "Female")
)
config1$feat_id_col <- "gene_id"               ## features[, feat_id_col] == rownames(exprs)
config1$obs_col <- "assay_id"                  ## unique id for observations; samples[, obs_col] == colnames(exprs)
config1$sample_col <- "sample_id"              ## id for samples in samps; not unique if tech reps;
    
## h0test.1:
out <- h0testr::f.load_data(config1)
state1 <- out$state
config1 <- out$config

## h0test.2:
state2 <- state1
config2 <- config1
out <- h0testr::f.normalize(state2, config2)
out <- h0testr::f.combine_reps(out$state, out$state)
state2 <- out$state
config2 <- out$config

## h0test.3:
state3 <- state2
config3 <- config2
out <- h0testr::f.filter(state3, config3)
state3 <- f.impute(out$state, out$config)
tbl <- h0testr::f.test(state3, config3)

