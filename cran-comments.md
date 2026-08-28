# cran-comments

## Test environments

- Windows 11, R 4.3.1 (local)

## R CMD check results

0 errors | 0 warnings | 3 notes

### Note: New submission

This is a new submission.

### Note: Suggests or Enhances not in mainstream repositories: prolfqua

`prolfqua` is one of several optional hypothesis-testing engines. It is available only from
GitHub (<https://github.com/wolski/prolfqua>), so it cannot be named in
`Additional_repositories:`, which requires a repository URL rather than a source tree.

It is used conditionally throughout, as required by the CRAN policy on Suggests:

- Nothing in the package attaches it. The only uses are `prolfqua::` calls inside
  `test_prolfqua()`, which begins by checking that the package can be loaded and, if not,
  stops with the command needed to install it.
- Its example is wrapped in `requireNamespace("prolfqua", quietly = TRUE)`, so on a machine
  without prolfqua the example runs to completion having done nothing. This was verified by
  substituting an absent package name into the guard: the example exits 0 with no output.
- The package builds, checks and runs its full default workflow with prolfqua absent. The
  default testing engine is `limma`, which is a hard dependency.

The same pattern covers every other optional engine (DEqMS, glmnet, impute, imputeLCMD,
lme4, lmerTest, missForest, msqrob2, pcaMethods, proDA, QFeatures, randomForest,
SummarizedExperiment, vsn); those are all in mainstream repositories and so are not flagged.

### Note: unable to verify current time

The check machine has no network time source. This is a property of the local environment,
not of the package.

### Note: examples with elapsed time > 5s

The flagged examples exercise the slower optional engines and the configuration-tuning
function, which fits many workflow combinations by design.
