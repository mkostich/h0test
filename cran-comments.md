# cran-comments

## Test environments

- Debian GNU/Linux 13 (trixie), R 4.6.1 (apptainer container, `--as-cran --no-manual`)
- Windows 11, R 4.3.1 (local)

## R CMD check results

0 errors | 0 warnings | 2 notes

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
lme4, lmerTest, lmtest, missForest, msqrob2, pcaMethods, proDA, QFeatures, randomForest,
SummarizedExperiment, vsn); those are all in mainstream repositories and so are not flagged.

### Note: examples with elapsed time > 5s

One example is flagged: `test_msqrob` at 12.2s user, 11.8s elapsed. It fits one mixed model
per gene through `msqrob2`, which is the slowest of the optional engines. The example is
already reduced to the smallest input that exercises both of its paths.
