## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new submission.

## Test environments

* local Windows 11 x64, R 4.5.2
* win-builder: R-devel and R-release
* macOS builder (mac.r-project.org/macbuilder)

## Notes for the reviewers

* The examples of `indepscan()` and `depanalysis()` run in under 2 seconds
  each, so none is wrapped in `\donttest{}`.
* `indepscan_null()` calls `set.seed()` only when the user passes its
  `seed` argument, to make a null table reproducible; by default it leaves
  the random number generator state untouched apart from drawing from it.
