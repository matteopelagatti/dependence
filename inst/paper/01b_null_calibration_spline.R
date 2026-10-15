## =====================================================================
##  01b_null_calibration_spline.R
##
##  The order calibration of Section S4.2 for the SPLINE basis: the largest
##  order p = q at which the chi-squared calibration of B_n keeps the size
##  within tolerance, for n = 25, ..., 1600.  It runs 01_null_calibration.R
##  unchanged except for the basis, so the design, seeds, safeguards and
##  reporting are identical.
##
##  Output: results/pq_calibration_v3_spline.csv and _spline_curve.csv.
##  The polynomial files, pq_calibration_v3*.csv, are not touched.
##
##  A BUILT-IN CHECK.  The data streams of that script are keyed on the
##  sample size and the block, never on the order, and up to order 3 the
##  spline basis is the polynomial basis (the package now builds it as one).
##  So for p <= 3 this run must reproduce the polynomial counts exactly --
##  provided it runs with the same number of cores as the polynomial run,
##  since the cores share out the draws.  Compare the size columns of the two
##  _curve.csv files at p <= 3: any difference means something else changed.
##
##  Requires the package version with the spline fixes of October 2026
##  (ortho_spline_basis_int() returning polynomials for p <= 3).
## =====================================================================
stopifnot("install the current package first: orders <= 3 must be polynomials" =
  identical(dependence::ortho_spline_basis_int(50L, 2L), dependence::chebyshev_basis(50L, 2L)))

.here <- local({
  a <- commandArgs(FALSE); m <- grep("^--file=", a, value = TRUE)
  if (length(m)) return(dirname(normalizePath(sub("^--file=", "", m[1L]))))
  for (i in rev(seq_len(sys.nframe()))) {
    f <- tryCatch(sys.frame(i)$ofile, error = function(e) NULL)
    if (is.character(f) && length(f) == 1L) return(dirname(normalizePath(f)))
  }
  if (requireNamespace("rstudioapi", quietly = TRUE) &&
      isTRUE(tryCatch(rstudioapi::isAvailable(), error = function(e) FALSE))) {
    p <- tryCatch(rstudioapi::getSourceEditorContext()$path, error = function(e) "")
    if (nzchar(p)) return(dirname(normalizePath(p)))
  }
  getwd()
})
.old <- options(dependence.calib.basis = "spline")
tryCatch(source(file.path(.here, "01_null_calibration.R"), local = FALSE),
         finally = options(.old))
