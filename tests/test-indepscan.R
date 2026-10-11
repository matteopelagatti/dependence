## Exactness checks for the order scan.  Plain stopifnot(), so R CMD check
## runs them without adding a test dependency.  Everything here is fast; the
## calibration and power study is in inst/paper/verify_indepscan.R.

library(dependence)
sp <- dependence:::scan_profile

## ------------------------------------------------------------------------
## 1.  The identity the scan rests on.  The k-th entry of the profile must
##     be the fixed-order statistic at p = q = k, for every basis, with and
##     without ties.  For the polynomial basis the profile comes from a
##     single pass over one cross-product; for the spline basis, which is
##     not nested, from one build per order.  Both must agree with
##     indeptest() to machine precision.
set.seed(20261011)
K <- 5L
mk <- list(
  function() { x <- runif(80);  list(x, sin(4 * pi * x) + 0.4 * rnorm(80)) },
  function() { x <- runif(60);  list(x, 4 * (x - 0.5)^2 + 0.3 * rnorm(60)) },
  function() { x <- abs(rnorm(120)); list(x, rnorm(120) * x) },
  function() { x <- sample(1:9, 100, TRUE); list(x, x + sample(1:9, 100, TRUE)) }
)
for (bas in c("poly", "spline"))
  for (tm in c("first", "last", "average"))
    for (g in mk) {
      d <- g()
      prof <- sp(d[[1]], d[[2]], K, basis = bas, ties = tm)
      ref <- vapply(seq_len(K), function(k) suppressWarnings(
        indeptest(d[[1]], d[[2]], p = k, q = k, basis = bas,
                  ties = tm)$P_stat), 0)
      stopifnot(length(prof) == K, all(is.finite(prof)),
                max(abs(prof - ref)) < 1e-9)
    }

## The polynomial basis is nested and the spline basis is not: this is the
## reason scan_profile() has two routes, so it is worth asserting.
B6 <- chebyshev_basis(100L, 6L)
stopifnot(max(abs(chebyshev_basis(100L, 3L) - B6[, 1:3])) < 1e-12)
S6 <- suppressWarnings(ortho_spline_basis_int(100L, 6L))
S4 <- suppressWarnings(ortho_spline_basis_int(100L, 4L))
stopifnot(max(abs(S4 - S6[, 1:4])) > 1e-6)

## ------------------------------------------------------------------------
## 2.  The table and the test agree on their own conventions.
tab <- indepscan_null(n = 50L, K = 3L, R = 400L, seed = 1)
x <- runif(50); y <- sin(4 * pi * x) + 0.4 * rnorm(50)
r <- indepscan(x, y, tab)
stopifnot(inherits(r, "indepscan"),
          length(r$profile) == 3L,
          length(r$cell.pvalues) == 3L,
          r$p.value > 0, r$p.value <= 1,
          all(r$cell.pvalues >= 0), all(r$cell.pvalues <= 1),
          r$k.min %in% 1:3,
          abs(sum(r$weights) - 1) < 1e-12,
          identical(r$k.min, which.min(r$cell.pvalues / r$weights)))
## The profile inside the object is the profile, not something else.
stopifnot(max(abs(r$profile - sp(x, y, 3L, "poly", "first"))) < 1e-12)

## K = 1 degenerates correctly: one cell, so the scan is that cell's test.
t1 <- indepscan_null(n = 40L, K = 1L, R = 200L, seed = 2)
r1 <- indepscan(runif(40), runif(40), t1)
stopifnot(length(r1$profile) == 1L,
          abs(r1$statistic - r1$cell.pvalues[1]) < 1e-12)

## Weights may be changed without retabulating, and only change the test.
rg <- indepscan(x, y, tab, weights = c(4, 2, 1))
stopifnot(max(abs(rg$profile - r$profile)) < 1e-12,
          max(abs(rg$cell.pvalues - r$cell.pvalues)) < 1e-12,
          abs(sum(rg$weights) - 1) < 1e-12)
## and "geometric" is the named form of 2^-k
stopifnot(max(abs(indepscan(x, y, tab, weights = "geometric")$weights -
                  indepscan(x, y, tab, weights = 2^-(1:3))$weights)) < 1e-12)

## The memoised null distribution must not change the answer.
stopifnot(abs(indepscan(x, y, tab)$p.value - r$p.value) < 1e-12)

## ------------------------------------------------------------------------
## 3.  Guards.  A mismatch between table and data must stop, never be
##     silently absorbed: the profile is read against the table column by
##     column, so a shortened or misaligned profile is a wrong p-value.
bad <- function(expr) stopifnot(inherits(try(expr, silent = TRUE), "try-error"))
bad(indepscan(runif(49), runif(49), tab))                     # wrong n
bad(indepscan(runif(50), runif(50), tab, basis = "spline"))   # wrong basis
bad(indepscan(runif(50), runif(49), tab))                     # unequal lengths
bad(indepscan(runif(50), runif(50)))                          # no table
bad(indepscan(runif(50), runif(50), list(n = 50)))            # not a table
bad(indepscan(runif(50), runif(50), tab, weights = c(1, -1, 1)))
bad(indepscan(runif(50), runif(50), tab, weights = c(1, 1)))
bad(indepscan(runif(50), runif(50), tab, weights = c(0, 1, 1)))
bad(indepscan_null(n = 10L, K = 10L, R = 10L))                # K >= n
bad(indepscan_null(n = 50L, K = 0L, R = 10L))                 # K < 1
## K beyond what the data can support: a binary x cannot carry order 3.
bad(indepscan(rep(0:1, 25), runif(50), tab))

## Ties in the data make the tabulated null approximate, and that is said.
tt <- tryCatch(indepscan(sample(1:6, 50, TRUE), runif(50), tab),
               warning = function(w) conditionMessage(w))
stopifnot(is.character(tt), grepl("ties", tt))

## ------------------------------------------------------------------------
## 4.  The printed output exists and mentions the p-value.
out <- capture.output(print(r))
stopifnot(length(out) > 3L, any(grepl("p", out)))

cat("indepscan: all exactness checks passed\n")
