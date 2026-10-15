## The normalisation of the bases, which is a contract rather than a detail.
##
## indeptest() and indepscan() take n^{-1} U'V to be the whitened cross-moment
## matrix whenever the ranks are untied, skipping the Cholesky step.  That is
## correct only if the basis satisfies n^{-1} B'B = I.  When ortho_spline_basis_int()
## did not -- its columns had unit sum of squares rather than unit mean square --
## every spline-basis statistic on continuous data was too small by a factor of n
## per margin and the test had size 0.000 at the nominal 5%, silently, because the
## construction is otherwise sound and nothing else looks wrong.
##
## The property is exact and costs nothing to check, so it is checked rather
## than relied upon.  Plain stopifnot(), so R CMD check runs it without a test
## dependency.

library(dependence)

## ------------------------------------------------------------------------
## 1.  Both bases, several sizes and orders: mean-zero and grid-orthonormal.
for (n in c(37L, 100L, 500L))
  for (p in c(1L, 2L, 4L, 7L))
    for (bf in c(chebyshev_basis, ortho_spline_basis_int)) {
      B <- suppressWarnings(bf(n, p))
      stopifnot(is.matrix(B), nrow(B) == n, all(is.finite(B)))
      stopifnot(ncol(B) == p)            # orders <= 3 are now polynomials,
      k <- p                             # not bs()'s silent three columns
      stopifnot(max(abs(colMeans(B))) < 1e-8,
                max(abs(crossprod(B) / n - diag(k))) < 1e-8)
    }

## 2.  The property the untied fast path actually uses: with ranks that are a
##     permutation of the grid, the design is orthonormal after permutation too.
set.seed(20261014)
for (bf in c(chebyshev_basis, ortho_spline_basis_int)) {
  n <- 200L
  B <- suppressWarnings(bf(n, 4L))
  r <- sample.int(n)
  stopifnot(max(abs(crossprod(B[r, , drop = FALSE]) / n - diag(ncol(B)))) < 1e-8)
}

## 3.  End to end: the statistic must not depend on the basis scale.  Doubling
##     a basis cannot change a canonical correlation, so a correctly scaled
##     basis and the Cholesky route must agree.  This is what failed before.
wc <- dependence:::whiten_cross
set.seed(1)
n <- 300L; p <- 3L
x <- runif(n); y <- sin(4 * pi * x) + 0.5 * rnorm(n)
rx <- rank(x, ties.method = "first"); ry <- rank(y, ties.method = "first")
for (bf in c(chebyshev_basis, ortho_spline_basis_int)) {
  B <- suppressWarnings(bf(n, p))
  U <- B[rx, , drop = FALSE]; V <- B[ry, , drop = FALSE]
  fast <- n * sum(wc(U, V, FALSE)^2)    # the untied path
  chol <- n * sum(wc(U, V, TRUE)^2)     # scale-free, hence the reference
  stopifnot(abs(fast - chol) < 1e-6 * max(1, chol))
}

## 4.  And the test itself is calibrated on both bases.  A short size check:
##     loose enough never to fail by chance, tight enough to catch a statistic
##     that is wrong by a factor of n, which gives size 0.
set.seed(99)
for (bas in c("poly", "spline")) {
  pv <- vapply(seq_len(600L), function(i)
    indeptest(runif(150), runif(150), p = 3L, q = 3L,
              basis = bas, ties = "first")$P_pvalue, 0)
  size <- mean(pv <= 0.05)
  stopifnot(size > 0.015, size < 0.11)   # se is 0.009; this is +-5 se
}

cat("basis scale: all checks passed\n")
