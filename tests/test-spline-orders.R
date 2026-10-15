## Spline bases at small and unequal orders, and on the codes of a factor.
##
## Three defects, all of which left the size or the meaning of the test wrong
## without any error:
##
##  - splines::bs() cannot build a cubic spline with fewer than three degrees
##    of freedom.  It silently built three, so an "order-1 spline" was the
##    leading column of a cubic basis, a bump, and p = q = 1 was not
##    Spearman's test.  Up to order 3 the spline space IS the polynomial space,
##    so the order-p spline is now the degree-p polynomial.
##
##  - The factor methods used that three-column basis whole while referring
##    the statistic to chi-squared(p q): a two-level factor against a numeric
##    variable was rejected 59% of the time at the nominal 5%, and two
##    two-level factors 100% of the time.  Above order 3, knots at quantiles
##    of the codes coincided whenever a level was frequent, and the basis was
##    rank-deficient in about half of the configurations tried.
##
##  - With p != q the numeric method built one basis of order max(p, q) and
##    gave the smaller side its first columns.  That is right for nested
##    polynomials and wrong for splines, whose leading columns live on the
##    lower ranks.
##
## Plain stopifnot(), so R CMD check runs it without a test dependency.

library(dependence)
fcb <- dependence:::factor_code_basis
set.seed(20261015)

## ------------------------------------------------------------------------
## 1.  Orders up to 3: the spline basis is the polynomial basis.
for (n in c(25L, 100L, 333L))
  for (p in 1:3)
    stopifnot(identical(ortho_spline_basis_int(n, p), chebyshev_basis(n, p)))
for (p in 4:6)                                  # and from 4 it is a spline
  stopifnot(ncol(ortho_spline_basis_int(100L, p)) == p,
            max(svd(crossprod(ortho_spline_basis_int(100L, p),
                              chebyshev_basis(100L, p)) / 100)$d) < 1 + 1e-12,
            min(svd(crossprod(ortho_spline_basis_int(100L, p),
                              chebyshev_basis(100L, p)) / 100)$d) < 0.999)

## 2.  So the two families give the same statistics up to order 3, and
##     p = q = 1 is Spearman's test: P_n = n r_S^2.
x <- rnorm(150); y <- x^2 + rnorm(150)
for (p in 1:3) {
  a <- indeptest(x, y, p = p, q = p, basis = "spline")
  b <- indeptest(x, y, p = p, q = p, basis = "poly")
  stopifnot(abs(a$P_stat - b$P_stat) < 1e-9, abs(a$B_stat - b$B_stat) < 1e-9)
}
rs <- cor(rank(x), rank(y))
for (bs_ in c("poly", "spline"))
  stopifnot(abs(indeptest(x, y, p = 1, q = 1, basis = bs_)$P_stat -
                150 * rs^2) < 1e-9)

## 3.  Unequal orders with splines: each variable has the basis of its own
##     order, not the leading columns of the larger one.
n  <- 120L; x <- runif(n); y <- sin(6 * x) + rnorm(n, sd = 0.3)
for (pq in list(c(4L, 6L), c(6L, 4L), c(2L, 5L), c(5L, 1L))) {
  p <- pq[1]; q <- pq[2]
  U <- ortho_spline_basis_int(n, p)[rank(x), , drop = FALSE]
  V <- ortho_spline_basis_int(n, q)[rank(y), , drop = FALSE]
  P_manual <- n * sum((crossprod(U, V) / n)^2)
  P_pkg    <- indeptest(x, y, p = p, q = q, basis = "spline", ties = "first")$P_stat
  stopifnot(abs(P_pkg - P_manual) < 1e-9)
}

## 4.  Factors.  The basis has exactly the requested number of columns, and at
##     the saturated order every basis gives Pearson's statistic.
for (L in 2:9) {
  for (rep in 1:5) {
    pr <- rexp(L)^3                           # very unequal level frequencies
    f  <- factor(sample.int(L, 300, TRUE, prob = pr), levels = seq_len(L))
    f  <- droplevels(f); c_ <- nlevels(f); if (c_ < 2) next
    for (b in c("poly", "spline")) {
      U <- fcb(as.integer(f), c_ - 1L, b)
      stopifnot(ncol(U) == c_ - 1L, max(abs(crossprod(U) - diag(c_ - 1L))) < 1e-8)
    }
    g  <- factor(sample.int(3, 300, TRUE))
    st <- vapply(c("poly", "spline", "dummy"), function(b)
            indeptest(f, g, basis = b)$P_stat, 0)
    stopifnot(diff(range(st)) < 1e-8 * max(1, abs(st[1])))
    stopifnot(abs(st[1] - suppressWarnings(chisq.test(table(f, g),
                                            correct = FALSE)$statistic)) < 1e-6)
  }
}
##     Below the saturated order the spline is a genuine choice, built with
##     knots equally spaced on the codes, and must be of full rank.
for (L in 6:10) {
  f <- factor(sample.int(L, 400, TRUE, prob = rexp(L)^3), levels = seq_len(L))
  f <- droplevels(f); c_ <- nlevels(f)
  if (c_ >= 6) for (k in 4:(c_ - 2L)) {
    U <- fcb(as.integer(f), k, "spline")
    stopifnot(ncol(U) == k, max(abs(crossprod(U) - diag(k))) < 1e-8)
  }
}
##     And a deficient basis is an error, not a basis padded with junk.
##     (Codes of a real factor are consecutive, so this cannot happen there;
##     a gap that leaves a whole B-spline without data shows the guard works.)
gappy <- rep(c(1:6, 1000), length.out = 210)
stopifnot(inherits(tryCatch(fcb(gappy, 5L, "spline", "z"), error = function(e) e),
                   "error"))

## 5.  Size, the property that was broken: two- and three-level factors with
##     splines, against a numeric variable and against a factor.  A coarse
##     check -- the defect gave 0.19 to 1.00 -- kept short for R CMD check.
R <- 300L; n <- 120L
for (L in 2:3) {
  rej <- c(fn = 0, ff = 0, nf = 0)
  for (r in seq_len(R)) {
    f <- factor(sample.int(L, n, TRUE)); g <- factor(sample.int(L, n, TRUE)); z <- runif(n)
    rej["fn"] <- rej["fn"] + (indeptest(f, z, basis_fct = "spline", basis_num = "spline")$P_pvalue <= 0.05)
    rej["ff"] <- rej["ff"] + (indeptest(f, g, basis = "spline")$P_pvalue <= 0.05)
    rej["nf"] <- rej["nf"] + (indeptest(z, f)$P_pvalue <= 0.05)       # defaults
  }
  stopifnot(all(rej / R < 0.10))
}

## 6.  numeric x factor with default arguments runs, and is the mirror image
##     of factor x numeric.
z <- runif(80); f <- factor(sample(c("a", "b", "c"), 80, TRUE))
a <- indeptest(z, f); b <- indeptest(f, z)
stopifnot(abs(a$P_stat - b$P_stat) < 1e-12, abs(a$B_stat - b$B_stat) < 1e-12)

cat("spline orders: all checks passed\n")
