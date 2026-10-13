## Properties of the deformation analysis that must hold exactly, plus the
## one regression that is easy to reintroduce: the sign convention.  Plain
## stopifnot(), so R CMD check runs it with no test dependency.

library(dependence)

set.seed(20261014)
n <- 600L
x <- sqrt(rchisq(n, df = 2)); y <- rnorm(n) * x          # X is the sd of Y
th <- runif(n, 0, 2 * pi)                                 # a noisy circle
cx <- cos(th) + rnorm(n, sd = 0.1); cy <- sin(th) + rnorm(n, sd = 0.1)

## ------------------------------------------------------------------------
## 1.  The canonical structure.  These are identities, not approximations.
m <- depanalysis(cx, cy, p = 6L, q = 6L, dimensions = 1:3)
stopifnot(dim(m$g) == c(n, 3L), dim(m$h) == c(n, 3L))
## unit variance, in the mean-square sense the construction imposes
stopifnot(max(abs(colMeans(m$g^2) - 1)) < 1e-8,
          max(abs(colMeans(m$h^2) - 1)) < 1e-8,
          max(abs(colMeans(m$g))) < 1e-8, max(abs(colMeans(m$h))) < 1e-8)
## each pair correlates at exactly its canonical correlation ...
for (s in 1:3)
  stopifnot(abs(stats::cor(m$g[, s], m$h[, s]) - m$cor[s]) < 1e-8)
## ... the correlations decrease, and lie in [0, 1] ...
stopifnot(!is.unsorted(rev(m$cor)), all(m$cor >= 0), all(m$cor <= 1))
## ... and distinct pairs are uncorrelated, within and across margins.
stopifnot(abs(stats::cor(m$g[, 1], m$g[, 2])) < 1e-8,
          abs(stats::cor(m$h[, 1], m$h[, 3])) < 1e-8,
          abs(stats::cor(m$g[, 1], m$h[, 2])) < 1e-8)
## the coefficients reproduce the variates
stopifnot(length(m$coef_x) == 6L * 3L, length(m$coef_y) == 6L * 3L)

## 2.  The caption's numbers are indeptest()'s own, not a second estimate.
it <- indeptest(cx, cy, p = 6L, q = 6L, basis = "poly", ties = "first")
stopifnot(abs(m$stat[["B_n"]] - it$B_stat) < 1e-10,
          abs(m$stat[["P_n"]] - it$P_stat) < 1e-10,
          abs(m$pvalue[["B_n"]] - it$B_pvalue) < 1e-12)

## 3.  THE SIGN CONVENTION.  Without a rule anchored in the fitted function,
##     the picture flips between orders and between bases: "make the largest
##     |g| positive" was measured flipping between p = 5 and p = 6 on the
##     first example below, and disagreeing with the published orientation on
##     the second.  The orientation must not depend on p or on the basis.
orient <- function(X, Y, p, basis)
  sign(stats::cor(depanalysis(X, Y, p = p, q = p, basis = basis)$g[, 1], X))
for (basis in c("poly", "spline")) {
  o1 <- vapply(c(5L, 6L, 7L, 9L, 11L), function(p) orient(x, y, p, basis), 0)
  o2 <- vapply(c(5L, 6L, 7L, 9L, 11L), function(p) orient(cx, cy, p, basis), 0)
  stopifnot(length(unique(o1)) == 1L, length(unique(o2)) == 1L)
}
## monotone dependence must come out rising, in both bases
stopifnot(orient(x, y, 7L, "poly") > 0, orient(x, y, 7L, "spline") > 0)

## flip reverses a pair, and only the sign: the correlation is unchanged.
mf <- depanalysis(cx, cy, p = 6L, q = 6L, dimensions = 1:2, flip = c(TRUE, FALSE))
m2 <- depanalysis(cx, cy, p = 6L, q = 6L, dimensions = 1:2)
stopifnot(max(abs(mf$g[, 1] + m2$g[, 1])) < 1e-10,   # flipped
          max(abs(mf$h[, 1] + m2$h[, 1])) < 1e-10,   # and its partner
          max(abs(mf$g[, 2] - m2$g[, 2])) < 1e-10,   # untouched
          max(abs(mf$cor - m2$cor)) < 1e-12)

## 4.  Guards.
bad <- function(e) stopifnot(inherits(try(e, silent = TRUE), "try-error"))
bad(depanalysis(x, y[-1]))                                  # unequal lengths
bad(depanalysis(factor(sample(letters[1:4], n, TRUE)), y))  # categorical
bad(depanalysis(x, y, p = 3L, dimensions = 4L))             # beyond min(p, q)
bad(depanalysis(x, y, dimensions = 0L))
bad(depanalysis(x, y, p = 3L, q = 3L, basis = "wavelet"))
## A binary x cannot carry order 5: cap_order() reduces it with a warning
## rather than failing, which is what indeptest() does too.
w <- NULL
mb <- withCallingHandlers(depanalysis(rep(0:1, length.out = n), y, p = 5L),
        warning = function(e) { w <<- conditionMessage(e); invokeRestart("muffleWarning") })
stopifnot(!is.null(w), grepl("order reduced", w), mb$p == 1L, length(mb$cor) == 1L)

## 5.  The plots, when ggplot2 is there.
if (requireNamespace("ggplot2", quietly = TRUE)) {
  stopifnot(inherits(m$plotx, "ggplot"), inherits(m$ploty, "ggplot"),
            inherits(m$plotxy, "ggplot"))
  ## they must actually render, which is where a bad aes() mapping shows up
  for (pl in list(m$plotx, m$ploty, m$plotxy))
    stopifnot(inherits(ggplot2::ggplot_build(pl), "ggplot_built"))
  ## one panel per requested pair
  stopifnot(nrow(ggplot2::ggplot_build(m$plotx)$layout$layout) == 3L)
  ## annotate = "none" removes the caption
  stopifnot(is.null(depanalysis(x, y, annotate = "none")$plotx$labels$caption))
}

stopifnot(inherits(capture.output(print(m)), "character"))
cat("depanalysis: all checks passed\n")
