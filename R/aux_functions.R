#' Grid-orthonormal spline basis for cubic splines on the first n integers
#'
#' It computes a mean-zero cubic B-spline basis with \eqn{p} degrees of
#' freedom and uniform knots, orthonormalised with respect to the uniform
#' measure on \eqn{\{1,\ldots,n\}}, so that
#' \deqn{n^{-1} B^\top B = I_p .}
#'
#' @details
#' The scaling is a contract, not a convention of convenience:
#' [chebyshev_basis()] returns a matrix in the same normalisation, and
#' [indeptest()] and [indepscan()] rely on it. When the ranks are untied
#' they take \eqn{S = n^{-1} U^\top V} to be the whitened cross-moment
#' matrix and skip the Cholesky step, which is correct only if the basis
#' satisfies the identity above.
#'
#' Before October 2026 this function returned `quick_qr()` unscaled, whose
#' columns have unit *sum* of squares rather than unit mean square, so that
#' \eqn{B^\top B = I} and \eqn{n^{-1} B^\top B = n^{-1} I}. Every statistic
#' computed from it on untied data was then too small by a factor of
#' \eqn{n} per margin, and the chi-squared test never rejected: measured
#' size was 0.000 at the nominal 5\% in all three affected paths
#' (numeric-numeric, factor-numeric and numeric-factor), against 0.049,
#' 0.046 and 0.055 after the correction. Tied data were never affected,
#' since those paths whiten by Cholesky, which is scale-free; nor was
#' [indepscan()], whose p-value comes from a simulated table in which a
#' common scale factor cancels. If you have saved results obtained with
#' `basis = "spline"` on continuous data, recompute them.
#'
#' For \eqn{p \le 3} the function returns [chebyshev_basis()]. A cubic
#' spline has at least three degrees of freedom, and with exactly three it
#' has no interior knot, so its span is the space of cubic polynomials and
#' every statistic is unchanged; below three there is no cubic spline, and the
#' natural meaning of an order-\eqn{p} spline is the polynomial of degree
#' \eqn{p}. Earlier versions passed \eqn{p < 3} to `splines::bs()`, which
#' silently used three degrees of freedom, and returned the leading columns of
#' that basis, which are not the linear and quadratic functions: with
#' `basis = "spline"` and \eqn{p = q = 1} the test was not a function of
#' Spearman's correlation, as it is now. The two families differ only from
#' \eqn{p = 4}, where the spline acquires its first interior knot.
#'
#' Note that [quick_qr()], used directly by the categorical methods of
#' [indeptest()], returns the \eqn{B^\top B = I} normalisation instead, and
#' those methods form \eqn{U^\top V} without dividing by \eqn{n}. The two
#' conventions are each self-consistent; this function must follow the
#' first, because that is what its callers assume.
#'
#' @param n Integer: the largest integer in the sequence 1, 2, ..., n
#' @param p Integer: the number of basis functions (degrees of freedom)
#'
#' @returns A \eqn{n\times p} matrix with mean-zero columns satisfying
#'   \eqn{n^{-1} B^\top B = I_p}.
#'
#' @seealso [chebyshev_basis()], which uses the same normalisation.
#'
#' @examples
#' M <- ortho_spline_basis_int(100, 3)
#' crossprod(M) / 100        # the identity, to machine precision
#' colMeans(M)               # zero
#'
#' @export
ortho_spline_basis_int <- function(n, p) {
  ## A cubic B-spline needs at least three degrees of freedom, and with exactly
  ## three it has no interior knot, so its span IS the cubic polynomials.  Below
  ## three there is no cubic spline at all: bs() silently returns three columns,
  ## and taking the first p of them gives a bump, not the linear or quadratic
  ## function.  An order-p spline with p <= 3 is the degree-p polynomial, so it
  ## is built as one; at p = 3 this changes no statistic, below 3 it makes p = 1
  ## Spearman's correlation, as it should be.
  if (p <= 3L) return(chebyshev_basis(n, p))
  ## The sqrt(n) is load-bearing; see the details section above.
  quick_qr(splines::bs(seq(n), df = p)) * sqrt(n)
}

#' Basis on the integer codes of a factor
#'
#' @description
#' Builds the design matrix for a factor coded by the integers of its levels,
#' in the [quick_qr()] normalisation (\eqn{U^\top U = I}) used by the
#' categorical methods of [indeptest()].
#'
#' Three rules, each of which fixes a defect of building the spline basis with
#' `splines::bs(z, df = k)` directly:
#' \itemize{
#'   \item For \eqn{k \le 3} the spline space is the polynomial space, and
#'     `bs()` cannot build fewer than three columns: it silently returned
#'     three, which the methods then used whole against \eqn{\chi^2_{kq}},
#'     so that a two- or three-level factor was tested with junk columns.
#'     Measured size at the nominal 5\% was 0.59 and 0.19 for factor against
#'     numeric, 1.00 and 0.99 for factor against factor.
#'   \item At the saturated order \eqn{k = c - 1}, any full-rank basis spans
#'     the centred indicators of the \eqn{c} levels, so every choice gives the
#'     same statistic, Pearson's; polynomials are used because they are always
#'     of full rank there.
#'   \item Otherwise the interior knots are equally spaced on the range of the
#'     codes, as `poly()` treats the levels as equally spaced. Knots at
#'     quantiles of the codes, the `bs()` default, coincide when a level is
#'     frequent, and left the basis rank-deficient in about half of the
#'     configurations tried.
#' }
#' A final rank check turns any remaining deficiency into an error rather
#' than a basis padded with arbitrary columns.
#'
#' @param z Integer codes of the factor levels.
#' @param k Order, at most the number of distinct codes minus one.
#' @param basis `"poly"` or `"spline"`.
#' @param label Name of the variable, for the error message.
#' @returns A matrix with `k` orthonormal, mean-zero columns.
#' @keywords internal
#' @noRd
factor_code_basis <- function(z, k, basis, label = "x") {
  k <- as.integer(k)
  m <- length(unique(z)) - 1L
  B <- if (basis == "poly" || k <= 3L || k >= m) {
    stats::poly(z, degree = k)
  } else {
    r  <- range(z)
    kn <- seq(r[1L], r[2L], length.out = k - 1L)[-c(1L, k - 1L)]
    splines::bs(z, knots = kn, degree = 3L, Boundary.knots = r)
  }
  B <- unclass(as.matrix(B))
  if (ncol(B) != k ||
      qr(B - rep(colMeans(B), each = nrow(B)))$rank < k)
    stop(sprintf(paste0("the %s basis on the codes of %s is rank-deficient at ",
                        "order %d; use basis = \"poly\" or \"dummy\", or a ",
                        "lower order."), basis, label, k), call. = FALSE)
  quick_qr(B)
}
