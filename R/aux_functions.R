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
  ## The sqrt(n) is load-bearing; see the details section above.
  quick_qr(splines::bs(seq(n), df = p)) * sqrt(n)
}
