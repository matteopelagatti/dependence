#' Orthonormal spline basis for cubic splines on the first n integers
#'
#' It computes an orthonormal mean-zero version of a cubic splines basis with
#' \eqn{p} degrees of freedom and uniform knots.
#'
#' @param n Integer: the largest integer in the sequence 1, 2, ..., n
#' @param p Integer: the number of basis functions (degrees of freedom)
#'
#' @returns A \eqn{n\times p} matrix.
#'
#' @examples
#' M <- spline_basis_int(100, 3)
#' crossprod(M)
#'
#' @export
ortho_spline_basis_int <- function(n, p) {
  quick_qr(splines::bs(seq(n), df = p))
}
