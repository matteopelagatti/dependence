// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace arma;

//' Discrete Chebyshev (Gram) polynomial basis
//'
//' Computes the orthonormal polynomial basis on the grid
//' {1/(n+1), 2/(n+1), ..., n/(n+1)} using the three-term
//' recurrence. Returns an n x p matrix U such that
//' (1/n) * t(U) %*% U = I_p exactly.
//' Cost: O(n*p) time, O(n*p) space.
//' No QR, no Vandermonde matrix needed.
//'
//' @param n positive integer with the number of points n
//' @param p positive integer with the degree of the polynomial
//'
//' @returns A \eqn{n \times p} matrix with the basis vectors.
//' @export
// [[Rcpp::export]]
arma::mat chebyshev_basis(const int n, const int p) {

  arma::mat U(n, p);

  // Work on the integer grid {1, 2, ..., n} rescaled
  // to {1/(n+1), ..., n/(n+1)}: the recurrence is the
  // same up to a trivial rescaling absorbed into the
  // normalization constants.
  // We use the integer grid internally for numerical
  // precision and normalize at the end.

  double mu = (n + 1.0) / 2.0; // mean of {1,...,n}

  // t_0(x) = 1 (constant, will be dropped since we
  // work with centred basis: but we keep it to drive
  // the recurrence and discard after)
  arma::vec t_prev(n), t_curr(n), t_next(n);
  t_prev.fill(1.0);              // t_0

  // t_1(x) = x - mu  (already centred)
  for (int i = 0; i < n; i++)
    t_curr(i) = (double)(i + 1) - mu; // t_1

  // Fill column 0 with t_1 (we skip the constant t_0)
  double norm0 = std::sqrt(arma::dot(t_curr, t_curr));
  U.col(0) = t_curr * (std::sqrt((double)n) / norm0);

  for (int k = 1; k < p; k++) {
    // Recurrence coefficients for discrete Chebyshev
    // polynomials on {1,...,n}:
    // alpha_k = mu  (by symmetry, same for all k)
    // beta_k  = k^2 * (n^2 - k^2) / (4 * (4k^2 - 1))
    double kd = (double)k;
    double beta = kd * kd * ((double)n * (double)n
                               - kd * kd)
      / (4.0 * (4.0 * kd * kd - 1.0));

      // t_{k+1}(x) = (x - mu)*t_k(x) - beta*t_{k-1}(x)
      for (int i = 0; i < n; i++) {
        double x = (double)(i + 1);
        t_next(i) = (x - mu) * t_curr(i)
          - beta   * t_prev(i);
      }

      // Normalize to obtain (1/n)||col||^2 = 1,
      // i.e. ||col||_2 = sqrt(n)
      double norm_k = std::sqrt(
        arma::dot(t_next, t_next));
      U.col(k) = t_next * (std::sqrt((double)n)
                             / norm_k);

      // Shift for next iteration
      t_prev = t_curr;
      t_curr = t_next;
  }

  return U;
}
