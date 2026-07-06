// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace arma;
//' Discrete Chebyshev (Gram) polynomial basis
//'
//' Computes the orthonormal polynomial basis on the grid
//' 1/(n+1), 2/(n+1), ..., n/(n+1) using the *normalized*
//' three-term recurrence. Returns an n x p matrix U such that
//' (1/n) * t(U) %*% U = I_p exactly.
//' Cost: O(n*p) time, O(n*p) space.
//' No QR, no Vandermonde matrix needed.
//'
//' Unlike the monic recurrence, whose values grow like (n/2)^k
//' and overflow double precision for large n (e.g. n = 1e6,
//' p = 62), the normalized recurrence propagates the
//' orthonormal polynomials themselves, whose grid values grow
//' only polynomially in the degree; the computed basis is
//' orthonormal to machine precision throughout the regime
//' p = o(sqrt(n)).
//'
//' @param n positive integer with the number of points n
//' @param p positive integer with the degree of the polynomial
//'          (must satisfy p < n)
//'
//' @returns A \eqn{n \times p} matrix with the basis vectors
//'          of degrees 1, ..., p (the constant is omitted).
//' @export
// [[Rcpp::export]]
arma::mat chebyshev_basis(const int n, const int p) {
   if (p >= n)
     Rcpp::stop("chebyshev_basis: p must be smaller than n");
   arma::mat U(n, p);
   const double nd = (double)n;
   const double mu = (nd + 1.0) / 2.0;       // mean of {1,...,n}
   // Normalization convention: each phi_k has squared norm n on
   // the grid, i.e. (1/n)*||phi_k||^2 = 1, so that
   // (1/n) t(U) U = I_p.
   // Normalized recurrence (phi_k orthonormal, beta_k the monic
   // recurrence coefficients):
   //   phi_0(x) = 1
   //   phi_1(x) = (x - mu) / sqrt(beta_1),  beta_1 = (n^2-1)/12
   //   phi_{k+1}(x) = { (x - mu) phi_k(x)
   //                    - sqrt(beta_k) phi_{k-1}(x) } / sqrt(beta_{k+1})
   //   beta_k = k^2 (n^2 - k^2) / {4 (4 k^2 - 1)}
   arma::vec x = arma::regspace<arma::vec>(1.0, nd); // {1,...,n}
   arma::vec phi_prev(n, arma::fill::ones);          // phi_0
   const double beta1 = (nd * nd - 1.0) / 12.0;
   arma::vec phi_curr = (x - mu) / std::sqrt(beta1); // phi_1
   U.col(0) = phi_curr;
   double sb_prev = std::sqrt(beta1);                // sqrt(beta_k), k = 1
   for (int k = 1; k < p; k++) {
     // beta_{k+1}: coefficient for the degree being created
     const double kd = (double)(k + 1);
     const double beta = kd * kd * (nd * nd - kd * kd)
       / (4.0 * (4.0 * kd * kd - 1.0));
     const double sb = std::sqrt(beta);
     arma::vec phi_next = ((x - mu) % phi_curr
                             - sb_prev * phi_prev) / sb;
     U.col(k) = phi_next;
     // Shift for next iteration
     phi_prev = phi_curr;
     phi_curr = phi_next;
     sb_prev  = sb;
   }
   // Optional polish: rescale each column by its computed norm
   // to remove the O(1e-14) floating-point drift, restoring
   // ||col||_2 = sqrt(n) exactly. Cost O(n*p).
   const double sn = std::sqrt(nd);
   for (int k = 0; k < p; k++)
     U.col(k) *= sn / arma::norm(U.col(k), 2);
   return U;
}
