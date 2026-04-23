// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

//' Quick QR on after de-meaning
//'
//' It computes the QR decomposition of a matrix B
//' after substracting the means of the columns of B.
//' It is much quicker than \code{qr.Q(qr(M))}.
//'
//' @param B matrix to be de-meaned and decomposed
//'
//' @returns The Q matrix of the decomposition.
// [[Rcpp::export]]
arma::mat quick_qr(arma::mat B) {
  // 1. Center the columns to have exactly mean zero
  // B.each_row() subtracts the column means from every row efficiently
  B.each_row() -= arma::mean(B, 0);

  arma::mat Q, R;

  // 2. Economic QR decomposition on the centered matrix
  arma::qr_econ(Q, R, B);

  // 3. Scale by sqrt(n) and return
  return Q;
}
