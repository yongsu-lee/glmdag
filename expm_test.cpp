#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
mat expm_test(mat X) {
  
  mat result = expmat(X);
  return(result);
  
}