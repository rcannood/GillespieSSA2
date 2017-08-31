#include <Rcpp.h>
using namespace Rcpp;

//[[Rcpp::export]]
NumericVector ssa_nutiling(NumericVector a, NumericMatrix nu, int j) {
  int M = nu.ncol();                  // Number of reaction channels in nu-tile
  int N = nu.nrow();                  // Number of states in nu tile
  int U = a.length() / M;             // Number of tessallations of nu tile
  int f = ceil((j / M) - 1);          // Frameshift factor
  int jp = j - f * M;                 // Relative reaction channel index
  NumericVector nu_j = rep(0.0, U*N);
  for (int i = 0; i < N; i++) {
    nu_j(i + f * N) = nu(i, jp);
  }
  return nu_j;                        // Return output
}
