#include <Rcpp.h>
#include "ssa.h"

using namespace Rcpp;

class SSA_etl : public SSA {
public:
  SSA_etl(double tau_) : SSA("ETL"), tau(tau_) {}

  double tau;

  // preallocated data structures
  IntegerVector k;
  void allocate(const int M, const int N) {
    k = IntegerVector(M);
  }

  void step_matrix(
      const NumericVector& state,
      const NumericVector& propensity,
      const IntegerMatrix& nu,
      double* dtime,
      NumericVector& dstate
  ) {
    int M = nu.ncol();
    int N = nu.nrow();

    // determine reaction firings
    for (int j = 0; j < M; j++) {
      k[j] = R::rpois(propensity[j] * tau);
    }

    // determine firing effects
    for (int i = 0; i < N; i++) {
      double sum = 0.0;
      for (int j = 0; j < M; j++) {
        sum += nu(i, j) * k[j];
      }
      dstate[i] = sum;
    }

    // tau leap
    *dtime = tau;
  }

  void step_vector(
      const NumericVector& state,
      const NumericVector& propensity,
      const IntegerVector& nu_row,
      const IntegerVector& nu_effect,
      double* dtime,
      NumericVector& dstate
  ) {
    int M = propensity.size();
    int N = dstate.size();

    // clear dstate
    for (int i = 0; i < N; i++) {
      dstate[i] = 0;
    }

    for (int j = 0; j < M; j++) {
      // determine reaction firing
      k[j] = R::rpois(propensity[j] * tau);

      // determine firing effect
      dstate[nu_row[j]] = nu_effect[j] * k[j];
    }

    // tau leap
    *dtime = tau;
  }
} ;

// [[Rcpp::export]]
SEXP make_ssa_etl(double tau) {
  SSA_etl *ssa = new SSA_etl(tau);
  XPtr<SSA_etl> ptr(ssa);
  return ptr;
}
