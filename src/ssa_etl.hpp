#ifndef DYNGEN_SSA_ETL_H
#define DYNGEN_SSA_ETL_H

#include <Rcpp.h>
#include "ssa.hpp"
#include "utils.hpp"

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

  void step(
      const NumericVector& state,
      const NumericVector& transition_rates,
      const IntegerMatrix& nu,
      double* dtime,
      NumericVector& dstate
  ) {
    int M = nu.ncol();
    int N = nu.nrow();

    for (int j = 0; j < M; j++) {
      k[j] = R::rpois(transition_rates[j] * tau);
    }

    for (int i = 0; i < N; i++) {
      double sum = 0.0;
      for (int j = 0; j < M; j++) {
        sum += nu(i, j) * k[j];
      }
      dstate[i] = sum;
    }

    *dtime = tau;
  }

  void step_single(
      const NumericVector& state,
      const NumericVector& transition_rates,
      const IntegerVector& nu_row,
      const IntegerVector& nu_effect,
      double* dtime,
      NumericVector& dstate
  ) {
    int M = transition_rates.size();

    int k;

    for (int j = 0; j < M; j++) {
      k = R::rpois(transition_rates[j] * tau);
      dstate[nu_row[j]] = nu_effect[j] * k;
    }

    *dtime = tau;
  }
} ;

// [[Rcpp::export]]
SEXP make_ssa_etl(double tau) {
  SSA_etl *ssa = new SSA_etl(tau);
  XPtr<SSA_etl> ptr(ssa);
  return ptr;
}

#endif


