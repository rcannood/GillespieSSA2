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

  void step(
      const NumericVector& state,
      const NumericVector& transition_rates,
      const IntegerMatrix& nu,
      double* dtime,
      NumericVector& dstate
  ) {
    NumericVector k = no_init(transition_rates.size());
    for (int i = 0; i < transition_rates.size(); i++) {
      k[i] = rpois(1, transition_rates[i] * tau)[0];
    }

    for (int i = 0; i < dstate.size(); i++) {
      double sum = 0.0;
      for (int j = 0; j < transition_rates.size(); j++) {
        sum += nu(i, j) * k[j];
      }
      dstate[i] = sum;
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


