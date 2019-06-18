#ifndef DYNGEN_SSA_DIRECT_H
#define DYNGEN_SSA_DIRECT_H

#include <Rcpp.h>
#include "ssa.hpp"
#include "utils.hpp"

using namespace Rcpp;

class SSA_direct : public SSA {
public:
  SSA_direct() : SSA("direct") {}

  void step(
      const NumericVector& state,
      const NumericVector& transition_rates,
      const IntegerMatrix& nu,
      double* dtime,
      NumericVector& dstate
  ) {
    // perform weighted sampling
    double sumtr = sum(transition_rates);
    double ran = R::runif(0, sumtr);
    int j = 0;
    while (ran > transition_rates[j]) {
      ran -= transition_rates[j];
      j++;
    }

    // perform a single reaction
    for (int i = 0; i < dstate.size(); i++) {
      dstate[i] = nu(i, j);
    }

    *dtime = -log(R::runif(0, 1)) / sumtr;
  }

  void step_single(
      const NumericVector& state,
      const NumericVector& transition_rates,
      const IntegerVector& nu_row,
      const IntegerVector& nu_effect,
      double* dtime,
      NumericVector& dstate
  ) {
    // perform weighted sampling
    double sumtr = sum(transition_rates);
    double ran = R::runif(0, sumtr);
    int j = 0;
    while (ran > transition_rates[j]) {
      ran -= transition_rates[j];
      j++;
    }

    // perform a single reaction
    dstate[nu_row[j]] = nu_effect[j];

    *dtime = -log(R::runif(0, 1)) / sumtr;
  }
} ;

// [[Rcpp::export]]
SEXP make_ssa_direct() {
  SSA_direct *ssa = new SSA_direct();
  XPtr<SSA_direct> ptr(ssa);
  return ptr;
}

#endif


