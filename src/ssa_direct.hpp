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
    int j = weighted_sample(transition_rates);

    for (int i = 0; i < dstate.size(); i++) {
      dstate(i) = nu(i, j);
    }

    *dtime = -log(runif(1, 0, 1)(0)) / sum(transition_rates);
  }
} ;

// [[Rcpp::export]]
SEXP make_ssa_direct() {
  SSA_direct *ssa = new SSA_direct();
  XPtr<SSA_direct> ptr(ssa);
  return ptr;
}

#endif


