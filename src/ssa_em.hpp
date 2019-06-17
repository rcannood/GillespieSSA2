#ifndef DYNGEN_SSA_EM_H
#define DYNGEN_SSA_EM_H

#include <Rcpp.h>
#include "ssa.hpp"

using namespace Rcpp;

class SSA_EM : public SSA {
public:
   SSA_EM(double tau_ = .01, double noise_strength_ = 2.0) : SSA("EM"), tau(tau_), noise_strength(noise_strength_) {}

  double tau;
  double noise_strength;

  void step(
      const NumericVector& state,
      const NumericVector& transition_rates,
      const IntegerMatrix& nu,
      double* dtime,
      NumericVector& dstate
  ) {
    *dtime = tau;

    for (int i = 0; i < dstate.size(); i++) {
      double out = 0.0;
      for (int j = 0; j < transition_rates.size(); j++) {
        out += nu(i, j) * transition_rates(j) * tau;
      }
      out += sqrt(abs(state(i))) * noise_strength * rnorm(1, 0.0, tau)(0);
      dstate(i) = out;
    }
  }
} ;

// [[Rcpp::export]]
SEXP make_ssa_em(double tau, double noise_strength) {
  SSA_EM *ssa = new SSA_EM(tau, noise_strength);
  XPtr<SSA_EM> ptr(ssa);
  return ptr;
}

#endif
