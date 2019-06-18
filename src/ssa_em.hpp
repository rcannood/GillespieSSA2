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
    int M = nu.ncol();
    int N = nu.nrow();

    for (int i = 0; i < N; i++) {
      double out = 0.0;

      // perform each reaction 'transition_rate' times
      for (int j = 0; j < M; j++) {
        out += nu(i, j) * transition_rates[j] * tau;
      }

      // add noise
      out += sqrt(abs(state[i])) * noise_strength * R::rnorm(0.0, tau);

      // save value
      dstate[i] = out;
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
    int N = state.size();

    for (int i = 0; i < N; i++) {
      dstate[i] = sqrt(abs(state[i])) * noise_strength * R::rnorm(0, tau);
    }

    for (int j = 0; j < M; j++) {
      dstate[nu_row[j]] += nu_effect[j] * transition_rates[j] * tau;
    }

    *dtime = tau;
  }
} ;

// [[Rcpp::export]]
SEXP make_ssa_em(double tau, double noise_strength) {
  SSA_EM *ssa = new SSA_EM(tau, noise_strength);
  XPtr<SSA_EM> ptr(ssa);
  return ptr;
}

#endif
