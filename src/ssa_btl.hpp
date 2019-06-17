#ifndef DYNGEN_SSA_BTL_H
#define DYNGEN_SSA_BTL_H

#include <Rcpp.h>
#include "ssa.hpp"
#include "utils.hpp"

using namespace Rcpp;

class SSA_BTL : public SSA {
public:
  SSA_BTL(double f_) : SSA("BTL"), f(f_) {}

  double f;

  void step(
      const NumericVector& state,
      const NumericVector& transition_rates,
      const IntegerMatrix& nu,
      double* dtime,
      NumericVector& dstate
  ) {

    // Calculate tau
    double tau = f / sum(transition_rates);
    if (tau > 1.0) tau = 1.0; // tau cannot be larger than unity

    bool coercing = false;

    for (int i = 0; i < dstate.size(); i++) {
      dstate[i] = 0.0;
    }

    // Loop over all reaction channels having propensity fun>0
    for (int j = 0; j < transition_rates.size(); j++) {

      if (transition_rates[j] > 0) {
        double k;
        double limiting = -1;
        double calc;
        for (int i = 0; i < dstate.size(); i++) {
          if (nu(i, j) < 0) {
            calc = (state[i] + dstate[i]) / -nu(i, j);
            if (limiting == -1 || calc < limiting) {
              limiting = calc;
            }
          }
        }
        if (limiting != -1) {
          double prob = transition_rates[j] * tau / limiting;
          if (prob > 1) {
            coercing = true;
            prob = 1;
          }
          k = rbinom(1, limiting, prob)[0];
        } else {
          k = rpois(1, transition_rates[j] * tau)[0];
        }

        // dstate += k * nu(_, j);
        for (int i = 0; i < dstate.size(); i++) {
          dstate[i] += k * nu(i, j);
        }
      }
    }

    if (coercing) warning("coerced p to unity - consider lowering f");

    *dtime = tau;

  }
} ;

// [[Rcpp::export]]
SEXP make_ssa_btl(double f) {
  SSA_BTL *ssa = new SSA_BTL(f);
  XPtr<SSA_BTL> ptr(ssa);
  return ptr;
}

#endif


