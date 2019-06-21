#include <Rcpp.h>
#include "ssa.h"

using namespace Rcpp;

class SSA_BTL : public SSA {
public:
  SSA_BTL(double f_) : SSA("BTL"), f(f_) {}

  double f;

  void step(
      const NumericVector& state,
      const NumericVector& propensity,
      const IntegerVector& nu_i,
      const IntegerVector& nu_p,
      const IntegerVector& nu_x,
      double* dtime,
      NumericVector& dstate
  ) {
    int M = propensity.size();

    // Calculate tau
    double tau = f / sum(propensity);
    if (tau > 1.0) tau = 1.0; // tau cannot be larger than unity

    bool coercing = false;

    // Loop over all reaction channels having propensity fun>0
    double limiting, limiting_test, prob;
    int k;

    int nu_pos;
    double nu_val;
    for (int j = 0; j < M; j++) {
      if (propensity[j] > 0) {
        limiting = -1;
        for (int i = nu_p[j]; i < nu_p[j+1]; i++) {
          nu_val = nu_x[i]
          if (nu_val < 0) {
            nu_pos = nu_i[i]
            limiting_test = (state[nu_pos] + dstate[nu_pos]) / -nu_val;
            if (limiting == -1 || limiting_test < limiting) {
              limiting = limiting_test;
            }
          }
        }
        if (limiting != -1) {
          prob = propensity[j] * tau / limiting;
          if (prob > 1) {
            coercing = true;
            prob = 1;
          }
          k = R::rbinom(limiting, prob);
        } else {
          k = R::rpois(propensity[j] * tau);
        }

        // determine firing effect
        for (int i = nu_p[j]; i < nu_p[j+1]; i++) {
          dstate[nu_i[i]] += nu_x[i] * k;
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
