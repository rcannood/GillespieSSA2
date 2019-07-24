#include <Rcpp.h>
#include "ssa_method.h"

using namespace Rcpp;

class SSA_BTL : public SSA_method {
public:
  SSA_BTL(double mean_firings_) : SSA_method("BTL"), mean_firings(mean_firings_) {}

  double mean_firings;

  void step(
      const NumericVector& state,
      const NumericVector& propensity,
      const IntegerVector& nu_i,
      const IntegerVector& nu_p,
      const IntegerVector& nu_x,
      double* dtime,
      NumericVector& dstate,
      NumericVector& firings
  ) {
    // TODO: could split up nu into positive and negative nus

    int M = propensity.size();

    // Calculate tau
    double tau = mean_firings / sum(propensity);
    if (tau > 1.0) tau = 1.0; // tau cannot be larger than one

    bool coercing = false;

    // Loop over all reaction channels having propensity fun>0
    double prob;
    int limiting, limiting_test, k;

    int nu_pos;
    double nu_val;
    int i, j;
    for (j = 0; j < M; j++) {
      if (propensity[j] > 0) {
        limiting = 1;
        for (i = nu_p[j]; i < nu_p[j+1]; i++) {
          nu_val = nu_x[i];
          if (nu_val < 0) {
            nu_pos = nu_i[i];
            limiting_test = (state[nu_pos] + dstate[nu_pos]) / -nu_val;
            if (limiting_test < limiting) {
              limiting = limiting_test;
            }
          }
        }
        if (limiting < 1) {
          prob = propensity[j] * tau / limiting;
          if (prob > 1) {
            coercing = true;
            prob = 1;
          }
          k = R::rbinom(limiting, prob);
        } else {
          k = R::rpois(propensity[j] * tau);
        }

        firings[j] += k;

        // TODO: Could determine firings first,
        // if prob > 1, mean_firings could be set to 1 for this iteration only

        // determine firing effect
        for (i = nu_p[j]; i < nu_p[j+1]; i++) {
          dstate[nu_i[i]] += nu_x[i] * k;
        }
      }
    }

    if (coercing) warning("coerced p to unity - consider lowering f");

    *dtime = tau;
  }
} ;

// [[Rcpp::export]]
SEXP make_ssa_btl(double mean_firings) {
  SSA_BTL *ssa = new SSA_BTL(mean_firings);
  XPtr<SSA_BTL> ptr(ssa);
  return ptr;
}
