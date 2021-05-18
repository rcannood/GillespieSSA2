#include <Rcpp.h>
#include "ssa_method.h"

using namespace Rcpp;

class SSA_ETL : public SSA_method {
public:
  SSA_ETL(double tau_) : SSA_method("ETL"), tau(tau_) {}

  double tau;

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
    RNGScope rngScope;

    int k;
    int i, j;

    for (j = 0; j < propensity.size(); j++) {
      // determine reaction firing
      k = R::rpois(propensity[j] * tau);

      firings[j] += k;

      // determine firing effect
      for (i = nu_p[j]; i < nu_p[j+1]; i++) {
        dstate[nu_i[i]] += nu_x[i] * k;
      }
    }

    // tau leap
    *dtime = tau;
  }
} ;


// [[Rcpp::export]]
SEXP make_ssa_etl(double tau) {
  SSA_ETL *ssa = new SSA_ETL(tau);
  XPtr<SSA_ETL> ptr(ssa);
  return ptr;
}
