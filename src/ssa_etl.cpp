#include <Rcpp.h>
#include "ssa.h"

using namespace Rcpp;

class SSA_etl : public SSA {
public:
  SSA_etl(double tau_) : SSA("ETL"), tau(tau_) {}

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
  SSA_etl *ssa = new SSA_etl(tau);
  XPtr<SSA_etl> ptr(ssa);
  return ptr;
}
