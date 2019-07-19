#include <Rcpp.h>
#include "ssa_method.h"

using namespace Rcpp;

class SSA_direct : public SSA_method {
public:
  SSA_direct() : SSA_method("direct") {}

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
    // perform weighted sampling
    double max = sum(propensity);
    double ran = R::runif(0, max);
    int j = 0;
    while (ran > propensity[j]) {
      ran -= propensity[j];
      j++;
    }

    firings[j]++;

    // perform a single reaction
    for (int i = nu_p[j]; i < nu_p[j+1]; i++) {
      dstate[nu_i[i]] = nu_x[i];
    }

    *dtime = -log(R::runif(0, 1)) / sum(propensity);
  }
} ;

// [[Rcpp::export]]
SEXP make_ssa_direct() {
  SSA_direct *ssa = new SSA_direct();
  XPtr<SSA_direct> ptr(ssa);
  return ptr;
}
