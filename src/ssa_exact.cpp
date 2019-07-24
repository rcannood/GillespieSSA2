#include <Rcpp.h>
#include "ssa_method.h"
#include "utils.h"

using namespace Rcpp;

class SSA_exact : public SSA_method {
public:
  SSA_exact() : SSA_method("exact") {}

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
    int j = gillespie::weighted_sample(propensity);

    firings[j]++;

    // perform a single reaction
    for (int i = nu_p[j]; i < nu_p[j+1]; i++) {
      dstate[nu_i[i]] = nu_x[i];
    }

    *dtime = -log(R::runif(0, 1)) / sum(propensity);
  }
} ;

// [[Rcpp::export]]
SEXP make_ssa_exact() {
  SSA_exact *ssa = new SSA_exact();
  XPtr<SSA_exact> ptr(ssa);
  return ptr;
}
