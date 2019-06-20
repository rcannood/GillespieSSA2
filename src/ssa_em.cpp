#include <Rcpp.h>
#include "ssa.h"

using namespace Rcpp;

class SSA_EM : public SSA {
public:
   SSA_EM(double tau_ = .01, double noise_strength_ = 2.0) : SSA("EM"), tau(tau_), noise_strength(noise_strength_) {}

  double tau;
  double noise_strength;

  void step(
      const NumericVector& state,
      const NumericVector& propensity,
      const IntegerVector& nu_i,
      const IntegerVector& nu_p,
      const IntegerVector& nu_x,
      double* dtime,
      NumericVector& dstate
  ) {
    // update state
    for (int j = 0; j < propensity.size(); j++) {
      for (int i = nu_p[j]; i < nu_p[j+1]; i++) {
        dstate[nu_i[i]] += nu_x[i] * propensity[j] * tau;
      }
    }

    // add noise
    for (int i = 0; i < state.size(); i++) {
      dstate[i] += sqrt(abs(state[i])) * noise_strength * R::rnorm(0.0, tau);
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
