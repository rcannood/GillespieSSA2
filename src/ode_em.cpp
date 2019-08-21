#include <Rcpp.h>
#include "ssa_method.h"

using namespace Rcpp;

class ODE_EM : public SSA_method {
public:
  ODE_EM(
    double tau_ = .01,
    double noise_strength_ = 2.0
  ) :
  SSA_method("EM"),
  tau(tau_),
  noise_strength(noise_strength_) {}

  double tau;
  double noise_strength;

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
    int i, j;

    // update state
    for (j = 0; j < propensity.size(); j++) {
      for (i = nu_p[j]; i < nu_p[j+1]; i++) {
        dstate[nu_i[i]] += nu_x[i] * propensity[j] * tau;
      }
      firings[j] += propensity[j] * tau;
    }

    // add noise
    for (i = 0; i < state.size(); i++) {
      dstate[i] += sqrt(fabs(state[i])) * noise_strength * R::rnorm(0.0, tau);
    }

    *dtime = tau;
  }
} ;

// [[Rcpp::export]]
SEXP make_ode_em(double tau, double noise_strength) {
  ODE_EM *ssa = new ODE_EM(tau, noise_strength);
  XPtr<ODE_EM> ptr(ssa);
  return ptr;
}

