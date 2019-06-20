#include <Rcpp.h>
#include "ssa.h"

using namespace Rcpp;

class SSA_EM : public SSA {
public:
   SSA_EM(double tau_ = .01, double noise_strength_ = 2.0) : SSA("EM"), tau(tau_), noise_strength(noise_strength_) {}

  double tau;
  double noise_strength;

  void step_matrix(
      const NumericVector& state,
      const NumericVector& propensity,
      const IntegerMatrix& nu,
      double* dtime,
      NumericVector& dstate
  ) {
    int M = nu.ncol();
    int N = nu.nrow();

    for (int i = 0; i < N; i++) {
      double out = 0.0;

      // perform each reaction 'propensity' times
      for (int j = 0; j < M; j++) {
        out += nu(i, j) * propensity[j] * tau;
      }

      // add noise
      out += sqrt(abs(state[i])) * noise_strength * R::rnorm(0.0, tau);

      // save value
      dstate[i] = out;
    }

    *dtime = tau;
  }

  void step_vector(
      const NumericVector& state,
      const NumericVector& propensity,
      const IntegerVector& nu_row,
      const IntegerVector& nu_effect,
      double* dtime,
      NumericVector& dstate
  ) {
    int M = propensity.size();
    int N = state.size();

    // perform each reaction 'propensity' times
    for (int i = 0; i < N; i++) {
      dstate[i] = sqrt(abs(state[i])) * noise_strength * R::rnorm(0, tau);
    }

    // add noise
    for (int j = 0; j < M; j++) {
      dstate[nu_row[j]] += nu_effect[j] * propensity[j] * tau;
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
