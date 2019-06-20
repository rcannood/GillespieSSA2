#include <Rcpp.h>
#include "ssa.h"

using namespace Rcpp;

class SSA_BTL : public SSA {
public:
  SSA_BTL(double f_) : SSA("BTL"), f(f_) {}

  double f;

  // preallocated data structures
  IntegerVector k;
  void allocate(const int M, const int N) {
    k = IntegerVector(M);
  }

  void step_matrix(
      const NumericVector& state,
      const NumericVector& propensity,
      const IntegerMatrix& nu,
      double* dtime,
      NumericVector& dstate
  ) {
    // Fetch commonly used values
    int M = nu.ncol();
    int N = nu.nrow();

    // Calculate tau
    double tau = f / sum(propensity);
    if (tau > 1.0) tau = 1.0; // tau cannot be larger than unity

    bool coercing = false;

    // determine reaction firing
    // if propensity is zero, the reaction can't fire
    // if the effect is positive, use standard ETL
    // else use binomial distribution to determine firings
    double limiting, limiting_test, prob;
    for (int j = 0; j < M; j++) {

      if (propensity[j] > 0) {
        limiting = -1;
        for (int i = 0; i < N; i++) {
          if (nu(i, j) < 0) {
            limiting_test = (state[i] + dstate[i]) / -nu(i, j);
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
          k[j] = R::rbinom(limiting, prob);
        } else {
          k[j] = R::rpois(propensity[j] * tau);
        }
      } else {
        k[j] = 0;
      }
    }

    // determine firing effect
    double sum;
    for (int i = 0; i < N; i++) {
      sum = 0.0;
      for (int j = 0; j < M; j++) {
        sum += k[j] * nu(i, j);
      }
      dstate[i] = sum;
    }

    if (coercing) warning("coerced p to unity - consider lowering f");

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

    // Calculate tau
    double tau = f / sum(propensity);
    if (tau > 1.0) tau = 1.0; // tau cannot be larger than unity

    bool coercing = false;

    // Loop over all reaction channels having propensity fun>0
    double limiting, prob;

    for (int j = 0; j < M; j++) {
      int i = nu_row[j];
      // determine reaction firing
      // if propensity is zero, the reaction can't fire
      // if the effect is positive, use standard ETL
      // else use binomial distribution to determine firings
      if (propensity[j] <= 0) {
        k[j] = 0;
      } else if (nu_effect[j] >= 0) {
        k[j] = R::rpois(propensity[j] * tau);
      } else {
        limiting = (state[i] + dstate[i]) / -nu_effect[j];
        prob = propensity[j] * tau / limiting;
        if (prob > 1) {
          coercing = true;
          prob = 1;
        }
        k[j] = R::rbinom(limiting, prob);
      }

      // determine firing effect
      dstate[i] = k[j] * nu_effect[j];
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
