#include <Rcpp.h>
#include "ssa.h"

using namespace Rcpp;

void SSA::allocate(const int M, const int N) {

}

void SSA::step_matrix(
    const NumericVector& state,
    const NumericVector& propensity,
    const IntegerMatrix& nu,
    double* dtime,
    NumericVector& dstate
) {
  stop("step_matrix() should have been overridden but wasn't!");
}

void SSA::step_vector(
    const NumericVector& state,
    const NumericVector& propensity,
    const IntegerVector& nu_row,
    const IntegerVector& nu_effect,
    double* dtime,
    NumericVector& dstate
) {
  stop("step_vector() should have been overridden but wasn't!");
}
