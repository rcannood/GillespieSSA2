#include <Rcpp.h>
#include "ssa.h"

using namespace Rcpp;

void SSA::allocate(const int M, const int N) {

}

void SSA::step(
    const NumericVector& state,
    const NumericVector& propensity,
    const IntegerVector& nu_i,
    const IntegerVector& nu_p,
    const IntegerVector& nu_x,
    double* dtime,
    NumericVector& dstate
) {
  stop("step() should have been overridden but wasn't!");
}
