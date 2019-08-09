#include <Rcpp.h>
#include "ssa_method.h"

using namespace Rcpp;

void SSA_method::step(
    const NumericVector& state,
    const NumericVector& propensity,
    const IntegerVector& nu_i,
    const IntegerVector& nu_p,
    const IntegerVector& nu_x,
    double* dtime,
    NumericVector& dstate,
    NumericVector& firings
) {
  stop("step() should have been overridden but wasn't!");
}
