#include <Rcpp.h>
#include <math.h>

#include "ssa_method.h"
#include "utils.h"

using namespace Rcpp;

// [[Rcpp::export]]
List test_ssa_method_cpp(
    SEXP ssa_alg,
    NumericVector& state,
    NumericVector& propensity,
    IntegerVector& nu_i,
    IntegerVector& nu_p,
    IntegerVector& nu_x
) {
  SSA_method *ssa_alg_ = XPtr<SSA_method>(ssa_alg);
  double dtime = 0;
  NumericVector dstate(state.size());
  NumericVector firings(propensity.size());
  ssa_alg_->step(state, propensity, nu_i, nu_p, nu_x, &dtime, dstate, firings);
  return List::create(
    _["dtime"] = dtime,
    _["dstate"] = dstate,
    _["firings"] = firings
  );
}

typedef void (*TR_FUN)(const NumericVector&, const NumericVector&, const double, NumericVector&, NumericVector&);

// [[Rcpp::export]]
List test_propensity_cpp(
    List propensity_funs,
    NumericVector& params,
    int buffer_size,
    int propensity_size,
    NumericVector& state,
    double sim_time
) {
  TR_FUN *prop_funs = new TR_FUN[propensity_funs.size()];
  for (int i = 0; i < propensity_funs.size(); i++) {
    prop_funs[i] = *XPtr<TR_FUN>(as<SEXP>(propensity_funs(i)));
  }
  NumericVector buffer = NumericVector(buffer_size);
  NumericVector propensity = NumericVector(propensity_size);

  for (int i = 0; i < propensity_funs.size(); i++) {
    prop_funs[i](state, params, sim_time, propensity, buffer);
  }

  delete[] prop_funs;

  return List::create(
    _["propensity"] = propensity,
    _["buffer"] = buffer
  );
}


