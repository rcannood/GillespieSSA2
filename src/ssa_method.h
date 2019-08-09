#pragma once

#include <Rcpp.h>

using namespace Rcpp;

class SSA_method {
public:
  SSA_method(std::string name_) : name(name_) {}
  std::string name;

  virtual ~SSA_method() {}

  virtual void step(
      const NumericVector& state,
      const NumericVector& propensity,
      const IntegerVector& nu_i,
      const IntegerVector& nu_p,
      const IntegerVector& nu_x,
      double* dtime,
      NumericVector& dstate,
      NumericVector& firings
  );
} ;
