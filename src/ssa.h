#pragma once

#include <Rcpp.h>

using namespace Rcpp;

class SSA {
public:
  SSA(std::string name_) : name(name_) {}
  std::string name;

  virtual ~SSA() {}

  virtual void allocate(const int M, const int N);

  virtual void step(
      const NumericVector& state,
      const NumericVector& propensity,
      const IntegerVector& nu_i,
      const IntegerVector& nu_p,
      const IntegerVector& nu_x,
      double* dtime,
      NumericVector& dstate
  );
} ;
