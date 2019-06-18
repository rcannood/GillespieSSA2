#ifndef DYNGEN_SSA_H
#define DYNGEN_SSA_H

#include <Rcpp.h>
#include "utils.hpp"

using namespace Rcpp;

class SSA {
public:
  SSA(std::string name_) : name(name_) {}
  std::string name;

  virtual ~SSA() {}

  virtual void allocate(const int M, const int N) {}

  virtual void step(
      const NumericVector& state,
      const NumericVector& transition_rates,
      const IntegerMatrix& nu,
      double* dtime,
      NumericVector& dstate
  ) {
    stop("step() should have been overridden but wasn't!");
  }

  virtual void step_single(
      const NumericVector& state,
      const NumericVector& transition_rates,
      const IntegerVector& nu_row,
      const IntegerVector& nu_effect,
      double* dtime,
      NumericVector& dstate
  ) {
    stop("step_single() should have been overridden but wasn't!");
  }
} ;


#endif
