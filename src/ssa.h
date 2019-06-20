#ifndef DYNGEN_SSA_H
#define DYNGEN_SSA_H

#include <Rcpp.h>

using namespace Rcpp;

class SSA {
public:
  SSA(std::string name_) : name(name_) {}
  std::string name;

  virtual ~SSA() {}

  virtual void allocate(const int M, const int N);

  virtual void step_matrix(
      const NumericVector& state,
      const NumericVector& propensity,
      const IntegerMatrix& nu,
      double* dtime,
      NumericVector& dstate
  );

  virtual void step_vector(
      const NumericVector& state,
      const NumericVector& propensity,
      const IntegerVector& nu_row,
      const IntegerVector& nu_effect,
      double* dtime,
      NumericVector& dstate
  );
} ;


#endif
