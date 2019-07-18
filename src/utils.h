#pragma once

#include <Rcpp.h>
using namespace Rcpp;

namespace gillespie{

template <typename T>
T resize_vector(const T& x, int n){
  int oldsize = x.size();
  if (n < oldsize) {
    oldsize = n;
  }
  T y(n);
  for( int i = 0; i < oldsize; i++) {
    y[i] = x[i];
  }
  return y;
}

template <typename T>
T resize_rows(const T& x, int n){
  int oldsize = x.nrow();
  if (n < oldsize) {
    oldsize = n;
  }
  T y(n, x.ncol());
  for( int i = 0; i < oldsize; i++) {
    for (int j = 0; j < x.ncol(); j++) {
      y(i, j) = x(i, j);
    }
  }
  return y;
}

int weighted_sample(const NumericVector& weight) {
  double max = sum(weight);
  double ran = R::runif(0, max);
  int j = 0;
  while (ran > weight[j]) {
    ran -= weight[j];
    j++;
  }
  return j;
}

} // namespace gillespie

// #endif
