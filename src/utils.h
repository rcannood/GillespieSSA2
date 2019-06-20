#ifndef DYNGEN_UTILS_H
#define DYNGEN_UTILS_H

#include <Rcpp.h>
using namespace Rcpp;

template <typename T>
T resize(const T& x, int n){
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

void fill_nu_vectors(
    const IntegerMatrix& nu,
    IntegerVector& nu_row,
    IntegerVector& nu_effect,
    bool* nu_single
) {
  for (int j = 0; j < nu.ncol() && *nu_single; j++) {
    for (int i = 0; i < nu.nrow(); i++) {
      if (nu(i, j) != 0) {
        if (nu_effect[j] == 0) {
          nu_effect[j] = nu(i, j);
          nu_row[j] = i;
        } else {
          *nu_single = false;
        }
      }
    }
  }
}

#endif
