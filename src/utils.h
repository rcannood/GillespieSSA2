#pragma once

#include <Rcpp.h>
using namespace Rcpp;

namespace gillespie{

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
