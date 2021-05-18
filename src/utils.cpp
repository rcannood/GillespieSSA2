#include <Rcpp.h>

using namespace Rcpp;

#include "utils.h"

int gillespie::weighted_sample(const NumericVector& weight) {
  RNGScope rngScope;

  double max = sum(weight);
  double ran = R::runif(0, max);
  int j = 0;
  while (ran > weight[j]) {
    ran -= weight[j];
    j++;
  }
  return j;
}
