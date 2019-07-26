#include <Rcpp.h>
#include <chrono>

using namespace Rcpp;

#include "utils.h"

int gillespie::weighted_sample(const NumericVector& weight) {
  double max = sum(weight);
  double ran = R::runif(0, max);
  int j = 0;
  while (ran > weight[j]) {
    ran -= weight[j];
    j++;
  }
  return j;
}

uint64_t gillespie::timems() {
  using namespace std::chrono;
  return duration_cast<milliseconds>(system_clock::now().time_since_epoch()).count();
}
