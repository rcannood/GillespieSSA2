#pragma once

#include <Rcpp.h>
using namespace Rcpp;

namespace gillespie{

int weighted_sample(const NumericVector& weight);

uint64_t timems();

} // namespace gillespie
