#pragma once

#include <Rcpp.h>
using namespace Rcpp;

namespace gillespie{

int weighted_sample(const NumericVector& weight);

} // namespace gillespie
