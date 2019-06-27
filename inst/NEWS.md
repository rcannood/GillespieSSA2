# fastgssa 0.2.0 (21-06-2019)

Complete rewrite of the package:

* The main SSA function and all SSA methods have been implemented in Rcpp.
* User-defined propensity functions get compiled to Rcpp at runtime.
* All SSA methods now assume the state-change matrix `nu` to be sparse.

# fastgssa 0.1.1 (05-01-2018)

* MINOR CHANGE: Added automated testing by travis.
* MINOR CHANGE: Fixes to documentation.
* MINOR CHANGE: Also output propensities.

# fastgssa 0.1.0 (31-08-2017)

Initial beta release of fastgssa:

* Major restructuring of GillespieSSA code.
* Optimise main algorithm code (e.g. do not save state variables directly to the global environment...).
* Implement SSA methods in Rcpp.
* Added `ssa_em()`, an Euler-Marumaya SDE method.
