# fastgssa 0.2.2 (12-07-2019)

* MINOR CHANGE: Renamed `ssa_em()` to `ode_em()`.

# fastgssa 0.2.1 (04-07-2019)

* BUG FIX: Fix isinf scope issue for Windows users.
* MINOR CHANGE: Move ggplot2 to Suggests.

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
* Added `ode_em()`, an Euler-Marumaya ODE method.
