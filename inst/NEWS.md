# gillespie 0.2.3 (unreleased)

* MAJOR CHANGE: Remove `nu` and `propensity_functions` from `ssa()`, instead
  expect a list of `reaction()` objects. This function provides a much more
  natural interface to specifying the effect and propensity of a reaction.

* MINOR CHANGE: Apply small allocation optimisations to `ssa_btl`, `ssa_etl` and `ode_em`.

# gillespie 0.2.2 (12-07-2019)

* MINOR CHANGE: Renamed `ssa_em()` to `ode_em()`.

# gillespie 0.2.1 (04-07-2019)

* BUG FIX: Fix isinf scope issue for Windows users.
* MINOR CHANGE: Move ggplot2 to Suggests.

# gillespie 0.2.0 (21-06-2019)

Complete rewrite of the package:

* The main SSA function and all SSA methods have been implemented in Rcpp.
* User-defined propensity functions get compiled to Rcpp at runtime.
* All SSA methods now assume the state-change matrix `nu` to be sparse.

# gillespie 0.1.1 (05-01-2018)

* MINOR CHANGE: Added automated testing by travis.
* MINOR CHANGE: Fixes to documentation.
* MINOR CHANGE: Also output propensities.

# gillespie 0.1.0 (31-08-2017)

Initial beta release of gillespie:

* Major restructuring of GillespieSSA code.
* Optimise main algorithm code (e.g. do not save state variables directly to the global environment...).
* Implement SSA methods in Rcpp.
* Added `ode_em()`, an Euler-Marumaya ODE method.
