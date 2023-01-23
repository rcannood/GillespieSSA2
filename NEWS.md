# GillespieSSA2 0.3.0

* MINOR CHANGE: Add `debug` parameter to print out propensity functions before compiling.

* MINOR CHANGE: Add spaces between terms in the reaction propensity C++ code to 
  avoid pasting together important key words (e.g. `if` and `else`).

# GillespieSSA2 0.2.10

* MINOR CHANGE: Turn array of propensity functions into vector of propensity functions.

# GillespieSSA2 0.2.9

* MINOR CHANGE: Update RcppExports.

# GillespieSSA2 0.2.8

* BUG FIX: RNG now uses RNGScope to properly get and put the RNG state before calling RNG functions (fixes #8, thanks @bbolker!).

* DOCUMENTATION: Added example for `compile_reactions()`.

# GillespieSSA2 0.2.7 (14-07-2020)

* BUG FIX: Renamed `autoplot.ssa()` to `plot_ssa()` for compatibility with R 4.0.

# GillespieSSA2 0.2.6 (14-03-2020)

* BUG FIX: Zap small negative propensity and state values.

* FEATURE `autoplot.ssa()`: Allow plotting the firings.

* DEBUG FEATURE `ssa()`: Allow returning the GillespieSSA2 instead of running it.

# GillespieSSA2 0.2.5 (21-08-2019)

* BUG FIX: Use `fabs()` instead of `abs()` to calculate the absolute value of a 
  floating point value.
  
* BUG FIX: Precompiling returns a list of compiled function pointers, instead of 
  a single function pointer that is secretly an array of function pointers.

# GillespieSSA2 0.2.4 (05-08-2019)

GillespieSSA2 is now on CRAN!

* MAJOR CHANGE: Split up Rcpp code to make separate parts easier to test.

* TESTING: Write unit tests for many of the functions.

* MINOR CHANGE: Renamed `ssa_direct()` to `ssa_exact()`.

* MINOR CHANGE: Store firings, buffer and propensity objects only if requested.

* BUG FIX: `limiting` variable in `ssa_btl()` should be an integer, not a double.

* MINOR CHANGE: Timer now has millisecond accuracy, instead of second.

# GillespieSSA2 0.2.3 (17-07-2019)

* MAJOR CHANGE: Remove `nu` and `propensity_functions` from `ssa()`, instead
  expect a list of `reaction()` objects. This function provides a much more
  natural interface to specifying the effect and propensity of a reaction.

* MINOR CHANGE: Apply small allocation optimisations to `ssa_btl`, `ssa_etl` and `ode_em`.

# GillespieSSA2 0.2.2 (12-07-2019)

* MINOR CHANGE: Renamed `ssa_em()` to `ode_em()`.

# GillespieSSA2 0.2.1 (04-07-2019)

* BUG FIX: Fix isinf scope issue for Windows users.

* MINOR CHANGE: Move ggplot2 to Suggests.

# GillespieSSA2 0.2.0 (21-06-2019)

Complete rewrite of the package:

* The main SSA function and all SSA methods have been implemented in Rcpp.

* User-defined propensity functions get compiled to Rcpp at runtime.

* All SSA methods now assume the state-change matrix `nu` to be sparse.

# GillespieSSA2 0.1.1 (05-01-2018)

* MINOR CHANGE: Added automated testing by travis.

* MINOR CHANGE: Fixes to documentation.

* MINOR CHANGE: Also output propensities.

# GillespieSSA2 0.1.0 (31-08-2017)

Initial beta release of GillespieSSA2:

* Major restructuring of GillespieSSA code.
* Optimise main algorithm code (e.g. do not save state variables directly to the global environment...).
* Implement SSA methods in Rcpp.
* Added `ode_em()`, an Euler-Marumaya ODE method.
