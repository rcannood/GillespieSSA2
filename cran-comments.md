# GillespieSSA2 0.2.6

* BUG FIX: Zap small negative propensity and state values.

* FEATURE `autoplot.ssa()`: Allow plotting the firings.

* DEBUG FEATURE `ssa()`: Allow returning the GillespieSSA2 instead of running it.
  
## Test environments
* local Fedora 31 install, (R 3.6.2)
* ubuntu 16.04 (on travis-ci; oldrelease, release, devel)
* Mac OS X (on travis-ci; release)
* win-builder (oldrelease, release, devel)
* Windows (on appveyor; release)

## R CMD check results

0 errors | 0 warnings | 0 notes
