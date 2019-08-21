# GillespieSSA2 0.2.5 

* BUG FIX: Use `fabs()` instead of `abs()` to calculate the absolute value of a 
  floating point value.
  
* BUG FIX: Precompiling returns a list of compiled function pointers, instead of 
  a single function pointer that is secretly an array of function pointers.
  
## Test environments
* local Fedora 30 install, (release)
* ubuntu 16.04 (on travis-ci; oldrelease, release, devel)
* Mac OS X (on travis-ci; release)
* win-builder (devel and release)
* Windows (on appveyor; release)

## R CMD check results

0 errors | 0 warnings | 0 notes
