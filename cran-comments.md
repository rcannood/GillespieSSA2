# GillespieSSA2 0.2.10

* MINOR CHANGE: Turn array of propensity functions into vector of propensity functions.
  
## Test environments
* local Fedora install (R release)
* ubuntu 16.04 (with Github Actions; oldrelease, release, devel)
* Mac OS X (with Github Actions; R release)
* Windows (with Github Actions; R release)
* win-builder (oldrelease, release, devel)

## R CMD check results

```
── R CMD check results ─────────────────────────^────── GillespieSSA2 0.2.10 ────
Duration: 3m 16.7s

❯ checking installed package size ... NOTE
    installed size is  6.7Mb
    sub-directories of 1Mb or more:
      doc    1.7Mb
      libs   4.4Mb

0 errors ✔ | 0 warnings ✔ | 1 note ✖

R CMD check succeeded
```


## Reverse dependencies

A reverse dependency check was run on all downstream dependencies.
(Summary at [revdep/README.md](revdep/README.md)). No new problems were found.

```
> revdepcheck::revdep_check(timeout = as.difftime(600, units = "mins"), num_workers = 30)
── CHECK ───────────────────────────────────────────────────────── 1 packages ──
✔ dyngen 1.0.3                           ── E: 1     | W: 0     | N: 0                                                 
OK: 1                                                                                                                
BROKEN: 0
Total time: 1 min
```
