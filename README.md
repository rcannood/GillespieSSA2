
<!-- README.md is generated from README.Rmd. Please edit that file -->

<a href="https://travis-ci.org/dynverse/GillespieSSA2"><img src="https://travis-ci.org/dynverse/GillespieSSA2.svg" align="left"></a>
<a href="https://codecov.io/gh/dynverse/GillespieSSA2"> [![AppVeyor
Build
Status](https://ci.appveyor.com/api/projects/status/github/dynverse/GillespieSSA2?branch=master&svg=true)](https://ci.appveyor.com/project/dynverse/GillespieSSA2)

# `GillespieSSA2`: Gillespie’s Stochastic Simulation Algorithm for impatient people.

**GillespieSSA2** is a fast, scalable, and versatile framework for
simulating large systems with Gillespie’s Stochastic Simulation
Algorithm (SSA). It is conceptually based on
[GillespieSSA](https://cran.r-project.org/web/packages/GillespieSSA/index.html),
but rewritten entirely in Rcpp with large scale systems in mind to make
it blazingly fast. The SSA methods currently implemented are: Exact,
Explicit tau-leaping (ETL), and the Binomial tau-leaping (BTL).

## Install

You can install the development version of GillespieSSA2 from GitHub
with

``` r
devtools::install_github("dynverse/GillespieSSA2", build_vignettes = TRUE)
```

## Examples

The following example models are available:

  - [Decaying-Dimerization Reaction Set (Gillespie,
    2001)](vignettes/decaying_dimer.md): `vignette("decaying_dimer",
    package="GillespieSSA2")`
  - [SIRS metapopulation model (Pineda-Krch,
    2008)](vignettes/epi_chain.md): `vignette("epi_chain",
    package="GillespieSSA2")`
  - [Linear Chain System (Cao et al., 2004)](vignettes/linear_chain.md):
    `vignette("linear_chain", package="GillespieSSA2")`
  - [Pearl-Verhulst Logistic Growth model (Kot,
    2001)](vignettes/logistic_growth.md): `vignette("logistic_growth",
    package="GillespieSSA2")`
  - [Lotka Predator-Prey model (Gillespie, 1977; Kot,
    2001)](vignettes/lotka_predator_prey.md):
    `vignette("lotka_predator_prey", package="GillespieSSA2")`
  - [Preparing your first SSA run with
    GillespieSSA2](vignettes/preparing_a_run.md):
    `vignette("preparing_a_run", package="GillespieSSA2")`
  - [Radioactive Decay model (Gillespie,
    1977)](vignettes/radioactive_decay.md):
    `vignette("radioactive_decay", package="GillespieSSA2")`
  - [Rosenzweig-MacArthur Predator-Prey model (Pineda-Krch et al.,
    2007)](vignettes/rm_predator_prey.md): `vignette("rm_predator_prey",
    package="GillespieSSA2")`
  - [Kermack-McKendrick SIR model (Brown & Rothery,
    1993)](vignettes/sir.md): `vignette("sir", package="GillespieSSA2")`

## Latest changes

Check out `news(package = "GillespieSSA2")` or [NEWS.md](inst/NEWS.md)
for a full list of
changes.

<!-- This section gets automatically generated from inst/NEWS.md, and also generates inst/NEWS -->

### Recent changes in GillespieSSA2 0.2.4 (26-07-2019)

  - MAJOR CHANGE: Split up Rcpp code to make separate parts easier to
    test.

  - TESTING: Write unit tests for many of the functions.

  - MINOR CHANGE: Renamed `ssa_direct()` to `ssa_exact()`.

  - MINOR CHANGE: Store firings, buffer and propensity objects only if
    requested.

  - BUG FIX: `limiting` variable in `ssa_btl()` should be an integer,
    not a double.

  - MINOR CHANGE: Timer now has millisecond accuracy, instead of second.

### Recent changes in GillespieSSA2 0.2.3 (17-07-2019)

  - MAJOR CHANGE: Remove `nu` and `propensity_functions` from `ssa()`,
    instead expect a list of `reaction()` objects. This function
    provides a much more natural interface to specifying the effect and
    propensity of a reaction.

  - MINOR CHANGE: Apply small allocation optimisations to `ssa_btl`,
    `ssa_etl` and `ode_em`.
