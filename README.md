
<!-- README.md is generated from README.Rmd. Please edit that file -->

<a href="https://travis-ci.org/dynverse/gillespie"><img src="https://travis-ci.org/dynverse/gillespie.svg" align="left"></a>
<a href="https://codecov.io/gh/dynverse/gillespie"> [![AppVeyor Build
Status](https://ci.appveyor.com/api/projects/status/github/dynverse/gillespie?branch=master&svg=true)](https://ci.appveyor.com/project/dynverse/gillespie)

# `gillespie`: Gillespie’s Stochastic Simulation Algorithm for impatient people.

**gillespie** is a fast, scalable, and versatile framework for
simulating large systems with Gillespie’s Stochastic Simulation
Algorithm (SSA). It is conceptually based on
[GillespieSSA](https://cran.r-project.org/web/packages/GillespieSSA/index.html),
but rewritten entirely in Rcpp with large scale systems in mind to make
it blazingly fast. The SSA methods currently implemented are: Direct,
Explicit tau-leaping (ETL), and the Binomial tau-leaping (BTL).

## Install

You can install the development version of gillespie from GitHub with

``` r
devtools::install_github("dynverse/gillespie", build_vignettes = TRUE)
```

## Examples

The following example models are available:

  - [Decaying-Dimerization Reaction Set (Gillespie,
    2001)](vignettes/decaying_dimer.md): `vignette("decaying_dimer",
    package="gillespie")`
  - [SIRS metapopulation model (Pineda-Krch,
    2008)](vignettes/epi_chain.md): `vignette("epi_chain",
    package="gillespie")`
  - [Linear Chain System (Cao et al., 2004)](vignettes/linear_chain.md):
    `vignette("linear_chain", package="gillespie")`
  - [Pearl-Verhulst Logistic growth model (Kot,
    2001)](vignettes/logistic_growth.md): `vignette("logistic_growth",
    package="gillespie")`
  - [Lotka predator-prey model (Gillespie, 1977; Kot,
    2001)](vignettes/lotka_predator_prey.md):
    `vignette("lotka_predator_prey", package="gillespie")`
  - [Preparing your first SSA run with
    gillespie](vignettes/preparing_a_run.md):
    `vignette("preparing_a_run", package="gillespie")`
  - [Radioactive decay model (Gillespie,
    1977)](vignettes/radioactive_decay.md):
    `vignette("radioactive_decay", package="gillespie")`
  - [Rosenzweig-MacArthur predator-prey model (Pineda-Krch et al.,
    2007)](vignettes/rm_predator_prey.md): `vignette("rm_predator_prey",
    package="gillespie")`
  - [Kermack-McKendrick SIR model (Brown & Rothery,
    1993)](vignettes/sir.md): `vignette("sir", package="gillespie")`

## Latest changes

Check out `news(package = "gillespie")` or [NEWS.md](inst/NEWS.md) for a
full list of
changes.

<!-- This section gets automatically generated from inst/NEWS.md, and also generates inst/NEWS -->

### Recent changes in gillespie 0.2.3 (17-07-2019)

  - MAJOR CHANGE: Remove `nu` and `propensity_functions` from `ssa()`,
    instead expect a list of `reaction()` objects. This function
    provides a much more natural interface to specifying the effect and
    propensity of a reaction.

  - MINOR CHANGE: Apply small allocation optimisations to `ssa_btl`,
    `ssa_etl` and `ode_em`.

### Recent changes in gillespie 0.2.2 (12-07-2019)

  - MINOR CHANGE: Renamed `ssa_em()` to `ode_em()`.
