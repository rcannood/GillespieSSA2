
  - [`GillespieSSA2`: Gillespie’s Stochastic Simulation Algorithm for
    impatient
    people.](#gillespiessa2-gillespies-stochastic-simulation-algorithm-for-impatient-people.)
      - [Install](#install)
      - [Examples](#examples)
      - [Latest changes](#latest-changes)
          - [Recent changes in GillespieSSA2 0.2.7
            (14-07-2020)](#recent-changes-in-gillespiessa2-0.2.7-14-07-2020)
          - [Recent changes in GillespieSSA2 0.2.6
            (14-03-2020)](#recent-changes-in-gillespiessa2-0.2.6-14-03-2020)

<!-- README.md is generated from README.Rmd. Please edit that file -->

![R-CMD-check](https://github.com/rcannood/GillespieSSA2/workflows/R-CMD-check/badge.svg)

# `GillespieSSA2`: Gillespie’s Stochastic Simulation Algorithm for impatient people.

**GillespieSSA2** is a fast, scalable, and versatile framework for
simulating large systems with Gillespie’s Stochastic Simulation
Algorithm (SSA). This package is the spiritual successor to the
GillespieSSA package originally written by Mario Pineda-Krch.

GillespieSSA2 has the following added benefits:

  - The whole algorithm is run in Rcpp which results in major speed
    improvements (\>100x). Even your propensity functions (reactions)
    are being compiled to Rcpp\!
  - Parameters and variables have been renamed to make them easier to
    understand.
  - Many unit tests try to ensure that the code works as intended.

The SSA methods currently implemented are: Exact (`ssa_exact()`),
Explicit tau-leaping (`ssa_etl()`), and the Binomial tau-leaping
(`ssa_btl()`).

## Install

You can install:

  - the latest released version from CRAN with
    
    ``` r
    install.packages("GillespieSSA2")
    ```

  - the latest development version from github with
    
    ``` r
    devtools::install_github("rcannood/GillespieSSA2", build_vignettes = TRUE)
    ```

If you encounter a bug, please file a minimal reproducible example on
the [issues](https://github.com/rcannood/GillespieSSA2/issues) page.

## Examples

The following example models are available:

  - [Introduction to GillespieSSA2](vignettes/an_introduction.md):  
    `vignette("an_introduction", package="GillespieSSA2")`
  - [Converting from GillespieSSA to
    GillespieSSA2](vignettes/converting_from_GillespieSSA.md):  
    `vignette("converting_from_GillespieSSA", package="GillespieSSA2")`
  - [Decaying-Dimerization Reaction Set (Gillespie,
    2001)](vignettes/decaying_dimer.md):  
    `vignette("decaying_dimer", package="GillespieSSA2")`
  - [SIRS metapopulation model (Pineda-Krch,
    2008)](vignettes/epi_chain.md):  
    `vignette("epi_chain", package="GillespieSSA2")`
  - [Linear Chain System (Cao et al.,
    2004)](vignettes/linear_chain.md):  
    `vignette("linear_chain", package="GillespieSSA2")`
  - [Pearl-Verhulst Logistic Growth model (Kot,
    2001)](vignettes/logistic_growth.md):  
    `vignette("logistic_growth", package="GillespieSSA2")`
  - [Lotka Predator-Prey model (Gillespie, 1977; Kot,
    2001)](vignettes/lotka_predator_prey.md):  
    `vignette("lotka_predator_prey", package="GillespieSSA2")`
  - [Radioactive Decay model (Gillespie,
    1977)](vignettes/radioactive_decay.md):  
    `vignette("radioactive_decay", package="GillespieSSA2")`
  - [Rosenzweig-MacArthur Predator-Prey model (Pineda-Krch et al.,
    2007)](vignettes/rm_predator_prey.md):  
    `vignette("rm_predator_prey", package="GillespieSSA2")`
  - [Kermack-McKendrick SIR model (Brown & Rothery,
    1993)](vignettes/sir.md):  
    `vignette("sir", package="GillespieSSA2")`

## Latest changes

Check out `news(package = "GillespieSSA2")` or [NEWS.md](NEWS.md) for a
full list of changes.

<!-- This section gets automatically generated from NEWS.md -->

### Recent changes in GillespieSSA2 0.2.7 (14-07-2020)

  - BUG FIX: Renamed `autoplot.ssa()` to `plot_ssa()` for compatibility
    with R 4.0.

### Recent changes in GillespieSSA2 0.2.6 (14-03-2020)

  - BUG FIX: Zap small negative propensity and state values.

  - FEATURE `autoplot.ssa()`: Allow plotting the firings.

  - DEBUG FEATURE `ssa()`: Allow returning the GillespieSSA2 instead of
    running it.
