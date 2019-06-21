
<!-- README.md is generated from README.Rmd. Please edit that file -->

<a href="https://travis-ci.org/dynverse/fastgssa"><img src="https://travis-ci.org/dynverse/fastgssa.svg" align="left"></a>
<a href="https://codecov.io/gh/dynverse/fastgssa"> [![AppVeyor Build
Status](https://ci.appveyor.com/api/projects/status/github/dynverse/fastgssa?branch=master&svg=true)](https://ci.appveyor.com/project/dynverse/fastgssa)

# `fastgssa`: Gillespie’s Stochastic Simulation Algorithm for impatient people.

**fastgssa** is a fast, scalable, and versatile framework for simulating
large systems with Gillespie’s Stochastic Simulation Algorithm (SSA). It
is conceptually based on
[GillespieSSA](https://cran.r-project.org/web/packages/GillespieSSA/index.html),
but rewritten entirely in Rcpp with large scale systems in mind to make
it blazingly fast. The SSA methods currently implemented are: Direct,
Explicit tau-leaping (ETL), and the Binomial tau-leaping (BTL).

## Install

You can install the development version of fastgssa from GitHub with

    ```R
    devtools::install_github("dynverse/fastgssa", build_vignettes = TRUE)
    ```

## Examples

To get started, have a look at one of the vignettes:

  - [Preparing a run](vignettes/preparing_a_run.md):
    `vignette("preparing_a_run", package = "fastgssa")`.

## Latest changes

Check out `news(package = "fastgssa")` or [NEWS.md](inst/NEWS.md) for a
full list of
changes.

<!-- This section gets automatically generated from inst/NEWS.md, and also generates inst/NEWS -->

### Recent changes in fastgssa 0.2.0 (unreleased)

Complete rewrite of the package: \* The main SSA function and all SSA
methods have been implemented in Rcpp. \* User-defined propensity
functions get compiled to Rcpp at runtime. \* All SSA methods now assume
the state-change matrix `nu` to be sparse.

### Recent changes in fastgssa 0.1.1 (05-01-2018)

  - MINOR CHANGE: Added automated testing by travis.
  - MINOR CHANGE: Fixes to documentation.
  - MINOR CHANGE: Also output propensities.
