
<!-- README.md is generated from README.Rmd. Please edit that file -->
<!-- badges: start -->

[![CRAN
Status](https://www.r-pkg.org/badges/version/GillespieSSA2)](https://cran.r-project.org/package=GillespieSSA2)
[![CRAN
Downloads](https://cranlogs.r-pkg.org/badges/GillespieSSA2)](https://cran.r-project.org/package=GillespieSSA2)
![R-CMD-check](https://github.com/rcannood/GillespieSSA2/workflows/R-CMD-check/badge.svg)
[![DOI](https://img.shields.io/badge/doi-10.1101/2020.02.06.936971-green)](https://doi.org/10.1101/2020.02.06.936971)
[![Coverage
Status](https://app.codecov.io/gh/rcannood/GillespieSSA2/branch/master/graph/badge.svg)](https://app.codecov.io/gh/rcannood/GillespieSSA2?branch=master)
<!-- badges: end -->

# `GillespieSSA2`: Gillespie’s Stochastic Simulation Algorithm for impatient people.

**GillespieSSA2** is a fast, scalable, and versatile framework for
simulating large systems with Gillespie’s Stochastic Simulation
Algorithm (SSA) (Cannoodt et al. 2021). This package is the spiritual
successor to the GillespieSSA package originally written by Mario
Pineda-Krch (Pineda-Krch 2008).

GillespieSSA2 has the following added benefits:

-   The whole algorithm is run in Rcpp which results in major speed
    improvements (>100x). Even your propensity functions (reactions) are
    being compiled to Rcpp!
-   Parameters and variables have been renamed to make them easier to
    understand.
-   Many unit tests try to ensure that the code works as intended.

The SSA methods currently implemented are: Exact (`ssa_exact()`),
Explicit tau-leaping (`ssa_etl()`), and the Binomial tau-leaping
(`ssa_btl()`).

## Install

You can install:

-   the latest released version from CRAN with

    ``` r
    install.packages("GillespieSSA2")
    ```

-   the latest development version from github with

    ``` r
    devtools::install_github("rcannood/GillespieSSA2", build_vignettes = TRUE)
    ```

If you encounter a bug, please file a minimal reproducible example on
the [issues](https://github.com/rcannood/GillespieSSA2/issues) page.

## Examples

The following example models are available:

-   Introduction to GillespieSSA2:  
    `vignette("an_introduction", package="GillespieSSA2")`
-   Converting from GillespieSSA to GillespieSSA2:  
    `vignette("converting_from_GillespieSSA", package="GillespieSSA2")`
-   Decaying-Dimerization Reaction Set:  
    `vignette("decaying_dimer", package="GillespieSSA2")`
-   SIRS metapopulation model:  
    `vignette("epi_chain", package="GillespieSSA2")`
-   Linear Chain System:  
    `vignette("linear_chain", package="GillespieSSA2")`
-   Pearl-Verhulst Logistic Growth model:  
    `vignette("logistic_growth", package="GillespieSSA2")`
-   Lotka Predator-Prey model:  
    `vignette("lotka_predator_prey", package="GillespieSSA2")`
-   Radioactive Decay model:  
    `vignette("radioactive_decay", package="GillespieSSA2")`
-   Rosenzweig-MacArthur Predator-Prey model:  
    `vignette("rm_predator_prey", package="GillespieSSA2")`
-   Kermack-McKendrick SIR model:  
    `vignette("sir", package="GillespieSSA2")`

## References

<div id="refs" class="references csl-bib-body hanging-indent">

<div id="ref-Cannoodt2021" class="csl-entry">

Cannoodt, Robrecht, Wouter Saelens, Louise Deconinck, and Yvan Saeys.
2021. “Spearheading Future Omics Analyses Using Dyngen, a Multi-Modal
Simulator of Single Cells.” *Nature Communications* 12 (1).
<https://doi.org/10.1038/s41467-021-24152-2>.

</div>

<div id="ref-PinedaKrch2008" class="csl-entry">

Pineda-Krch, Mario. 2008. “GillespieSSA: Implementing the Stochastic
Simulation Algorithm in r.” *Journal of Statistical Software* 25 (12).
<https://doi.org/10.18637/jss.v025.i12>.

</div>

</div>
