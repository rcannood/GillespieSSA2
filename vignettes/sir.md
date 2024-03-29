Kermack-McKendrick SIR model
================

<!-- github markdown built using 
rmarkdown::render("vignettes/sir.Rmd", output_format = "github_document")
-->

The Kermack-McKendrick SIR model (Brown, Rothery, et al. 1993) is
defined as

    dS/dt = -beta*N*S
    dI/dt = beta*N*S - gamma*I
    dR/dt = gamma*I

Note that simulations of this model can generate in all zero propensity,
if the first reaction is a recovery of the single ‘Infected’ individual.

Define parameters

``` r
library(GillespieSSA2)
sim_name <- "Kermack-McKendrick SIR model"
params <- c(beta = .001, gamma = .1)
final_time <- 100
initial_state <- c(S = 500, I = 1, R = 0)
```

Define reactions

``` r
reactions <- list(
  reaction("beta * S * I", c(S = -1, I = +1), name = "transmission"),
  reaction("gamma * I", c(I = -1, R = +1), name = "recovery")
)
```

Run simulations with the Exact method

``` r
set.seed(1)
out <- ssa(
  initial_state = initial_state,
  reactions = reactions,
  params = params,
  final_time = final_time,
  method = ssa_exact(),
  sim_name = sim_name
) 
plot_ssa(out)
```

<img src="sir_files/figure-gfm/exact-1.png" width="100%" />

Run simulations with the Explict tau-leap method

``` r
set.seed(1)
out <- ssa(
  initial_state = initial_state,
  reactions = reactions,
  params = params,
  final_time = final_time,
  method = ssa_etl(),
  sim_name = sim_name
) 
plot_ssa(out)
```

<img src="sir_files/figure-gfm/etl-1.png" width="100%" />

Run simulations with the Binomial tau-leap method

``` r
set.seed(2)
out <- ssa(
  initial_state = initial_state,
  reactions = reactions,
  params = params,
  final_time = final_time,
  method = ssa_btl(),
  sim_name = sim_name
) 
plot_ssa(out)
```

<img src="sir_files/figure-gfm/btl-1.png" width="100%" />

## References

<div id="refs" class="references csl-bib-body hanging-indent">

<div id="ref-Brown1993" class="csl-entry">

Brown, David, Peter Rothery, et al. 1993. *Models in Biology:
Mathematics, Statistics and Computing.* John Wiley & Sons Ltd.

</div>

</div>
