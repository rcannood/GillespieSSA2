Preparing your first SSA run with gillespie
================
2019-07-17

<!-- github markdown built using 
rmarkdown::render("vignettes/preparing_a_run.Rmd", output_format = "github_document")
-->

In order to invoke `ssa()`, the stochastic model needs at least three
components:

  - the initial state vector (`initial_state`),
  - the reactions (`reactions`),
  - the final time of the simulation (`final_time`).

The initial state vector defines the population sizes in all the states
at \(t = 0\). For example, for a system with two species `prey` and
`predators` where both have an initial population size of 1000, the
initial state vector is defined as follows.

``` r
library(gillespie)
initial_state <- c(prey = 1000, predators = 1000)
```

The reactions define the change in the number of individuals that can
occur at any given timepoint during the simulation. During an
infinitesimal period of time, the reaction can occur with a probability
defined by its propensity function. For example, a system with
abovementioned species could have three reactions; one in which the prey
population grows, one in which the predator population grows by feasting
on the prey, and one in which the predator population diminishes. The
matrix could then be defined as follows.

``` r
params <- c(c1 = 10, c2 = 0.01, c3 = 10)
reactions <- list(
  #        ↓ propensity function      ↓ effects                        ↓ name for reaction
  reaction(~c1 * prey,                c(prey = +1),                    name = "prey_up"),
  reaction(~c2 * prey * predators,    c(prey = -1, predators = +1),    name = "predation"),
  reaction(~c3 * predators,           c(predators = -1),               name = "pred_down")
)
```

The simulation can be started by calling the `ssa()` function.

``` r
out <- 
  ssa(
    initial_state = initial_state,
    reactions = reactions,
    params = params,
    method = ssa_direct(),
    final_time = 5,
    census_interval = .001,
    verbose = TRUE
  )
#> Running SSA direct with console output every 1 seconds
#> Start time: CURRTIME
#> walltime: 0, simtime: 0
#> SSA finished!
```

``` r
print(out$stats)
#>   method stop_simtime stop_extinction stop_negative_state stop_zero_prop
#> 1 direct         TRUE           FALSE               FALSE          FALSE
#>   stop_walltime walltime_start walltime_end walltime_elapsed num_steps
#> 1         FALSE     1563380511   1563380511                0    151518
#>     dtime_mean     dtime_sd
#> 1 3.299954e-05 1.538995e-07
```

``` r
ssa_plot(out)
```

<img src="preparing_a_run_files/figure-gfm/unnamed-chunk-5-1.png" width="100%" />
