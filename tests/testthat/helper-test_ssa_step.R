test_ssa_step <- function(
  ssa_method,
  state,
  propensity,
  nu
) {
  sim <- new(SSA_simulation)
  sim$initialise(
    0, # num functions
    NULL, # pointer
    numeric(), # params
    0, # buffer size
    state, # initial state
    nu@i,
    nu@p,
    nu@x,
    0, # census_interval
    FALSE, # log_propensity
    FALSE, # log_firings
    FALSE, # log_buffer
    FALSE, # stop_on_neg_state
    0, # final_time
    0, # max_walltime
    NA_character_, # sim_name
    FALSE, # verbose
    100 # console_interval
  )
  sim$reset()
  sim$propensity <- propensity
  configure_method(ssa_method, sim)

  step_fun <- function() {
    sim$state <- state

    sim$make_step()

    list(
      firings = sim$dfirings[],
      dstate = sim$dstate[],
      dtime = sim$dtime[]
    )
  }

  list(
    sim = sim,
    step_fun = step_fun
  )
}
