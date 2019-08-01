context("ssa simulation")

params <- c(
  "a" = 1,
  "b" = 2,
  "c" = 3,
  "d" = 4.5321
)
state <- c(
  "a_very_long_state_name_value_is_just_to_try_and_see_if_it_works" = 14,
  "short_one" = 7,
  "x" = 1,
  "y" = 0,
  "z" = 0
)
reactions <- list(
  reaction("a", c(x = +1), "react1"),
  reaction("b", c(x = +2), "react2"),
  reaction("c", c(y = +3), "react3"),
  reaction("d", c(z = +4), "react4"),
  reaction("a_very_long_state_name_value_is_just_to_try_and_see_if_it_works", c(x = +5, short_one = +1), "react5"),
  reaction("short_one", c(short_one = -1, z = 1), "react6"),
  reaction("x", c(x = -1), "react7"),
  reaction("y", c(y = -1), "react8"),
  reaction("z", c(z = -1), "react9"),
  reaction(
    propensity = "a * b * a_very_long_state_name_value_is_just_to_try_and_see_if_it_works * short_one / c / d",
    effect = c(x = 3, y = 1),
    name = "react10"
  ),
  reaction(
    propensity = "buf1 = a + 1; buf2 = buf1 + 1; buf3 = buf2 + 1; buf4 = buf3 + 1; buf4 + 1",
    effect = c(y = +1, z = +1),
    name = "react11"
  ),
  reaction(
    propensity = ~ (a * b + c) / d,
    effect = c(y = +1, z = +1),
    name = "react12"
  )
)

expected_prop <- with(as.list(c(params, state)), c(
  params,
  state,
  a * b * a_very_long_state_name_value_is_just_to_try_and_see_if_it_works * short_one / c / d,
  a + 5,
  (a * b + c) / d
)) %>%
  unname()

comp_reac <- compile_reactions(
  reactions = reactions,
  state_ids = names(state),
  params = params,
  hardcode_params = FALSE
)

test_that("perform manual steps", {
  sim <- create_simulation(
    compiled_reactions = comp_reac,
    params = params,
    method_ptr = ssa_exact()$factory(),
    initial_state = state,
    final_time = 1,
    log_buffer = TRUE,
    log_propensity = TRUE,
    log_firings = TRUE,
    census_interval = .01,
    max_walltime = 10,
    stop_on_neg_state = FALSE,
    sim_name = "test",
    verbose = TRUE,
    console_interval = 1
  )
  sim$reset()

  # check variables before simulation
  expect_equal(sim$all_zero_propensity, FALSE)
  expect_equal(sim$all_zero_state, FALSE)
  expect_equal(sim$buffer, 2:5)
  expect_equal(sim$census_interval, 0.01)
  expect_equal(sim$console_interval, 1)
  expect_equal(sim$dfirings, rep(0, length(reactions)))
  expect_equal(sim$dstate, rep(0, length(state)))
  expect_equal(sim$dtime, 0)
  expect_equal(sim$dtime_mean, 0)
  expect_equal(sim$dtime_sd, 0)
  expect_equal(sim$final_time, 1)
  expect_equal(sim$firings, rep(0, length(reactions)))
  expect_equal(sim$firings_mean, 0)
  expect_equal(sim$firings_sd, 0)
  expect_equal(sim$initial_state, state)
  expect_equal(sim$log_buffer, TRUE)
  expect_equal(sim$log_firings, TRUE)
  expect_equal(sim$log_propensity, TRUE)
  expect_equal(sim$max_walltime, 10)
  expect_equal(sim$negative_propensity, FALSE)
  expect_equal(sim$negative_state, FALSE)
  expect_equal(sim$nu_i, comp_reac$state_change@i)
  expect_equal(sim$nu_p, comp_reac$state_change@p)
  expect_equal(sim$nu_x, comp_reac$state_change@x)
  expect_equal(sim$num_functions, comp_reac$num_functions)
  expect_equal(sim$num_steps, 0)
  expect_equal(nrow(sim$output_buffer), 10)
  expect_equal(ncol(sim$output_buffer), 4)
  expect_equal(sim$output_buffer[1,], 2:5)
  expect_true(all(sim$output_buffer[-1,] == 0))
  expect_true(all(sim$output_firings == 0))
  expect_equal(sim$output_nexti, 1)
  expect_equal(nrow(sim$output_propensity), 10)
  expect_equal(ncol(sim$output_propensity), length(reactions))
  expect_equal(sim$output_propensity[1, ], expected_prop)
  expect_true(all(sim$output_propensity[-1, ] == 0))
  expect_equal(nrow(sim$output_state), 10)
  expect_equal(ncol(sim$output_state), length(state))
  expect_equal(sim$output_state[1, ], state %>% unname)
  expect_true(all(sim$output_state[-1, ] == 0))
  expect_true(all(sim$output_time == 0))
  expect_equal(sim$params, params)
  expect_equal(sim$propensity, expected_prop)
  expect_equal(sim$sim_name, "test")
  expect_equal(sim$sim_time, 0)
  expect_equal(sim$state, state %>% unname)
  expect_equal(sim$stop_on_neg_state, FALSE)
  expect_equal(sim$verbose, TRUE)

  sim$make_step()

  expect_length(sim$dfirings, length(reactions))
  expect_equal(sum(sim$dfirings == 1), 1L)
  expect_equal(sum(sim$dfirings == 0), length(reactions) - 1)

  expect_length(sim$dstate, length(state))
  expect_is(sim$dstate, "numeric")
  ix <- which(sim$dfirings == 1)
  expect_equal(sim$dstate, comp_reac$state_change[, ix])

  expect_length(sim$dtime, 1)
  expect_is(sim$dtime, "numeric")
  expect_gt(sim$dtime, 0)
  expect_lte(sim$dtime, sum(expected_prop))

  expect_equal(sim$dtime_mean, sim$dtime)
  expect_equal(sim$firings_mean, 1)

  sim$calculate_propensity()

  # try to find out whether code might cause segfaults to occur
  while (
    sim$num_steps < 1000 &&
    !sim$negative_propensity &&
    !sim$negative_state &&
    !sim$all_zero_propensity &&
    !sim$all_zero_state
  ) {
    sim$make_step()
    sim$calculate_propensity()
  }

  expect_false(sim$negative_propensity)
  expect_false(sim$negative_state)
  expect_false(sim$all_zero_propensity)
  expect_false(sim$all_zero_state)
})




test_that("perform simulation with Rcpp function", {
  sim <- create_simulation(
    compiled_reactions = comp_reac,
    params = params,
    method_ptr = ssa_exact()$factory(),
    initial_state = state,
    final_time = 1,
    log_buffer = TRUE,
    log_propensity = TRUE,
    log_firings = TRUE,
    census_interval = .01,
    max_walltime = 10,
    stop_on_neg_state = FALSE,
    sim_name = "test",
    verbose = TRUE,
    console_interval = 1
  )
  sim$reset()

  # let simulator run freely
  expect_output(
    {
      sim$run()
    },
    regexp = "Running SSA exact"
  )


  # check simulator values again
  expect_equal(sim$all_zero_propensity, FALSE)
  expect_equal(sim$all_zero_state, FALSE)
  expect_equal(sim$buffer, 2:5)
  expect_equal(sim$census_interval, 0.01)
  expect_equal(sim$console_interval, 1)
  expect_equal(sum(sim$dfirings), 1)
  expect_equal(sum(sim$dfirings != 0), 1)
  ix <- which(sim$dfirings != 0)
  exp_dstate <- comp_reac$state_change[,ix]
  expect_equal(sim$dstate, exp_dstate)
  expect_gt(sim$dtime, 0)
  expect_gt(sim$dtime_mean, 0)
  expect_gt(sim$dtime_sd, 0)
  expect_gt(sim$final_time, .99)
  expect_equal(sim$firings, rep(0, length(reactions)))
  expect_equal(sim$firings_mean, 1)
  expect_equal(sim$firings_sd, 0)
  expect_equal(sim$initial_state, state)
  expect_equal(sim$log_buffer, TRUE)
  expect_equal(sim$log_firings, TRUE)
  expect_equal(sim$log_propensity, TRUE)
  expect_equal(sim$max_walltime, 10)
  expect_equal(sim$negative_propensity, FALSE)
  expect_equal(sim$negative_state, FALSE)
  expect_equal(sim$nu_i, comp_reac$state_change@i)
  expect_equal(sim$nu_p, comp_reac$state_change@p)
  expect_equal(sim$nu_x, comp_reac$state_change@x)
  expect_equal(sim$num_functions, comp_reac$num_functions)
  expect_gte(sim$num_steps, 0)
  expect_gte(sim$output_nexti, 1)
  expect_equal(nrow(sim$output_buffer), sim$output_nexti)
  expect_equal(ncol(sim$output_buffer), 4)
  expect_true(all(apply(sim$output_buffer, 1, function(x) all(x == 2:5))))
  expect_equal(nrow(sim$output_firings), sim$output_nexti)
  expect_equal(ncol(sim$output_firings), length(reactions))
  expect_true(all(rowSums(sim$output_firings[-1,]) >= 1))
  expect_equal(nrow(sim$output_propensity), sim$output_nexti)
  expect_equal(ncol(sim$output_propensity), length(reactions))
  expect_equal(sim$output_propensity[1, ], expected_prop)
  expect_equal(nrow(sim$output_state), sim$output_nexti)
  expect_equal(ncol(sim$output_state), length(state))
  expect_equal(sim$output_state[1, ], state %>% unname)
  expect_equal(length(sim$output_time), sim$output_nexti)
  expect_true(all(diff(sim$output_time) > 0))
  expect_equal(sim$params, params)

  exp_state <- state + (comp_reac$state_change %*% colSums(sim$output_firings))[,1]

  exp_prop <- with(as.list(c(params, exp_state)), c(
    params,
    exp_state,
    a * b * a_very_long_state_name_value_is_just_to_try_and_see_if_it_works * short_one / c / d,
    a + 5,
    (a * b + c) / d
  )) %>%
    unname()
  expect_equal(sim$propensity, exp_prop)
  expect_equal(sim$sim_name, "test")
  expect_gte(sim$sim_time, 1)
  expect_equal(sim$state, exp_state %>% unname)
  expect_equal(sim$stop_on_neg_state, FALSE)
  expect_equal(sim$verbose, TRUE)

  stats <- sim$get_statistics()
  expect_is(stats, "data.frame")
})

test_that("perform simulation with R function", {
  expect_output(
    {
      out <- ssa(
        initial_state = state,
        reactions = reactions,
        final_time = 1,
        params = params,
        method = ssa_exact(),
        census_interval = .01,
        max_walltime = 10,
        stop_on_neg_state = TRUE,
        log_propensity = TRUE,
        log_firings = TRUE,
        log_buffer = TRUE,
        verbose = TRUE,
        console_interval = 1,
        sim_name = "test"
      )
    },
    regexp = "Running SSA exact"
  )

  # check simulator values again
  expect_equal(ncol(out$buffer), 4)
  expect_true(all(apply(out$buffer, 1, function(x) all(x == 2:5))))
  expect_equal(nrow(out$firings), nrow(out$buffer))
  expect_equal(ncol(out$firings), length(reactions))
  expect_true(all(rowSums(out$firings[-1,]) >= 1))
  expect_equal(nrow(out$propensity), nrow(out$buffer))
  expect_equal(ncol(out$propensity), length(reactions))
  expect_equal(out$propensity[1, ], expected_prop %>% setNames(comp_reac$reaction_ids))
  expect_equal(nrow(out$state), nrow(out$buffer))
  expect_equal(ncol(out$state), length(state))
  expect_equal(out$state[1, ], state)
  expect_equal(length(out$time), nrow(out$buffer))
  expect_true(all(diff(out$time) > 0))

  exp_state <- state + (comp_reac$state_change %*% colSums(out$firings))[,1]

  exp_prop <- with(as.list(c(params, exp_state)), c(
    params,
    exp_state,
    a * b * a_very_long_state_name_value_is_just_to_try_and_see_if_it_works * short_one / c / d,
    a + 5,
    (a * b + c) / d
  )) %>% setNames(comp_reac$reaction_ids)
  expect_equal(out$propensity[nrow(out$propensity),], exp_prop)
  expect_gte(out$time[[nrow(out$buffer)]], 1)
  expect_equal(out$state[nrow(out$state), ], exp_state)

  stats <- out$stats
  expect_is(stats, "data.frame")

  g <- autoplot.ssa(
    out,
    state = TRUE,
    propensity = TRUE,
    buffer = TRUE
  )
  expect_is(g, "ggplot")


  g <- autoplot.ssa(
    out,
    state = FALSE,
    propensity = FALSE,
    buffer = TRUE,
    geom = "step"
  )
  expect_is(g, "ggplot")
})


test_that("simulation when a negative propensity is generated", {
  params <- c("a" = -1)
  state <- c("x" = 1)
  reactions <- list(
    reaction("a", c(x = +1), "react1")
  )
  comp_reac <- compile_reactions(
    reactions = reactions,
    state_ids = names(state),
    params = params
  )
  sim <- create_simulation(
    compiled_reactions = comp_reac,
    params = params,
    initial_state = state,
    final_time = 1
  )
  sim$state <- state
  sim$calculate_propensity()
  expect_equal(sim$propensity, -1)
  expect_true(sim$negative_propensity)
})

test_that("simulation when a negative state value is generated", {
  params <- c("a" = 1)
  state <- c("x" = 1)
  reactions <- list(
    reaction("a", c(x = -1000), "react1")
  )
  comp_reac <- compile_reactions(
    reactions = reactions,
    state_ids = names(state),
    params = params
  )
  sim <- create_simulation(
    compiled_reactions = comp_reac,
    params = params,
    initial_state = state,
    final_time = 1
  )
  sim$state <- state
  sim$calculate_propensity()
  sim$make_step()
  expect_equal(sim$state %>% unname, -999)
  expect_true(sim$negative_state)
})
