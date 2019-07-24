context("compile reactions")

params <- c(
  "a" = 1,
  "b" = 2,
  "c" = 3,
  "d" = 4.5321
)
state <- c(
  "a_very_long_state_name_value_is_just_to_try_and_see_if_it_works" = 12.34,
  "short_one" = 7,
  "x" = 1,
  "y" = 0,
  "z" = 0
)
reactions <- list(
  reaction("a", c(x = +1, y = +1), "react1"),
  reaction("b", c(x = +1, y = -1), "react2"),
  reaction("c", c(x = +1, z = +1), "react3"),
  reaction("d", c(x = +1, z = -1), "react4"),
  reaction("a_very_long_state_name_value_is_just_to_try_and_see_if_it_works", c(x = +1), "react5"),
  reaction("short_one", c(x = -1, y = +1), "react6"),
  reaction("x", c(x = -1, y = -1), "react7"),
  reaction("y", c(x = -1, z = +1), "react8"),
  reaction("z", c(x = -1, z = -1), "react9"),
  reaction(
    propensity = "a * b * c * d / a_very_long_state_name_value_is_just_to_try_and_see_if_it_works / short_one",
    effect = c(x = -1),
    name = "react10"
  ),
  reaction(
    propensity = "buf1 = a + 1; buf2 = buf1 + 1; buf3 = buf2 + 1; buf4 = buf3 + 1; buf4 + 1",
    effect = c(y = +1, z = +1),
    name = "react11"
  ),
  reaction(
    propensity = ~ (a * b + c) / d,
    effect = c(y = +1, z = -1),
    name = "react12"
  )
)

test_that("compilation works", {
  comp_reac <- compile_reactions(
    reactions = reactions,
    state_ids = names(state),
    params = params,
    hardcode_params = FALSE
  )

  # check state change matrix
  nu <- comp_reac$state_change
  expect_equal(nrow(nu), length(state))
  expect_equal(ncol(nu), length(reactions))

  expected_p <- reactions %>% map_int(~length(.$effect)) %>% cumsum %>% c(0, .)
  expect_equal(nu@p, expected_p)
  expected_i <- reactions %>% map(~names(.$effect)) %>% unlist() %>% match(names(state)) - 1
  expect_equal(nu@i, expected_i)
  expected_x <- reactions %>% map(~.$effect) %>% unlist() %>% unname()
  expect_equal(expected_x, nu@x)

  # check reaction ids
  expected_reac_ids <- map_chr(reactions, "name")
  expect_equal(comp_reac$reaction_ids, expected_reac_ids)

  # check other params
  expect_equal(comp_reac$buffer_ids, paste0("buf", 1:4))
  expect_equal(comp_reac$buffer_size, 4)
  expect_equal(comp_reac$hardcode_params, FALSE)

  # check funs (have to start up a simulation object for this)
  out <- test_propensity_calculation(
    comp_reac,
    params,
    state,
    0
  )
  expected_prop <- with(as.list(c(params, state)), c(
    params,
    state,
    a * b * c * d / a_very_long_state_name_value_is_just_to_try_and_see_if_it_works / short_one,
    a + 5,
    (a * b + c) / d
  )) %>%
    unname()
  expect_equal(out$propensity, expected_prop, tolerance = .001)
  expect_equal(out$buffer, 2:5, tolerance = .001)
})


test_that("compilation works when params are hardcoded", {
  comp_reac <- compile_reactions(
    reactions = reactions,
    state_ids = names(state),
    params = params,
    hardcode_params = TRUE
  )
  out <- test_propensity_calculation(
    comp_reac,
    params,
    state,
    0
  )
  expected_prop <- with(as.list(c(params, state)), c(
    params,
    state,
    a * b * c * d / a_very_long_state_name_value_is_just_to_try_and_see_if_it_works / short_one,
    a + 5,
    (a * b + c) / d
  )) %>%
    unname()
  expect_equal(out$propensity, expected_prop, tolerance = .001)
  expect_equal(out$buffer, 2:5, tolerance = .001)
})
