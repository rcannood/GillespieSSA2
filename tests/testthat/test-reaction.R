context("testing reactions")

test_that("reaction works", {
  reac <- reaction(
    propensity = "propensity",
    effect = c(value = -10L),
    name = "my_reaction"
  )

  expect_equal(reac$propensity, "propensity")
  expect_equal(reac$effect, c(value = -10L))
  expect_equal(reac$name, "my_reaction")

  reac <- reaction(
    propensity = ~new_propensity + 10,
    effect = c(value = -10L),
    name = "my_reaction"
  )

  expect_equal(reac$propensity, "new_propensity + 10")
  expect_equal(reac$effect, c(value = -10L))
  expect_equal(reac$name, "my_reaction")

})


test_that("reaction fails gracefully", {
  expect_error(
    reaction(propensity = 10, effect = c(value = -10L)),
    regexp = "propensity"
  )

  expect_error(
    reaction(propensity = "value", effect = c(-10L)),
    regexp = "names"
  )

  expect_error(
    reaction(propensity = "value", effect = c(value = -10L), name = 10),
    regexp = "name"
  )

  expect_error(
    reaction("value", c("par name" = -10L), "react"),
    regexp = "name"
  )
})
