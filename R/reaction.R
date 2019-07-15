reaction <- function(
  propensity,
  ...,
  effect = NULL,
  name = NULL
) {
  # check propensity
  assert_that(
    is_formula(propensity) || is.character(propensity)
  )

  # collect effects
  effect <- c(..., effect)
  assert_that(
    length(effect) > 0,
    !is.null(names(effect)),
    all(names(effect) != "")
  )

  # return output
  out <- list(
    propensity = propensity,
    effect = effect,
    name = name
  )
  class(out) <- "fastgssa_reactions"
  out
}
