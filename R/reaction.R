#' Define a reaction
#'
#' During an SSA simulation, at any infinitesimal time interval,
#' a reaction will occur with a probability defined according to its
#' propensity. If it does, then it will change the state vector according
#' to its effects.
#'
#' It is possible to use 'buffer' values in order to speed up the computation
#' of the propensity functions. For instance, instead of `"(c3 * s1) / (1 + c3 * c1)"`,
#' it is possible to write `"buf = c3 * s1; buf / (buf + 1)"` instead.
#'
#' @param propensity `[character/formula]` A character or formula representation of the propensity function, written in C++.
#' @param effect `[named integer vector]` The change in state caused by this reaction.
#' @param name `[character]` A name for this reaction (Optional). May only contain characters matching `[A-Za-z0-9_]`.
#'
#' @return `[SSA_reaction]` This object describes a single reaction as part of an SSA simulation. It contains the following member values:
#'
#' * `r[["propensity"]]`: The propensity function as a character.
#' * `r[["effect"]]`: The change in state caused by this reaction.
#' * `r[["name"]]`: The name of the reaction, `NA_character_` if no name was provided.
#'
#' @importFrom rlang is_formula is_integerish
#'
#' @export
#' @examples
#' #        propensity                        effect
#' reaction(~ c1 * s1,                          c(s1 = -1))
#' reaction("c2 * s1 * s1",                     c(s1 = -2, s2 = +1))
#' reaction("buf = c3 * s1; buf / (buf + 1)",   c(s1 = +2))
reaction <- function(
  propensity,
  effect,
  name = NA_character_
) {
  # check propensity
  assert_that(
    is_formula(propensity) || is.character(propensity)
  )
  if (is_formula(propensity)) {
    propensity <- as.character(propensity)[[2]]
  }

  # check effects
  assert_that(
    is_integerish(effect),
    length(effect) > 0,
    !is.null(names(effect)),
    all(grepl("^[A-Za-z_0-9]+$", names(effect)))
  )

  effect <- set_names(as.integer(effect), names(effect))

  # check name
  assert_that(
    is_scalar_character(name),
    is.na(name) || grepl("^[A-Za-z_0-9]+$", name)
  )

  # return output
  out <- list(
    propensity = propensity,
    effect = effect,
    name = name
  )
  class(out) <- "SSA_reaction"
  out
}

#' Print various SSA objects
#' @param x An SSA reaction or SSA method
#' @param ... Not used
#' @export
#'
#' @rdname print_ssa
print.SSA_reaction <- function(x, ...) {
  effect <- x[["effect"]]
  effect_str <-
    if (length(effect) > 0) {
      paste0(names(x[["effect"]]), ": ", ifelse(x[["effect"]] > 0, "+", ""), x[["effect"]], collapse = ", ")
    } else {
      "<none>"
    }
  cat(
    "Reaction: ", ifelse(!is.na(x[["name"]]), x[["name"]], ""), "\n",
    " - Propensity: ", x[["propensity"]], "\n",
    " - Effects: ", effect_str, "\n",
    sep = ""
  )
}

# TODO:
# Replace reaction propensity functions with rstan?
# http://tobiasmadsen.com/2017/03/31/automatic_differentiation_in_r/
