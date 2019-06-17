
ssa_method <- function(name, params, factory) {
  l <- lst(
    name,
    params,
    factory
  )
  class(l) <- "fastgssa::ssa_method"
  l
}

#' Euler-Maruyama method (EM)
#'
#' Euler-Maruyama method implementation
#'
#' @param h h parameter
#' @param noise_strength noise_strength parameter
#'
#' @importFrom stats rnorm
#'
#' @export
ssa_em <- function(h = 0.01, noise_strength = 2) {
  ssa_method(
    name = "em",
    params = lst(h, noise_strength),
    factory = function() {
      make_ssa_em(h, noise_strength)
    }
  )
}


#' Direct method (D)
#'
#' Direct method implementation of the \acronym{SSA} as described by Gillespie (1977).
#'
#' @return an object of to be used by \code{\link{ssa}}.
#' @export
ssa_direct <- function() {
  ssa_method(
    name = "direct",
    params = list(),
    factory = function() {
      make_ssa_direct()
    }
  )
}
