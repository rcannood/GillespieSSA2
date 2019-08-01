#' Euler-Maruyama method (EM)
#'
#' Euler-Maruyama method implementation of the \acronym{ODE}.
#'
#' @param tau tau parameter
#' @param noise_strength noise_strength parameter
#'
#' @return an object of to be used by [ssa()].
#'
#' @export
ode_em <- function(tau = 0.01, noise_strength = 2) {
  ssa_method(
    name = "EM",
    class = "ODE_EM",
    params = list(tau = tau, noise_strength = noise_strength),
    factory = make_ode_em
  )
}
