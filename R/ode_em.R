#' Euler-Maruyama method (EM)
#'
#' Euler-Maruyama method implementation of the \acronym{ODE}.
#'
#' @param tau tau parameter
#' @param noise_strength noise_strength parameter
#'
#' @export
ode_em <- function(tau = 0.01, noise_strength = 2) {
  ssa_method(
    name = "EM",
    class = "ODE_EM",
    params = lst(tau, noise_strength)
  )
}

configure_method.ODE_EM <- function(method, simulation) {
  simulation$use_ode_em(method$params$tau, method$params$noise_strength)
}
