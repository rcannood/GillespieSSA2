#' Explicit tau-leap method (ETL)
#'
#' Explicit tau-leap method implementation of the \acronym{SSA} as described by Gillespie (2001).
#' Note that this method does not attempt to select an appropriate value for tau, nor does it
#' implement estimated-midpoint technique.
#'
#' @param tau the step-size (default 0.3).
#'
#' @return an object of to be used by [ssa()].
#'
#' @references Gillespie D.T. 2001. Approximate accelerated stochastic simulation of chemically reacting systems. J. Chem. Phys. 115:1716-1733. \doi{10.1063/1.1378322}.
#'
#' @export
ssa_etl <- function(tau = .3) {
  ssa_method(
    name = "ETL",
    params = list(tau = tau),
    factory = make_ssa_etl
  )
}
