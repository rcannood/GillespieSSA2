#' Explicit tau-leap method (ETL)
#'
#' Explicit tau-leap method implementation of the \acronym{SSA} as described by Gillespie (2001).
#'
#' @param tau the step-size (default 0.3).
#'
#' @return an object of to be used by [ssa()].
#'
#' @references Gillespie D.T. 2001. Approximate accelerated stochastic simulation of chemically reacting systems. J. Chem. Phys. 115:1716-1733. \url{http://dx.doi.org/10.1063/1.1378322 }
#'
#' @export
ssa_etl <- function(tau = .3) {
  ssa_method(
    name = "ETL",
    params = lst(tau)
  )
}

configure_method.SSA_ETL <- function(method, simulation) {
  simulation$use_ssa_etl(method$params$tau)
}
