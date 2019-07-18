#' Direct method (D)
#'
#' Direct method implementation of the \acronym{SSA} as described by Gillespie (1977).
#'
#' @return an object of to be used by [ssa()].
#'
#' @references Gillespie D.T. 1977. Exact stochastic simulation of coupled chemical reactions. J. Phys. Chem. 81:2340. \url{http://dx.doi.org/10.1021/j100540a008}
#'
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
