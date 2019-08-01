#' Explicit tau-leap method (ETL)
#'
#' Explicit tau-leap method implementation of the \acronym{SSA} as described by Gillespie (2001).
#'
#' @param tau the step-size (default 0.3).
#'
#' @return an object of to be used by [ssa()].
#'
#' @references Gillespie (2001) <doi:10.1063/1.1378322>
#'
#' @export
ssa_etl <- function(tau = .3) {
  ssa_method(
    name = "ETL",
    params = list(tau = tau),
    factory = function() {
      make_ssa_etl(tau)
    }
  )
}
