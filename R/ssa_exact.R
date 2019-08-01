#' Exact method
#'
#' Exact method implementation of the \acronym{SSA} as described by Gillespie (1977).
#'
#' @return an object of to be used by [ssa()].
#'
#' @references Gillespie (1977) <doi:10.1021/j100540a008>
#'
#' @export
ssa_exact <- function() {
  ssa_method(
    name = "exact",
    params = list(),
    factory = make_ssa_exact
  )
}
