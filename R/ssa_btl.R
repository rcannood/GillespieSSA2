#' Binomial tau-leap method (BTL)
#'
#' Binomial tau-leap method implementation of the \acronym{SSA} as described by Chatterjee et al. (2005).
#'
#' @param f coarse-graining factor (see page 4 in Chatterjee et al. 2005).
#'
#' @return an object of to be used by [ssa()].
#'
#' @references Chatterjee A., Vlachos D.G., and Katsoulakis M.A. 2005. Binomial distribution based tau-leap accelerated stochastic simulation. J. Chem. Phys. 122:024112. \url{http://dx.doi.org/10.1063/1.1833357}
#'
#' @export
ssa_btl <- function(f = 10) {
  assert_that(
    is.numeric(f),
    f > 1
  )
  ssa_method(
    name = "BTL",
    params = lst(f),
    factory = function() {
      make_ssa_btl(f)
    }
  )
}
