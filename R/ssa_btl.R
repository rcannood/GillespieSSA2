#' Binomial tau-leap method (BTL)
#'
#' Binomial tau-leap method implementation of the \acronym{SSA} as described by Chatterjee et al. (2005).
#'
#' @param mean_firings A coarse-graining factor of how many firings will occur at each iteration on average.
#'   Depending on the propensity functions, a value for `mean_firings` will result in warnings generated
#'   and a loss of accuracy.
#'
#' @return an object of to be used by [ssa()].
#'
#' @references Chatterjee A., Vlachos D.G., and Katsoulakis M.A. 2005. Binomial distribution based tau-leap accelerated stochastic simulation. J. Chem. Phys. 122:024112. \doi{10.1063/1.1833357}.
#'
#'
#' @export
ssa_btl <- function(mean_firings = 10) {
  assert_that(
    is.numeric(mean_firings),
    mean_firings >= 1
  )
  ssa_method(
    name = "BTL",
    params = list(mean_firings = mean_firings),
    factory = make_ssa_btl
  )
}
