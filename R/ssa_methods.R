
ssa_method <- function(name, params, factory) {
  l <- lst(
    name,
    params,
    factory
  )
  class(l) <- "gillespie::ssa_method"
  l
}

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
    params = lst(tau, noise_strength),
    factory = function() {
      make_ode_em(tau, noise_strength)
    }
  )
}


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
    params = lst(tau),
    factory = function() {
      make_ssa_etl(tau)
    }
  )
}

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
  ssa_method(
    name = "BTL",
    params = lst(f),
    factory = function() {
      make_ssa_btl(f)
    }
  )
}
