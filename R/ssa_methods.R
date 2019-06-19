
ssa_method <- function(name, params, factory) {
  l <- lst(
    name,
    params,
    factory
  )
  class(l) <- "fastgssa::ssa_method"
  l
}

#' Euler-Maruyama method (EM)
#'
#' Euler-Maruyama method implementatio of the \acronym{SDE} as described by Euler and Maruyama (?).
#'
#' @param tau tau parameter
#' @param noise_strength noise_strength parameter
#'
#' @importFrom stats rnorm
#'
#' @export
ssa_em <- function(tau = 0.01, noise_strength = 2) {
  ssa_method(
    name = "EM",
    params = lst(tau, noise_strength),
    factory = function() {
      make_ssa_em(tau, noise_strength)
    }
  )
}


#' Direct method (D)
#'
#' Direct method implementation of the \acronym{SSA} as described by Gillespie (1977).
#'
#' @return an object of to be used by [ssa()].
#'
#' @references Gillespie (1977)
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
#' @references Gillespie (2001)
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
#' @references Chatterjee et al. (2005)
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
