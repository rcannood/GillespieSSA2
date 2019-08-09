#' \pkg{GillespieSSA2}: Gillespie's Stochastic Simulation Algorithm for impatient people.
#'
#' \pkg{GillespieSSA2} is a fast, scalable, and versatile framework for simulating large systems with
#' Gillespie's Stochastic Simulation Algorithm (\acronym{SSA}). This package is the spiritual successor to
#' the \pkg{GillespieSSA} package originally written by Mario Pineda-Krch.
#'
#' GillespieSSA2 has the following added benefits:
#'
#' * The whole algorithm is run in Rcpp which results in major speed improvements (>100x).
#' Even your propensity functions (reactions) are being compiled to Rcpp!
#' * Parameters and variables have been renamed to make them easier to understand.
#' * Many unit tests try to ensure that the code works as intended.
#'
#' The SSA methods currently implemented are: Exact ([ssa_exact()]), Explicit tau-leaping ([ssa_etl()]),
#' and the Binomial tau-leaping ([ssa_btl()]).
#'
#' @name GillespieSSA2
#' @aliases GillespieSSA2-package GillespieSSA2
#' @docType package
#'
#' @section The stochastic simulation algorithm:
#' The stochastic simulation algorithm (\acronym{SSA}) is a procedure for constructing
#' simulated trajectories of finite populations in continuous time.
#' If \eqn{X_i(t)} is the number of individuals in population \eqn{i}
#' (\eqn{i = 1,\ldots,N}{i = 1,...,N}) at time \eqn{t},
#' the \acronym{SSA} estimates the state vector
#' \eqn{ \mathbf{X}(t) \equiv (X_1(t),\ldots,X_N(t)) }{ X(t) = (X_1(t),...,X_N(t))},
#' given that the system initially (at time \eqn{t_0})
#' was in state \eqn{\mathbf{X}(t_0) = \mathbf{x_0}}{X(t_0) = x_0}.
#'
#' Reactions are single instantaneous events changing at least one of the populations (e.g.
#' birth, death, movement, collision, predation, infection, etc).
#' These cause the state of the system to change over time.
#'
#' The \acronym{SSA} procedure samples the time \eqn{\tau}{tau}
#' to the next reaction \eqn{R_j} (\eqn{j = 1,\ldots,M}{j = 1,...,M})
#' and updates the system state \eqn{\mathbf{X}(t)}{X(t)} accordingly.
#'
#' Each reaction \eqn{R_j} is characterized mathematically by two quantities;
#' its state-change vector \eqn{\bm{\nu_j}} and its propensity function \eqn{a_j(\mathbf{x})}.
#' The state-change vector is defined as \eqn{\bm{\nu}_j \equiv ( \nu_{1j},\ldots,\nu_{Nj} )}{nu_j = (nu_1j,...,nu_Nj)},
#' where \eqn{ \nu_{ij} }{nu_ij} is the change in the number of individuals in
#' population \eqn{i} caused by one reaction of type \eqn{j}.
#' The propensity function is defined as \eqn{a_j(\mathbf{x})}, where
#' \eqn{a_j(\mathbf{x})dt}{a_j(x)dt} is the probability that a particular
#' reaction \eqn{j} will occur in the next infinitesimal time interval
#' \eqn{\left[t,t+dt\right]}{[t,t+dt]}.
#'
#' @section Contents of this package:
#'
#' * [ssa()]: The main entry point for running an \acronym{SSA} simulation.
#' * [autoplot.ssa()]: A standard visualisation for generating an overview plot fo the output.
#' * [ssa_exact()], [ssa_etl()], [ssa_btl()]: Different \acronym{SSA} algorithms.
#' * [ode_em()]: An \acronym{ODE} algorithm.
#' * [compile_reactions()]: A function for precompiling the reactions.
#'
#' @seealso [ssa()] for more explanation on how to use \pkg{GillespieSSA2}
#'
#' @useDynLib GillespieSSA2
#'
#' @importFrom assertthat assert_that
#' @importFrom purrr %>% map map_chr map_df map_int keep discard set_names walk

NULL

# re-enable this on the next release of rlang
# @importFrom rlang %|%

# remove this on the next release of rlang
`%|%` <- function(x, y) {
  ifelse(is.na(x), y, x)
}
