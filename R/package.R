#' Gillespie's Stochastic Simulation Algorithm for impatient people.
#'
#' Package description and overview of basic SSA theory
#'
#' \pkg{fastgssa} is a fast and versatile framework for various Monte Carlo
#' implementations of the stochastic simulation algorithm (\acronym{SSA}).
#' It is conceptually based on \pkg{GillespieSSA}, but rewritten entirely in
#' Rcpp in order to make \pkg{fastgssa} blazingly fast.
#' The SSA methods currently implemented are: the Direct, the Explicit
#' tau-leaping (\acronym{ETL}), and the Binomial tau-leaping (\acronym{BTL}) method.
#'
#' @name fastgssa
#' @aliases fastgssa-package fastgssa
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
#' * [ssa_plot()]: A standard visualisation for generating an overview plot fo the output.
#' * [ssa_direct()], [ssa_etl()], [ssa_btl()]: Different \acronym{SSA} algorithms.
#' * [ssa_em()]: Not actually an \acronym{SSA} algorithm.
#' * [compile_propensity_functions()]: A function for precompiling the propensity functions.
#'
#' @seealso [ssa()] for more explanation on how to use \pkg{fastgssa}
#'
#' @useDynLib fastgssa
#'
#' @importFrom dplyr first
#' @importFrom tibble lst
#' @importFrom assertthat assert_that
#' @importFrom purrr %>% map map_df map_chr map_lgl map_int map_dbl keep discard invoke pmap map2 map2_df set_names imap
NULL




# @references \itemize{
#   \item Brown D. and Rothery P. 1993. Models in biology: mathematics, statistics, and computing. John Wiley & Sons. ISBN: 0471933228.
#   \item Cao Y., Li H., and Petzold L. 2004. Efficient formulation of the stochastic simulation algorithm for chemically reacting systems. J. Chem. Phys. 121:4059-4067. \url{http://dx.doi.org/10.1063/1.1778376 }
#   \item Cao Y., Gillespie D.T., and Petzold L.R. 2006. Efficient step size selection for the tau-leaping method. J. Chem. Phys. 124:044109. \url{http://dx.doi.org/10.1063/1.2159468}
#   \item Cao Y., Gillespie D.T., and Petzold L.R. 2007. Adaptive explicit tau-leap method with automatic tau selection. J. Chem. Phys. 126:224101. \url{http://dx.doi.org/10.1063/1.2745299 }
#   \item Chatterjee A., Vlachos D.G., and Katsoulakis M.A. 2005. Binomial distribution based tau-leap accelerated stochastic simulation. J. Chem. Phys. 122:024112. \url{http://dx.doi.org/10.1063/1.1833357}
#   \item Gillespie D.T. 1977. Exact stochastic simulation of coupled chemical reactions. J. Phys. Chem. 81:2340. \url{http://dx.doi.org/10.1021/j100540a008}
#   \item Gillespie D.T. 2001. Approximate accelerated stochastic simulation of chemically reacting systems. J. Chem. Phys. 115:1716-1733. \url{http://dx.doi.org/10.1063/1.1378322 }
#   \item Gillespie D.T. 2007. Stochastic simulation of chemical kinetics. Annu. Rev. Chem. 58:35 \url{http://dx.doi.org/10.1146/annurev.physchem.58.032806.104637}
#   \item Kot M. 2001. Elements of mathematical ecology. Cambridge University Press. \url{http://dx.doi.org/10.2277/052180213X}
#   \item Pineda-Krch M. 2008. Implementing the stochastic simulation algorithm in R. Submitted to the Journal of Statistical Software 25(12): 1-18. \url{http://www.jstatsoft.org/v25/i12}
#   \item Pineda-Krch M., Blok H.J., Dieckmann U., and Doebeli M. 2007. A tale of two cycles --- distinguishing quasi-cycles and limit cycles in finite predator-prey populations. Oikos 116:53-64. \url{http://dx.doi.org/10.1111/j.2006.0030-1299.14940.x}
# }
