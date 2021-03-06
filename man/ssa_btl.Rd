% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ssa_btl.R
\name{ssa_btl}
\alias{ssa_btl}
\title{Binomial tau-leap method (BTL)}
\usage{
ssa_btl(mean_firings = 10)
}
\arguments{
\item{mean_firings}{A coarse-graining factor of how many firings will occur at each iteration on average.
Depending on the propensity functions, a value for \code{mean_firings} will result in warnings generated
and a loss of accuracy.}
}
\value{
an object of to be used by \code{\link[=ssa]{ssa()}}.
}
\description{
Binomial tau-leap method implementation of the \acronym{SSA} as described by Chatterjee et al. (2005).
}
\references{
Chatterjee A., Vlachos D.G., and Katsoulakis M.A. 2005. Binomial distribution based tau-leap accelerated stochastic simulation. J. Chem. Phys. 122:024112. \doi{10.1063/1.1833357}.
}
