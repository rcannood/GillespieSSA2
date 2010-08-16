\name{ssa.otl.diag}
\alias{ssa.otl.diag}

\title{Optimized tau-leap method (OTL) for nu-diagonalized systems}
\description{Optimized tau-leap method for nu-diagonalized systems.}

\usage{
ssa.otl.diag(x,a,nu_tile,hor,nc,epsilon,dtf,nd)
}

\arguments{
\item{x}{state vector.}
\item{a}{vector of evaluated propensity functions.} 
\item{nu_tile}{state-change matrix.} 
\item{hor}{highest order reaction vector (one entry per species in \code{x})}
\item{nc}{number of critical reactions threshold parameter.}
\item{epsilon}{error control parameter.}
\item{dtf}{Direct method threshold factor for temporarily suspending the \code{OTL} method.} 
\item{nd}{number of Direct method steps to perform during an \code{OTL} suspension.}
}

\details{Performs one time step using the Explicit tau-leap method. It is usually called from within \code{\link{ssa}}, but can be invoked directly, see \code{\link{ssa.otl}} for Examples.}

\value{A list with three elements, 1) the time leap (\code{tau}) and 2) the realized state change vector (\code{nu_j}), and 3) a boolean value (\code{suspendedTauLeapMethod}) indicating if the simulation should revert to the Direct method for \code{nd} time steps.}

\note{Third order-reactions (\eqn{S_1+S_2+S_3 \rightarrow \ldots}{S_1 + S_2 + S_3 ---> ...}) are not supported currently since they are approximations to sets of coupled first- and second-order reactions). See Cao et al. (2006) for more details.}

\seealso{
  \code{\link{ssa.otl}},
}

\examples{
## Not intended to be invoked stand alone
}

\keyword{misc}
\keyword{datagen}
\keyword{ts} 
