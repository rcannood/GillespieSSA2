% $Id: ssa.otl.Rd 76 2007-08-19 04:30:03Z pineda $
\name{ssa.otl}
\alias{ssa.otl}

\title{Optimized tau-leap method (OTL)}
\description{Optimized tau-leap method implementation of the \acronym{SSA} as described by Cao et al. (2006).}

\usage{
ssa.otl(x,a,nu,hor,nc,epsilon,dtf,nd)
}

\arguments{
\item{x}{state vector.}
\item{a}{vector of evaluated propensity functions.} 
\item{nu}{state-change matrix.} 
\item{hor}{highest order reaction vector (one entry per species in \code{x})}
\item{nc}{number of critical reactions threshold parameter.}
\item{epsilon}{error control parameter.}
\item{dtf}{Direct method threshold factor for temporarily suspending the \code{OTL} method.} 
\item{nd}{number of Direct method steps to perform during an \code{OTL} suspension.}
}

\details{Performs one time step using the Explicit tau-leap method. Intended to be invoked by \code{\link{ssa}}.}

\value{A list with three elements, 1) the time leap (\code{tau}) and 2) the realized state change vector (\code{nu_j}), and 3) a boolean value (\code{suspendedTauLeapMethod}) indicating if the simulation should revert to the Direct method for \code{nd} time steps.}

\references{Cao et al. (2006)}

\note{Third order-reactions (\eqn{S_1+S_2+S_3 \rightarrow \ldots}{S_1 + S_2 + S_3 ---> ...}) are not supported currently since they are approximations to sets of coupled first- and second-order reactions). See Cao et al. (2006) for more details.}

\seealso{
  \link{GillespieSSA-package},
  \code{\link{ssa}},
  \link{ssa.references}
}

\examples{
## Not intended to be invoked stand alone
}

\keyword{misc}
\keyword{datagen}
\keyword{ts} 
