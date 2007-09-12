% $Id: ssa.etl.Rd 76 2007-08-19 04:30:03Z pineda $
\name{ssa.etl}
\alias{ssa.etl}

\title{Explicit tau-leap method (ETL)}
\description{Explicit tau-leap method implementation of the \acronym{SSA} as described by Gillespie (2001).}
\usage{
ssa.etl(a,nu,tau)
}

\arguments{
\item{a}{vector of evaluated propensity functions.} 
\item{nu}{state-change matrix.} 
\item{tau}{tau-leap.}
}

\details{Performs one time step using the Explicit tau-leap method. Intended to be invoked by \code{\link{ssa}}.}}

\value{A list with two elements, 1) the time leap (\code{tau}) and 2) the realized state change vector (\code{nu_j}).}

\references{Gillespie (2001)}

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
