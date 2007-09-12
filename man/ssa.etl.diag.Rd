% $Id: ssa.etl.Rd 76 2007-08-19 04:30:03Z pineda $
\name{ssa.etl.diag}
\alias{ssa.etl.diag}

\title{Explicit tau-leap method (ETL) for nu-diagonalized systems}
\description{Explicit tau-leap method for nu-diagonalized systems.}
\usage{
ssa.etl.diag(a,nu_tile,tau)
}

\arguments{
\item{a}{vector of evaluated propensity functions.} 
\item{nu_tile}{state-change matrix.} 
\item{tau}{tau-leap.}
}

\details{Performs one time step using the Explicit tau-leap method. Intended to be invoked by \code{\link{ssa}}.}}

\value{A list with two elements, 1) the time leap (\code{tau}) and 2) the realized state change vector (\code{nu_j}).}

\seealso{
  \code{\link{ssa.etl}},
}

\examples{
## Not intended to be invoked stand alone
}

\keyword{misc}
\keyword{datagen}
\keyword{ts} 
