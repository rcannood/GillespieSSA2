% $Id: ssa.d.Rd 76 2007-08-19 04:30:03Z pineda $
\name{ssa.d.diag}
\alias{ssa.d.diag}

\title{Direct method (D) for nu-diagonalized systems}
\description{Direct method for nu-diagonalized systems.}
\usage{
ssa.d.diag(a,nu)
}

\arguments{
\item{a}{vector of evaluated propensity functions.} 
\item{nu}{state-change matrix.} 
}

\details{Performs one time step using the Direct method. Intended to be invoked by \code{\link{ssa}}.}

\value{A list with two elements, 1) the time leap (\code{tau}) and 2) the realized state change vector (\code{nu_j}).}

\seealso{
  \code{\link{ssa.d}}
}

\examples{
## Not intended to be invoked stand alone
}

\keyword{misc}
\keyword{datagen}
\keyword{ts} 
