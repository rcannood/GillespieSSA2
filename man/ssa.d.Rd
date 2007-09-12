% $Id: ssa.d.Rd 76 2007-08-19 04:30:03Z pineda $
\name{ssa.d}
\alias{ssa.d}

\title{Direct method (D)}
\description{Direct method implementation of the \acronym{SSA} as described by Gillespie (1977).}
\usage{
ssa.d(a,nu)
}

\arguments{
\item{a}{vector of evaluated propensity functions.} 
\item{nu}{state-change matrix.} 
}

\details{Performs one time step using the Direct method. Intended to be invoked by \code{\link{ssa}}.}

\value{A list with two elements, 1) the time leap (\code{tau}) and 2) the realized state change vector (\code{nu_j}).}

\references{Gillespie (1977)}

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
