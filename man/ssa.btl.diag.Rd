% $Id: ssa.btl.Rd 76 2007-08-19 04:30:03Z pineda $
\name{ssa.btl.diag}
\alias{ssa.btl.diag}

\title{Binomial tau-leap method (BTL) for nu-diagonalized systems}
\description{Binomial tau-leap method for nu-diagonalized systems}
\usage{
ssa.btl.diag(a,nu_tile,x,f)
}

\arguments{
\item{a}{vector of evaluated propensity functions.} 
\item{nu_tile}{state-change matrix.} 
\item{x}{state vector.}
\item{f}{coarse-graining factor (see page 4 in Chatterjee et al. 2005).}
}

\details{Performs one time step using the Binomial tau-leap method. Intended to be invoked by \code{\link{ssa}}.}

\value{A list with two elements, 1) the time leap (\code{tau}) and 2) the realized state change vector (\code{nu_j}).}


\seealso{
  \code{\link{ssa.btl}},
}

\examples{
## Not intended to be invoked stand alone.
}

\keyword{misc}
\keyword{datagen}
\keyword{ts} 
