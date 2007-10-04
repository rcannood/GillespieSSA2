% $Id: ssa.btl.Rd 155 2007-10-04 06:19:46Z pineda $
\name{ssa.btl}
\alias{ssa.btl}

\title{Binomial tau-leap method (BTL)}
\description{Binomial tau-leap method implementation of the \acronym{SSA} as described by Chatterjee et al. (2005).}
\usage{
ssa.btl(a,nu,x,f)
}

\arguments{
\item{a}{vector of evaluated propensity functions.} 
\item{nu}{state-change matrix.} 
\item{x}{state vector.}
\item{f}{coarse-graining factor (see page 4 in Chatterjee et al. 2005).}
}

\details{Performs one time step using the Binomial tau-leap method. Intended to be invoked by \code{\link{ssa}}.}

\value{A list with two elements, 1) the time leap (\code{tau}) and 2) the realized state change vector (\code{nu_j}).}

\references{Chatterjee et al. (2005)}

\seealso{
  \link{GillespieSSA-package},
  \code{\link{ssa}}
}

\examples{
## Not intended to be invoked stand alone.
}

\keyword{misc}
\keyword{datagen}
\keyword{ts} 
