\name{ssa.nutiling}
\alias{ssa.nutiling}

\title{Direct method nu-diagonalization mapping}
\description{Auxiliary function for \code{ssa.d.diag} performing virtual mapping of nu-diagonalized systems.}
\usage{
ssa.nutiling(a,nu,j)
}

\arguments{
\item{a}{vector of evaluated propensity functions.} 
\item{nu}{state-change matrix.}
\item{j}{Reaction index to map} 
}

\value{The virtual realized state change vector (\code{nu_j}).}

\seealso{
  \code{\link{ssa.d.diag}}
  \code{\link{ssa.d}}
}

\examples{
## Not intended to be invoked stand alone
}

\keyword{misc}
\keyword{datagen}
\keyword{ts} 
