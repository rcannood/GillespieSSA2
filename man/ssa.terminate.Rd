\name{ssa.terminate}
\alias{ssa.terminate}

\title{Terminates a simulation that was invoked using ssa}
\description{Terminates a simulation that was invoked using ssa.}
\usage{
ssa.terminate(args,out.rxn,tf,method,maxWallTime,verbose)
}

\arguments{
  \item{args}{list of arguments and their values passed to \code{\link{ssa}}.}
  
  \item{out.rxn}{list object as returned from \code{\link{ssa.run}}.}
  
  \item{tf}{final time.}
  
  \item{method}{ssa method to use.}

  \item{maxWallTime}{maximum simulation wall time duration (in seconds) that the simulation is allowed to run for before terminated.}
  
  \item{verbose}{boolean value indicating if some basic simulation statistics should be displayed on the console.} 
} 

\details{Terminates an invocation of the \code{link{ssa}} wrapper function. Returns the same list object as \code{\link{ssa}}. This function is called from within \code{\link{ssa}} and is not intended to be invoked stand alone.}

\seealso{
  \code{\link{ssa}}
}

\examples{
## Not intended to be invoked stand alone
}

\keyword{misc}
\keyword{datagen}
\keyword{ts} 
