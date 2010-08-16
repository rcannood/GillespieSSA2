\name{ssa.check.method}
\alias{ssa.check.method}

\title{Validates consistency of the system definition}
\description{Validates consistency of the system definition.}
\usage{
ssa.check.method(x0,a,nu,method,tau,f)
}

\arguments{
  \item{x0}{numerical vector of initial states where the component elements must be named using the same notation as the corresponding state variable in the propensity vector, \code{a}.}
  
  \item{a}{character vector of propensity functions where state variables correspond to the names of the elements in \code{x0}.}
  
  \item{nu}{numerical matrix of change if the number of individuals in each state (rows) caused by a single reaction of any given type (columns).}
  
  \item{method}{text string indicating the \acronym{SSA} method to use, the valid options are: \code{D} --- Direct method (default method), \code{ETL} - Explicit tau-leap, \code{BTL} --- Binomial tau-leap, or \code{OTL} --- Optimized tau-leap.}
  
  \item{tau}{step size for the \code{ETL} method (\eqn{>0}).}
  
  \item{f}{coarse-graining factor for the \code{BTL} method (\eqn{>1}) where a higher value results in larger step-size.}
} 

\details{Performs a few basic consistency checks the defined system, e.g. that the number of rows and columns in the state-change matrix and the number of elements in the initial state vector and the vector of propensity functions are consistent. This function is called from within \code{\link{ssa}} and is not intended to be invoked stand alone.}

\seealso{
  \code{\link{ssa}}
  \code{\link{ssa.check.args}}
}

\examples{
## Not intended to be invoked stand alone
}

\keyword{misc}
\keyword{datagen}
\keyword{ts} 
