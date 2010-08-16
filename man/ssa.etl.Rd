\name{ssa.etl}
\alias{ssa.etl}

\title{Explicit tau-leap method (ETL)}
\description{Explicit tau-leap method implementation of the \acronym{SSA} as described by Gillespie (2001). It is usually called from within \code{\link{ssa}}, but can be invoked directly.}
\usage{
ssa.etl(a = stop("missing propensity vector (a)"), 
       nu = stop("missing state-change matrix (nu)"),
      tau = stop("missing step size (tau)"))
}

\arguments{
\item{a}{vector of evaluated propensity functions.} 
\item{nu}{state-change matrix.} 
\item{tau}{tau-leap.}
}

\details{Performs one time step using the Explicit tau-leap method. Intended to be invoked by \code{\link{ssa}}.}

\value{A list with two elements, 1) the time leap (\code{tau}) and 2) the realized state change vector (\code{nu_j}).}

\references{Gillespie (2001)}

\seealso{
  \link{GillespieSSA-package},
  \code{\link{ssa}}
}

\examples{
a = function(parms,x){
 b <- parms[1]
 d <- parms[2]
 K <- parms[3]
 N <- x[1]
 return(c(b*N , N*b + (b-d)*N/K))
} 
parms <- c(2,1,1000,500)
x <- 500
nu <- matrix(c(+1, -1),ncol=2)
t <- 0
for (i in seq(100)) {
  out <- ssa.etl(a(parms,x),nu,tau=0.3)
  x <- x + out$nu_j
  t <- t + 1
  cat("t:",t,", x:",x,"\n")
}
}

\keyword{misc}
\keyword{datagen}
\keyword{ts} 
