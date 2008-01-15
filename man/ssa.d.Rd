% $Id: ssa.d.Rd 213 2008-01-14 20:39:21Z pineda $
\name{ssa.d}
\alias{ssa.d}

\title{Direct method (D)}
\description{Direct method implementation of the \acronym{SSA} as described by Gillespie (1977). It is usually called from within \code{\link{ssa}}, but can be invoked directly.}
\usage{
ssa.d(a = stop("missing propensity vector (a)"), 
     nu = stop("missing state-change matrix (nu)"))
}

\arguments{
\item{a}{vector of evaluated propensity functions.} 
\item{nu}{state-change matrix.} 
}

\details{Performs one time step using the Direct method.}

\value{A list with two elements, 1) the time leap (\code{tau}) and 2) the realized state change vector (\code{nu_j}).}

\references{Gillespie (1977)}

\seealso{
  \link{GillespieSSA-package},
  \code{\link{ssa}}
}

\examples{
## Logistic growth model
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
  out <- ssa.d(a(parms,x),nu)
  x <- x + out$nu_j
  t <- t + 1
  cat("t:",t,", x:",x,"\n")
}
}

\keyword{misc}
\keyword{datagen}
\keyword{ts} 
