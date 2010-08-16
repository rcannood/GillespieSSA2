\name{ssa.btl}
\alias{ssa.btl}

\title{Binomial tau-leap method (BTL)}
\description{Binomial tau-leap method implementation of the \acronym{SSA} as described by Chatterjee et al. (2005). It is usually called from within \code{\link{ssa}}, but can be invoked directly.}
\usage{
ssa.btl(x = stop("missing state vector (x)"),
        a = stop("missing propensity vector (a)"),
       nu = stop("missing state-change matrix (nu)"),
        f = stop("missing coarse-graining factor (f)"))
}

\arguments{
\item{x}{state vector.}
\item{a}{vector of evaluated propensity functions.} 
\item{nu}{state-change matrix.} 
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
  out <- ssa.btl(x,a(parms,x),nu,f=10)  
  x <- x + out$nu_j
  t <- t + 1
  cat("t:",t,", x:",x,"\n")
}
}

\keyword{misc}
\keyword{datagen}
\keyword{ts} 
