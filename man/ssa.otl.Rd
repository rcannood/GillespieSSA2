% $Id: ssa.otl.Rd 213 2008-01-14 20:39:21Z pineda $
\name{ssa.otl}
\alias{ssa.otl}

\title{Optimized tau-leap method (OTL)}
\description{Optimized tau-leap method implementation of the \acronym{SSA} as described by Cao et al. (2006). It is usually called from within \code{\link{ssa}}, but can be invoked directly.}

\usage{
ssa.otl(x = stop("missing state vector (x)"),
        a = stop("missing propensity vector (a)"),
       nu = stop("missing state-change matrix (nu)"),
      hor = stop("missing highest order reaction vector (hot)"),
       nc = stop("missing critical reactions threshold parameter (nc)"),
  epsilon = stop("missing error control parameter"),
      dtf = stop("missing direct method threshold factor (dtf)"),
       nd = stop("missing OTL suspension duration parameter (nd)"))
}

\arguments{
\item{x}{state vector.}
\item{a}{vector of evaluated propensity functions.} 
\item{nu}{state-change matrix.} 
\item{hor}{highest order reaction vector (one entry per species in \code{x})}
\item{nc}{number of critical reactions threshold parameter.}
\item{epsilon}{error control parameter.}
\item{dtf}{Direct method threshold factor for temporarily suspending the \code{OTL} method.} 
\item{nd}{number of Direct method steps to perform during an \code{OTL} suspension.}
}

\details{Performs one time step using the Explicit tau-leap method. Intended to be invoked by \code{\link{ssa}}.}

\value{A list with three elements, 1) the time leap (\code{tau}) and 2) the realized state change vector (\code{nu_j}), and 3) a boolean value (\code{suspendedTauLeapMethod}) indicating if the simulation should revert to the Direct method for \code{nd} time steps.}

\references{Cao et al. (2006)}

\note{Third order-reactions (\eqn{S_1+S_2+S_3 \rightarrow \ldots}{S_1 + S_2 + S_3 ---> ...}) are not supported currently since they are approximations to sets of coupled first- and second-order reactions). See Cao et al. (2006) for more details.}

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
  out <- ssa.otl(x,a(parms,x),nu,hor=1,nc=10,epsilon=0.03,dtf=10,nd=100)  
  x <- x + out$nu_j
  t <- t + 1
  cat("t:",t,", x:",x,"\n")
}
}

\keyword{misc}
\keyword{datagen}
\keyword{ts} 
