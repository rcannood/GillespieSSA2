% $Id: ssa.d.Rd 76 2007-08-19 04:30:03Z pineda $
\name{ssa.check.args}
\alias{ssa.check.args}

\title{Validates the arguments for the ssa wrapper function}
\description{Validates the arguments for the ssa wrapper function.}
\usage{
ssa.check.args(x0,a,nu,tf,method,tau,f,epsilon,nc,hor,dtf,nd,
               ignoreNegativeState,consoleInterval,
               censusInterval,verbose)
}

\arguments{
  \item{x0}{numerical vector of initial states where the component elements must be named using the same notation as the corresponding state variable in the propensity vector, \code{a}.}
  
  \item{a}{character vector of propensity functions where state variables correspond to the names of the elements in \code{x0}.}
  
  \item{nu}{numerical matrix of change if the number of individuals in each state (rows) caused by a single reaction of any given type (columns).}
  
  \item{tf}{final time.}
  
  \item{method}{text string indicating the \acronym{SSA} method to use, the valid options are: \code{D} --- Direct method (default method), \code{ETL} - Explicit tau-leap, \code{BTL} --- Binomial tau-leap, or \code{OTL} --- Optimized tau-leap.}
  
  \item{tau}{step size for the \code{ETL} method (\eqn{>0}).}
  
  \item{f}{coarse-graining factor for the \code{BTL} method (\eqn{>1}) where a higher value results in larger step-size.}
  
  \item{epsilon}{accuracy control parameter for the \code{OTL} method (\eqn{>0}).}
  
  \item{nc}{critical firing threshold for the \code{OTL} method (positive integer).}
  
  \item{hor}{numerical vector of the highest order reaction for each species where \eqn{\mathtt{hor} \in \{1,2,22\}}{hor=(1,2,22)}. Setting \code{hor=NaN} uses the default \code{hor=rep(22,N)} where \code{N} is the number of species (See page 6 in Cao et al. 2006). Unless \code{hor=NaN} the number of elements must equal the number of states \eqn{N}. Only applicable in the \code{OTL} method.} 
  
  \item{dtf}{\code{D} method threshold factor for the \code{OTL} method. The \code{OTL} method is suspended if \code{tau} it estimates is smaller than the \code{dtf} multiple of the \code{tau} that the \code{D} method would have used (i.e. \eqn{\tau_{\mathtt{OTL}} < \mathtt{dtf} \times \tau_{\mathtt{D}}}{tau_OTL<dtf*\tau_D}) (See step 3, page 3 in Cao et al. 2006).}
  
  \item{nd}{number of single-reaction steps performed using the Direct method during \code{otl} suspension (See step 3, page 3, Cao et al. 2006).}
  
  \item{ignoreNegativeState}{boolean object indicating if negative state values should be ignored (this can occur in the \code{etl} method). If \code{ignoreNegativeState=TRUE} the simulation finishes gracefully when encountering a negative population size (i.e. does not throw an error). If \code{ignoreNegativeState=FALSE} the simulation stops with an error message when encountering a negative population size.}
  
  \item{consoleInterval}{(approximate) interval at which \code{ssa} produces simulation status output on the console (assumes \code{verbose=TRUE}). If \code{consoleInterval=0} console output is generated each time step (or tau-leap). If \code{consoleInterval=Inf} no console output is generated. Note, \code{verbose=FALSE} disables all console output. \strong{Console output drastically slows down simulations.}} 

  \item{censusInterval}{(approximate) interval between recording the state of the system. If \code{censusInterval=0} \eqn{(t,x)} is recorded at each time step (or tau-leap). If \code{censusInterval=Inf} only \eqn{(t_0,x_0)}{(t0,x0)} and \eqn{(t_f,x_t)}{(tf,xf)} is recorded. Note, the size of the time step (or tau-leaps) ultimatelly limits the interval between subsequent recordings of the system state since the state of the system cannot be recorded at a finer time interval the size of the time steps (or tau-leaps).}  

  \item{verbose}{boolean object indicating if the status of the simulation simulation should be displayed on the console. If \code{verbose=TRUE} the elapsed wall time and \eqn{(t,x)} is displayed on the console every \code{consoleInterval} time step and a brief summary is displayed at the end of the simulation. If \code{verbose=FALSE} the simulation runs \emph{entirely} silent (overriding \code{consoleInterval}). \strong{Verbose runs drastically slows down simulations.}} 
} 

\details{Performs basic type checking of many of the arguments passed to the \code{\link{ssa}} wrapper function. Note that no logical checking is currently performed, e.g. which arguments are required with which method (see \code{\link{ssa.check.method}}). This function is called from within \code{\link{ssa}} and is not intended to be invoked stand alone.}

\seealso{
  \code{\link{ssa}}
  \code{\link{ssa.check.method}}
}

\examples{
## Not intended to be invoked stand alone
}

\keyword{misc}
\keyword{datagen}
\keyword{ts} 
