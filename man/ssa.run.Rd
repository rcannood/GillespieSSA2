% $Id: ssa.d.Rd 76 2007-08-19 04:30:03Z pineda $
\name{ssa.run}
\alias{ssa.run}

\title{Higher-level interface to the method functions}
\description{Higher-level interface to the method functions.}
\usage{
ssa.run(x0,a,nu,parms,tf,method,tau,f,epsilon,nc,hor,dtf,nd,
        ignoreNegativeState,consoleInterval,censusInterval,
        verbose,maxWallTime)
}

\arguments{
  \item{x0}{initial states vector.}
  
  \item{a}{vector of propensity functions.}
  
  \item{nu}{state-change matrix.}
  
  \item{parms}{vector of model parameters.}
  
  \item{tf}{final time.}
  
  \item{method}{ssa method to use.}
  
  \item{tau}{step size for the \code{ETL} method (\eqn{>0}).}
  
  \item{f}{coarse-graining factor for the \code{BTL} method (\eqn{>1}) where a higher value results in larger step-size.}
  
  \item{epsilon}{accuracy control parameter for the \code{OTL} method (\eqn{>0}).}
  
  \item{nc}{critical firing threshold for the \code{OTL} method (positive integer).}
  
  \item{hor}{numerical vector of the highest order reaction for each species where \eqn{\mathtt{hor} \in \{1,2,22\}}{hor=(1,2,22)}. Only applicable in the \code{OTL} method.} 
  
  \item{dtf}{\code{D} method threshold factor for the \code{OTL} method. The \code{OTL} method is suspended if \code{tau} it estimates is smaller than the \code{dtf} multiple of the \code{tau} that the \code{D} method would have used (i.e. \eqn{\tau_{\mathtt{OTL}} < \mathtt{dtf} \times \tau_{\mathtt{D}}}{tau_OTL<dtf*\tau_D}) (See step 3, page 3 in Cao et al. 2006).}
  
  \item{nd}{number of single-reaction steps performed using the Direct method during \code{otl} suspension (See step 3, page 3, Cao et al. 2006).}
  
  \item{ignoreNegativeState}{boolean object indicating if negative state values should be ignored (this can occur in the \code{etl} method). If \code{ignoreNegativeState=TRUE} the simulation finishes gracefully when encountering a negative population size (i.e. does not throw an error). If \code{ignoreNegativeState=FALSE} the simulation stops with an error message when encountering a negative population size.}
  
  \item{consoleInterval}{(approximate) interval at which \code{ssa} produces simulation status output on the console (assumes \code{verbose=TRUE}). If \code{consoleInterval=0} console output is generated each time step (or tau-leap). If \code{consoleInterval=Inf} no console output is generated. Note, \code{verbose=FALSE} disables all console output. \strong{Console output drastically slows down simulations.}} 

  \item{censusInterval}{(approximate) interval between recording the state of the system. If \code{censusInterval=0} \eqn{(t,x)} is recorded at each time step (or tau-leap). If \code{censusInterval=Inf} only \eqn{(t_0,x_0)}{(t0,x0)} and \eqn{(t_f,x_t)}{(tf,xf)} is recorded. Note, the size of the time step (or tau-leaps) ultimatelly limits the interval between subsequent recordings of the system state since the state of the system cannot be recorded at a finer time interval the size of the time steps (or tau-leaps).}  

  \item{verbose}{boolean object indicating if the status of the simulation simulation should be displayed on the console. If \code{verbose=TRUE} the elapsed wall time and \eqn{(t,x)} is displayed on the console every \code{consoleInterval} time step and a brief summary is displayed at the end of the simulation. If \code{verbose=FALSE} the simulation runs \emph{entirely} silent (overriding \code{consoleInterval}). \strong{Verbose runs drastically slows down simulations.}} 
  
  \item{maxWallTime}{maximum wall time duration (in seconds) that the simulation is allowed to run for before terminated. This option is usefull, in particular, for systems that can end up growing uncontrollably.}
} 

\details{Invokes a specific method function until the termination criteria are fullfilled. Updates the state vector, time, and re-evaluates the propensity functions inbetween time steps. Also collects simulation data and returns it as a list object. This function is called from within \code{\link{ssa}} and is not intended to be invoked stand alone.}

\value{Returns a list object with the following elements,
  \item{timeSeries}{a numerical matrix object of the simulation time series where the first column is the time vector and subsequent columns are the state frequencies.} 
  \item{eval_a}{vector of the evaluated propensity functions.}
  \item{elapsedWallTime}{elapsed wall time in seconds.}
  \item{startWallTime}{start wall clock time (YYYY-mm-dd HH:MM:SS)}.
  \item{endWallTime}{end wall clock time (YYYY-mm-dd HH:MM:SS).}
  \item{stepSize}{vector of step sizes (i.e. time increments).}
  \item{nSuspendedTauLeaps}{number of steps performed using the Direct method due to \code{OTL} suspension (only applicable for the \code{OTL} method).}
}

\seealso{
  \code{\link{ssa}}
}

\examples{
## Not intended to be invoked stand alone
}

\keyword{misc}
\keyword{datagen}
\keyword{ts} 
