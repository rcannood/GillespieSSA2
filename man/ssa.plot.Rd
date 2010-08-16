\name{ssa.plot}
\alias{ssa.plot}

\title{Simple plotting of ssa output}
\description{Provides basic functionally for simple and quick time series plot of simulation output from \code{\link{ssa}}.
}
\usage{
ssa.plot(        out = stop("requires simulation output object"),
                file = "ssaplot",
                  by = 1,
           plot.from = 2,
             plot.to = dim(out$data)[2],
             plot.by = 1,
          show.title = TRUE,
         show.legend = TRUE)
}
\arguments{
\item{out}{data object returned from \code{\link{ssa}}.}
\item{file}{name of the output file (only applicable if \code{plot.device!="x11"}.}
\item{by}{time increment in the plotted time series}
\item{plot.from}{first population to plot the time series for (see note)}
\item{plot.to}{last population to plot the time series for (see note)}
\item{plot.by}{increment in the sequence of populations to plot the time series for (see note)}
\item{show.title}{boolean object indicating if the plot should display a title}
\item{show.legend}{boolean object indicating if the legend is displayed}
}

\note{The options \code{by}, \code{plot.from}, \code{plot.to}, and \code{plot.by} can be used to plot a sparser sequence of data points. To plot the population sizes using a larger time interval the \code{by} option can be set, e.g. to plot only every 10th time point \code{by=10}. To plot only specific populations the \code{plot.from}, \code{plot.to}, and \code{plot.by} options can be set to subset the state vector. Note that the indexing of the populations is based on the \eqn{(t,\mathbf{X})}{(t,X)} vector, i.e. the first column is the time vector while the first population is index by 2 and the last population by \eqn{N+1}. Display of a plot title above the plot and legend is optional (and are set with the arguments {{{show.title}}} and {{{show.legend}}}. Above the plot panel miscellaneous information for the simulation are displayed, i.e. method, elapsed wall time, number of time steps executed, and the number of time steps per data point.}

\seealso{
  \link{GillespieSSA-package},
  \code{\link{ssa}}
}

\examples{
\dontrun{
## Define the Kermack-McKendrick SIR model and run once using the Direct method
parms <- c(beta=.001, gamma=.100)
x0 <- c(S=500, I=1, R=0)                      # Initial state vector
nu <- matrix(c(-1,0,1,-1,0,1),nrow=3,byrow=T) # State-change matrix
a  <- c("beta*S*I", "gamma*I")                # Propensity vector
tf <- 100                                     # Final time
simName <- "Kermack-McKendrick SIR"
out <- ssa(x0,a,nu,parms,tf,method="D",simName,verbose=TRUE,consoleInterval=1)
}

\dontrun{
## Basic ssa plot
ssa.plot(out) 
}

\dontrun{
# Plot only the infectious class
ssa.plot(out,plot.from=3,plot.to=3)
}

\dontrun{
## Multipanel plot using different SSA methods
layout(matrix(seq(4),ncol=4,byrow=TRUE))

## Using the Direct method
ssa.plot(out)

## Run and plot results using the ETL method
out <- ssa(x0,a,nu,parms,tf=100,method="ETL,simName="Kermack-McKendrick SIR")
ssa.plot(out,show.title=FALSE,show.legend=FALSE)

## Run and plot results using the BTL method
out <- ssa(x0,a,nu,parms,tf=100,method="BTL,simName="Kermack-McKendrick SIR")
ssa.plot(out,show.title=FALSE,show.legend=FALSE)

## Run and plot results using the OTL method
out <- ssa(x0,a,nu,parms,tf=100,method="OTL,simName="Kermack-McKendrick SIR")
ssa.plot(out,show.title=FALSE,show.legend=FALSE)
}
}

\keyword{misc}
\keyword{datagen}
\keyword{ts}
\keyword{device}
\keyword{utilities}
\keyword{hplot}
