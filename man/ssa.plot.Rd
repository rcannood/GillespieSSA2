% $Id: ssa.plot.Rd 146 2007-09-12 22:04:21Z pineda $
\name{ssa.plot}
\alias{ssa.plot}

\title{Simple plotting of ssa output}
\description{Provides basic functionally for simple and quick time series plot of simulation output from \code{\link{ssa}}.
}
\usage{
ssa.plot(        out = stop("requires simulation output object"),
         plot.device = "x11",
                file = "ssaplot",
                  by = 1,
           plot.from = 2,
             plot.to = dim(out$data)[2],
             plot.by = 1)
}
\arguments{
\item{out}{data object returned from \code{\link{ssa}}.}
\item{plot.device}{character string indicating the device where plot should rendered, valid options are \code{"pdf","png", "jpeg", "bmp"}}
\item{file}{name of the output file (only applicable if \code{plot.device!="x11"}.}
\item{by}{time increment in the plotted time series}
\item{plot.from}{first population to plot the time series for (see note)}
\item{plot.to}{last population to plot the time series for (see note)}
\item{plot.by}{increment in the sequence of populations to plot the time series for (see note)}
}

\value{Rendered time series plot on \code{plot.device}.}

\note{When rendering time series with a large number of data points \code{plot.device="pdf"}, can be slow and can generate very large files compared to bitmaps (e.g. \code{plot.device="png"}, \code{plot.device="jpeg"}, and \code{plot.device="bmp"}). The options \code{by}, \code{plot.from}, \code{plot.to}, and \code{plot.by} can be used to plot a sparser sequense of data points. To plot the population sizes using a larger time interval the \code{by} option can be set, e.g. to plot only every 10th time point \code{by=10}. To plot only specific populations the \code{plot.from}, \code{plot.to}, and \code{plot.by} options can be set to subset the state vector. Note that the indexing of the populations is based on the \eqn{(t,\mathbf{X})}{(t,X)} vector, i.e. the first column is the time vector while the first population is index by 2 and the last population by \eqn{N+1}.}

\seealso{
  \link{GillespieSSA-package},
  \code{\link{ssa}}
}

\examples{
\dontrun{
## Plot to x11 device
out <- sir()
ssa.plot(out)
}

\dontrun{
## Plot to pdf device 
out <- lotka("ETL",tau=0.01)
ssa.plot(out,plot.device="jpeg")
}

\dontrun{
## Kermack-McKendrick SIR model 
x0  <- c(S=499,I=1,R=0)
a   <- c("0.001*{S}*{I}","0.1*{I}") 
nu  <- matrix(c(-1,0,+1,-1,0,+1),nrow=3,byrow=T)
out <- ssa(x0,a,nu,tf=100,simName="Kermack-McKendrick SIR")
ssa.plot(out)

# Plot only the infectious class
ssa.plot(out,plot.from=3,plot.to=3)
}
}

\keyword{misc}
\keyword{datagen}
\keyword{ts}
\keyword{device}
\keyword{utilities}
\keyword{hplot}
