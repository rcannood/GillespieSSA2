% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ssa_plot.R
\name{plot_ssa}
\alias{plot_ssa}
\title{Simple plotting of ssa output}
\usage{
plot_ssa(
  ssa_out,
  state = TRUE,
  propensity = FALSE,
  buffer = FALSE,
  firings = FALSE,
  geom = c("point", "step")
)
}
\arguments{
\item{ssa_out}{Data object returned by \code{\link[=ssa]{ssa()}}.}

\item{state}{Whether or not to plot the state values.}

\item{propensity}{Whether or not to plot the propensity values.}

\item{buffer}{Whether or not to plot the buffer values.}

\item{firings}{Whether or not to plot the reaction firings values.}

\item{geom}{Which geom to use, must be one of \code{"point"}, \code{"step"}.}
}
\description{
Provides basic functionally for simple and quick time series plot of simulation output from \code{\link[=ssa]{ssa()}}.
}
