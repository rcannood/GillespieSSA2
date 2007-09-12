% $Id: GillespieSSA.Rd 76 2007-08-19 04:30:03Z pineda $
\name{GillespieSSA-package}
\alias{GillespieSSA-package}
\docType{package}

\concept{Direct method}
\concept{Explicit tau-leap method}
\concept{Binomial tau-leap method}
\concept{Optimized tau-leap method}
\concept{Poisson}
\concept{stochastic simulation algorithm}
\concept{Gillespie}
\concept{ecology}
\concept{epidemiology}
\concept{evolution}
\concept{biology}

\title{Gillespie Stochastic Simulation Algorithm package}

\description{Package description and overview of basic SSA theory}

\details{\pkg{GillespieSSA} is a versatile and extensible framework for stochastic simulation in \R and provides a simple interface to a number implementations of the stochastic simulation algorithm (\acronym{SSA}). The methods currently implemented are: the Direct method (\acronym{D}), Explicit tau-leaping (\acronym{ETL}), Binomial tau-leaping (\acronym{BTL}), and Optimized tau-leaping (\acronym{OTL}). The package also provides a library of ecological, epidemiological, and evolutionary continuous-time (demo) models that can easily be customized and extended. Currently the following models are included, Decaying-Dimerization Reaction Set, Linear Chain System, single-species logistic growth model, Lotka predator-prey model, Rosenzweig-MacArthur predator-prey model, and Kermack-McKendrick \acronym{SIR} model.}
  
\section{The stochastic simulation algorithm}{
The stochastic simulation algorithm (\acronym{SSA}) is a procedure for constructing simulated trajectories of finite populations in continuous time. If \eqn{X_i(t)} is the number of individuals in population \eqn{i} (\eqn{i=1,\ldots,N}{i=1,...,N}) at time \eqn{t} the \acronym{SSA} estimates the state vector \eqn{ \mathbf{X}(t) \equiv (X_1(t),\ldots,X_N(t)) }{ X(t) = (X_1(t),...,X_N(t))}, given that the system initially (at time \eqn{t_0}) was in state \eqn{\mathbf{X}(t_0)=\mathbf{x_0}}{X(t_0)=x_0}. Reactions, single instantaneous events changing at least one of the populations (e.g. birth, death, movement, collision, predation, infection, etc), cause the state of the system to change over time. The \acronym{SSA} procedure samples the time \eqn{\tau}{tau} to the next reaction \eqn{R_j} (\eqn{j=1,\ldots,M}{j=1,...,M}) and updates the system state \eqn{\mathbf{X}(t)}{X(t)} accordingly. Each reaction \eqn{R_j} is characterized mathematically by two quantities; its state-change vector \eqn{\bm{\nu}_j \equiv ( \nu_{1j},\ldots,\nu_{Nj} )}{nu_j = (nu_1j,...,nu_Nj)}, where \eqn{ \nu_{ij} }{nu_ij} is the change in the number of individuals in population \eqn{i} caused by one reaction of type \eqn{j} and its propensity function \eqn{a_j(\mathbf{x})}{a_j(x)}, where \eqn{a_j(\mathbf{x})dt}{a_j(x)dt} is the probability that a particular reaction \eqn{j} will occur in the next infinitesimal time interval \eqn{\left[t,t+dt\right]}{[t,t+dt]}. 
}
  
\section{SSA implementations}{
There are numerous exact Monte Carlo procedures implementing the \acronym{SSA}. Perhaps the simplest is the Direct method of Gillespie (1977. The Direct method is an exact continuous-time numerical realization of the corresponding stochastic time-evolution equation. Because the Direct method simulates one reaction at a time it is often, however, computationally too slow for practical applications. 
  
Approximate implementations of the \acronym{SSA} sacrifices exactness for large improvements in computational efficiency. The most common technique used is tau-leaping where reaction-bundles are attempted in coarse-grained time increments \eqn{\tau}{tau}. Speed-ups of several orders of magnitude compared to the Direct method are common. Tau-leaping must be used with care, however, as it is not as foolproof as the Direct method.}

\section{How to cite this package}{
\preformatted{
  author = Mario Pineda-Krch,
  title  = GillespieSSA: a stochastic simulation package for R,
  year   = 2007, 
  url    = http://pineda-krch.com/GillespieSSA
}
}

\section{Acknowledgements}{Thanks to the following people who contributed with valuable feedback, suggestions, and code: Heinrich zu Dohna (UC Davis)}

\section{License}{This package is distributed under the terms of the \acronym{GNU} General Public License (\acronym{GPL}) version 3 (or newer). Copyright 2007 Mario Pineda-Krch.\cr
\cr
This file is part of the \R package \pkg{GillespieSSA}.\cr
\cr
\pkg{GillespieSSA} is free software; you can redistribute it and/or modify it under the terms of the \acronym{GNU} General Public License as published by the Free Software Foundation; either version 3 of the License, or (at your option) any later version.\cr
\cr
\pkg{GillespieSSA} is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the \acronym{GNU} General Public License for more details.\cr
\cr
You should have received a copy of the \acronym{GNU} General Public License along with this program. If not, see <\url{http://www.gnu.org/licenses/}>.
} 

\author{Mario Pineda-Krch (\url{http://pineda-krch.com})}

\seealso{
  \code{\link{ssa}},
  \code{\link{ssa.d}},
  \code{\link{ssa.etl}},
  \code{\link{ssa.btl}},
  \code{\link{ssa.otl}},
  \code{\link{ssa.plot}},
  \link{ssa.references}
}

\keyword{package}
\keyword{distribution}
