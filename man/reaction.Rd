% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/reaction.R
\name{reaction}
\alias{reaction}
\title{Define a reaction}
\usage{
reaction(propensity, effect, name = NA_character_)
}
\arguments{
\item{propensity}{\verb{[character/formula]} A character or formula representation of the propensity function, written in C++.}

\item{effect}{\verb{[named integer vector]} The change in state caused by this reaction.}

\item{name}{\verb{[character]} A name for this reaction (Optional). May only contain characters matching \verb{[A-Za-z0-9_]}.}
}
\value{
\verb{[SSA_reaction]} This object describes a single reaction as part of an SSA simulation. It contains the following member values:
\itemize{
\item \code{r[["propensity"]]}: The propensity function as a character.
\item \code{r[["effect"]]}: The change in state caused by this reaction.
\item \code{r[["name"]]}: The name of the reaction, \code{NA_character_} if no name was provided.
}
}
\description{
During an SSA simulation, at any infinitesimal time interval,
a reaction will occur with a probability defined according to its
propensity. If it does, then it will change the state vector according
to its effects.
}
\details{
It is possible to use 'buffer' values in order to speed up the computation
of the propensity functions. For instance, instead of \code{"(c3 * s1) / (1 + c3 * c1)"},
it is possible to write \code{"buf = c3 * s1; buf / (buf + 1)"} instead.
}
\examples{
#        propensity                        effect
reaction(~ c1 * s1,                          c(s1 = -1))
reaction("c2 * s1 * s1",                     c(s1 = -2, s2 = +1))
reaction("buf = c3 * s1; buf / (buf + 1)",   c(s1 = +2))
}
