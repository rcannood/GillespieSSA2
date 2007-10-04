# $Id: ssa.etl.diag.R 155 2007-10-04 06:19:46Z pineda $

`ssa.etl.diag` <- function(a, nu_tile, tau) {
  MU <- length(a)         # Toto nr of reaction channels
  k  <- rpois(MU,(a*tau)) # Nr of firings per channel
  M  <- dim(nu_tile)[2]   # Nr of reaction channel per patch (nu_tile)
  U  <- MU/M              # Nr of tilings
  nu_j <- NULL
  for(f in (seq(U)-1))
    nu_j <- c(nu_j, rowSums(matrix(rep(k[1:M+f*M],dim(nu_tile)[1]),byrow=TRUE,ncol=M)*nu_tile))
  return(list(tau=tau, nu_j=nu_j))
}

