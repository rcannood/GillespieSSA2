# $Id: ssa.nutiling.R 103 2007-08-25 01:24:52Z pineda $

ssa.nutiling <- function(a,nu,j) {
  M  <- dim(nu)[2]        # Number of reaction channels in nu tile
  N  <- dim(nu)[1]        # Number of states in nu tile
  U  <- length(a)/M       # Number of tessallations of nu tile
  f  <- ceiling((j/M)-1)  # Frameshift factor
  jp <- j-f*M             # Relative reaction channel index
  nu_jp <- nu[,jp] 
  nu_j <- c(rep(0,f*N),   # Leading zeros 
            nu_jp,        # Relative state-change matrix
            rep(0,(U*N-(f*N+N)))) # Lagging zeros
  return(nu_j)
}
