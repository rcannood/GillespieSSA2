# $Id: ssa.d.R 212 2008-01-14 14:57:13Z pineda $

`ssa.d` <-
function(a = stop("missing propensity vector (a)"), 
        nu = stop("missing state-change matrix (nu)")) {
  j    <- sample(seq(length(a)), size=1, prob=a)
  nu_j <- nu[,j]
  tau  <- -log(runif(1))/sum(a)
  return(list(tau=tau, nu_j=nu_j))
}