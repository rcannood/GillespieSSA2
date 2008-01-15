# $Id: ssa.d.diag.R 211 2008-01-11 23:31:17Z pineda $

`ssa.d.diag` <-
function(a, nu) {
  j    <- sample(seq(length(a)), size=1, prob=a)
  nu_j <- ssa.nutiling(a,nu,j)
  tau  <- -log(runif(1))/sum(a)
  return(list(tau=tau, nu_j=nu_j))
}
