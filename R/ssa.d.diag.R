# $Id: ssa.d.diag.R 101 2007-08-24 03:14:56Z pineda $

`ssa.d.diag` <-
function(a, nu) {
  j    <- sample(seq_len(length(a)), size=1, prob=a)
  nu_j <- ssa.nutiling(a,nu,j)
  tau  <- -log(runif(1))/sum(a)
  return(list(tau=tau, nu_j=nu_j))
}
