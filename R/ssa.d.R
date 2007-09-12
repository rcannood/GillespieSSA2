# $Id: ssa.d.R 101 2007-08-24 03:14:56Z pineda $

`ssa.d` <-
function(a, nu) {
  j    <- sample(seq_len(length(a)), size=1, prob=a)
  nu_j <- nu[,j]
  tau  <- -log(runif(1))/sum(a)
  return(list(tau=tau, nu_j=nu_j))
}