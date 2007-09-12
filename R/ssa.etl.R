# $Id: ssa.etl.R 120 2007-08-31 05:23:50Z pineda $

`ssa.etl` <-
function(a, nu, tau) {
  M <- length(a)
  k <- rpois(M,(a*tau))
  # MPK: It should not correct to call the realized state-change vector nu_j 
  # since, in leap methods, actually is not reaction specific. In Pineda-Krch 
  # (JSS ms) it is refered to as the \bm{\tilde{nu}}
  return(list(tau=tau, nu_j=rowSums(matrix(rep(k,dim(nu)[1]),byrow=TRUE,ncol=M)*nu)))
}

