# $Id: ssa.etl.R 212 2008-01-14 14:57:13Z pineda $

`ssa.etl` <-
function(a = stop("missing propensity vector (a)"), 
        nu = stop("missing state-change matrix (nu)"),
       tau = stop("missing step size (tau)")) {
  M <- length(a)
  k <- rpois(M,(a*tau))
  # MPK: Strictly speaking it is not correct to call the realized state-change
  # vector nu_j here since, in leap methods, actually is not reaction specific.
  # In Pineda-Krch (JSS ms) it is refered to as the \bm{\tilde{nu}}
  return(list(tau=tau, nu_j=rowSums(matrix(rep(k,dim(nu)[1]),byrow=TRUE,ncol=M)*nu)))
}

