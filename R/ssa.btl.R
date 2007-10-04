# $Id: ssa.btl.R 155 2007-10-04 06:19:46Z pineda $

`ssa.btl` <-
function(a,   # Vector of evaluated propensity functions 
         nu,  # State-change matrix 
         x,   # State vector
         f) { # Coarse-graining factor (see p.4 in Chatterjee et al. (2005))

  coercing <- FALSE

  # Calculate tau
  tau <- f/sum(a)   # Calculate the current tau
  if (tau>1) tau<-1 # Tau cannot be larger than unity!

  M <- length(a)    # Number of reaction channels
  tilde_x <- x    
  nu_j <- matrix(rep(0,length(x)))
  
  # Loop over all reaction channels having propensity fun>0 
  for (j in seq_len(M)[a>0]) {    
    if (any(nu[,j]<0)) { # do this if there are limiting reactions
      mask <- nu[,j]<0
      L <- min(floor(tilde_x[mask]/abs(nu[mask,j])))
      if (a[j]*tau>L) {
        p <- 1
        coercing <- TRUE
      }
      else p <- a[j]*tau/L  
      k <- rbinom(1,L,p)
    } 
    else { # do this if there are no limiting reactions
      k <- rpois(1,(a[j]*tau))
    }

    # Update tilde_x for the current reaction j
    tmp_nu_j <- matrix(rep(k,dim(nu)[1]), byrow=TRUE, ncol=1)*nu[,j]
    tilde_x <- tilde_x + tmp_nu_j

    # Record the current state change (it is returned by ssa.btl)
    nu_j <- nu_j + tmp_nu_j    
  } # for()

  # Throw a warning message if p was coerced to unity. Coercing implies too 
  # large steps-size due to a high coarse-graining factor (f)
  if(coercing) warning("coerced p to unity - consider lowering f")

  return(list(tau=tau, nu_j=nu_j))
}

