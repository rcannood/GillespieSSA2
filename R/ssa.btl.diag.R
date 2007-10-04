# $Id: ssa.btl.diag.R 155 2007-10-04 06:19:46Z pineda $

`ssa.btl.diag` <-
function(a,        # Vector of evaluated propensity functions 
         nu_tile,  # State-change matrix 
         x,        # State vector
         f) {      # Coarse-graining factor (see p.4 in Chatterjee et al. (2005))

  coercing <- FALSE

  # Calculate tau
  tau <- f/sum(a)   # Calculate the current tau
  if (tau>1) tau<-1 # Tau cannot be larger than unity!

  M <- dim(nu_tile)[2] # Number of reaction channels per nu-tile
  N <- dim(nu_tile)[1] # Number of states per nu-tile
  MU <- length(a)      # Toto nr of reaction channels
  U <- MU/M            # Nr of tilings
  tilde_x <- x    
  nu_j <- rep(0,(N*U))
  
  # Identify potential limiting reactions
  #mask <- apply(nu_tile,2,function(x) any(x<0))
  
  # Loop over all reaction channels having a non-zero (>0) propensity fun 
  for (j in seq_len(U*M)[a>0]) {
    f <- ceiling((j/M)-1)
    jp <- j-f*M  # Intra-patch reaction channel index (j->jp)
    x1 <- 1+f*N
    x2 <- 1+f*N+(N-1)
    if (any(nu_tile[,jp]<0)) {  # Do this if there are limiting reactions
      mask <- nu_tile[,jp]<0    # Which species has the limiting reaction
      tilde_xt <- tilde_x[x1:x2]
      L <- min(floor(tilde_xt[mask]/abs(nu_tile[mask,jp])))
      if (a[j]*tau>L) {
        p <- 1
        coercing <- TRUE
      }
      else {
        p <- a[j]*tau/L  
        k <- rbinom(1,L,p)
      } 
    } else { # do this if there are no limiting reactions
      k <- rpois(1,(a[j]*tau))
    }

    # Update tilde_x for the current reaction j
    tmp_nu_j <- rep(k,dim(nu_tile)[1])*nu_tile[,jp]
    tilde_x[x1:x2] <- tilde_x[x1:x2] + tmp_nu_j
    
    # Record the current state change
    nu_j[x1:x2] <- nu_j[x1:x2] + tmp_nu_j    
    if (any(is.na(nu_j))) browser() # MPK: Just in case!
  } # for()

  # Throw a warning message if p was coerced to unity. Coercing implies too 
  # large steps-size due to a high coarse-graning factor (f)
  if(coercing) warning("coerced p to unity - consider lowering f")
  return(list(tau=tau, nu_j=nu_j))
}

