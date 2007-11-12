#===============================================================================
# Linear Chain System (Cao et al., 2004)
#===============================================================================

# The Linear Chain System consists of M chain reactions with
# M+1 species as follows:
# S_1 --c1--> S_2 --c2-->...--cM--> S_(M+1)

library(GillespieSSA)

# Rate parameter
parms <- c(c=1)

# Number of chain reactions
M <- 50

# Initial state vector
x0 <- c(1000, rep(0,M)) 
names(x0) <- paste("x",seq(M+1),sep="") 
  
# State-change matrix
nu <- matrix(rep(0,(M*(M+1))),ncol=M)
diag(nu) <- -1
diag(nu[2:M,]) <- +1
nu[M+1,M] <- +1

# Propensity vector
a <- as.vector(paste("c*{x",seq(M),"}",sep=""))

tf <- 10 # Final time
simName <- "Linear Chain System"

# Run the simulations

# Direct method
out <- ssa(x0,a,nu,parms,tf,method="D",simName,maxWallTime=5)
ssa.plot(out)

# Explict tau-leap method
out <- ssa(x0,a,nu,parms,tf,method="ETL",simName,tau=0.1,maxWallTime=5)
ssa.plot(out)

# Binomial tau-leap method
out <- ssa(x0,a,nu,parms,tf,method="BTL",simName,f=50,maxWallTime=5)
ssa.plot(out)

# Optimal tau-leap method
out <- ssa(x0,a,nu,parms,tf,method="OTL",simName,maxWallTime=5)
ssa.plot(out)
