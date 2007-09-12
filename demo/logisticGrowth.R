# $Id: logisticGrowth.R 73 2007-08-18 23:39:50Z pineda $

#===============================================================================
# Logistic growth (Pearl-Verhulst model) (Kot, 2001)
#===============================================================================

# The logistic growth model is given by dN/dt=rN(1-N/K) 
# where N is the number (density) of indviduals at time t, K 
# is the carrying capacity of the population, r is the 
# intrinsic growth rate of the population. We assume r=b-d 
# where b is the per capita p.c. birth rate and d is the 
# p.c. death rate. 
#
# This model consists of two reaction channels,
# N ---b--->  N + N
# N ---d'---> 0
# where d'=d+(b-d)N/K. The propensity functions are a_1=bN 
# and a_2=d'N.

# Parameters
assign("b", b <- 2, env=.GlobalEnv)
assign("d", d <- 1, env=.GlobalEnv)
assign("K", K <- 1000, env=.GlobalEnv)

x0 <- c(N=500)                          # Initial state vector
nu <- matrix(c(+1,-1),ncol=2)           # State-change matrix
a  <- c("b*{N}", "(d+(b-d)*{N}/K)*{N}") # Propensity vector
tf <- 10                                # Final time
simName <- "Logistic growth" 

# Run the simulations

# Direct method
out <- ssa(x0,a,nu,tf,method="D",simName,maxWallTime=5)
ssa.plot(out) 
 
# Explict tau-leap method
out <- ssa(x0,a,nu,tf,method="ETL",simName,tau=0.03,maxWallTime=5)
ssa.plot(out) 

# (Don't run: wrong results)
# Run the simulation using the Binomial tau-leap method 
# out <- ssa(x0,a,nu,tf,method="BTL",simName,f=5,maxWallTime=5)
# ssa.plot(out) 

# Optimized tau-leap method
out <- ssa(x0,a,nu,tf,method="OTL",simName,maxWallTime=5)
ssa.plot(out)
