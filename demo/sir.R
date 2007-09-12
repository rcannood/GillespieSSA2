# $Id: sir.R 73 2007-08-18 23:39:50Z pineda $

#===============================================================================
# Kermack-McKendric SIR model (Brown & Rothery, 1993)
#===============================================================================

# The Kermack-McKendrick SIR model is defined as
# dS/dt = -beta*N*S
# dI/dt = beta*N*S - gamma*I
# dR/dt = gamma*I
#
# This model has two reaction channels,
#  S + I --beta--> I
#      I --gamma-> R

# Define parameters
assign("beta",  beta  <- .001, env=.GlobalEnv)
assign("gamma", gamma <- .100, env=.GlobalEnv)

# Define system
x0  <- c(S=500, I=1, R=0)                      # Initial state vector
nu  <- matrix(c(-1,0,1,-1,0,1),nrow=3,byrow=T) # State-change matrix
a   <- c("beta*{S}*{I}", "gamma*{I}")          # Propensity vector
tf <- 100                                      # Final time
simName <- "Kermack-McKendrick SIR"

# Run the simulations

# Direct method
out <- ssa(x0,a,nu,tf,method="D",simName)
ssa.plot(out) 

# Explicit tau-leap method
out <- ssa(x0,a,nu,tf,method="ETL",simName)
ssa.plot(out) 

# Binomial tau-leap method
out <- ssa(x0,a,nu,tf,method="BTL",simName)
ssa.plot(out) 

# Optimized tau-leap method
out <- ssa(x0,a,nu,tf,method="OTL",simName)
ssa.plot(out) 
