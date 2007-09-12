# $Id: lotka.R 73 2007-08-18 23:39:50Z pineda $

#===============================================================================
# Lotka predator-prey model (Gillespie, 1977; Kot, 2001)
#===============================================================================

# This version of the Lotka predator-prey model is given by
# dY1/dt = c1*Y1 - c2*Y1*Y2
# dY2/dt = c2*Y1*Y2 - c3*Y2
# consisting of the three reaction channels,,
#      Y1 --c1--> Y1 + Y1 
# Y1 + Y2 --c2--> Y2 + Y2 
#      Y1 --c3--> 0

# Define parameters
assign("c1", c1 <- 10,  env=.GlobalEnv)
assign("c2", c2 <- .01, env=.GlobalEnv)
assign("c3", c3 <- 10,  env=.GlobalEnv)

# Define system
x0 <- c(Y1=1000, Y2=1000)                           # Initial state vector
nu <- matrix(c(+1, -1, 0, 0, 1, -1),nrow=2,byrow=T) # State-change matrix
a  <- c("c1*{Y1}", "c2*{Y1}*{Y2}","c3*{Y2}")        # Propensity vector  
tf <- 10                                            # Final time
simName <- "Lotka predator-prey model"

# Run the simulations 

# Direct method 
out <- ssa(x0,a,nu,tf,method="D",simName="Lotka predator-prey model",maxWallTime=10)
ssa.plot(out)

# Explict tau-leap method
out <- ssa(x0,a,nu,tf,method="ETL",simName,tau=0.002)
ssa.plot(out)

# Don't run: gives wrong results
# Binomial tau-leap method
#out <- ssa(x0,a,nu,tf,method="BTL",simName,f=100)
#ssa.plot(out)

# Optimized tau-leap method
out <- ssa(x0,a,nu,tf,method="OTL",simName,epsilon=0.1)
ssa.plot(out)


