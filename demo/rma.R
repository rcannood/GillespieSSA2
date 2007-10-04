#===============================================================================
# Rosenzweig-MacArthur predator-prey model (Pineda-Krch et al., 2007)
#===============================================================================

# The Rosenzweig-MacArthur predator-prey is defined as
# dN/dt = r(1-N/K - alpha/(1+wN))NP
# dP/dt = c*alpha/(1+wN))NP
# 
# The model has five reaction channels,
#     N --BN--> N + N    Prey birh
#     N --DN--> 0        Prey death 
# N + P --NP--> P + P    Predation
#     P --BP--> P + P    Predator birth
#     P --DP--> 0        Predator death
#
# where
# BN = b
# DN = d+(b-d)N/K
# NP = alpha/(1+wN)
# BP = c*alpha/(1+wN)N
# DP = g
#
# Propensity functions:
# a1 = b * N
# a2 = (d+(b-d)N/K) * N
# a3 = alpha/(1+wN) * N * P
# a4 = c*alpha/(1+wN) * N * P
# a5 = g * P

# Define parameters 
# (B in Figure 1 in Pineda-Krch et al. 2007)
parms <- c(b=2, d=1, K=1000, alpha=0.005, 
           w=0.0025, c=2, g=2)       # Parameters
x0  <- c(N=500, P=500)               # Initial state vector
nu  <- matrix(c(+1, -1, -1,  0,  0,  # State-change matrix
                 0,  0,  0, +1, -1),     
                 nrow=2,byrow=TRUE) 
a   <- c("b*{N}",                    # Propensity vector
         "(d+(b-d)*{N}/K)*{N}",
         "alpha/(1+w*{N})*{N}*{P}",
         "c*alpha/(1+w*{N})*{N}*{P}",
         "g*{P}")   

tf <- 100
simName <- "Rosenzweig-MacArthur predator-prey model"

# Run the simulations

# Direct method
out <- ssa(x0,a,nu,parms,tf,method="D",simName,maxWallTime=10)
ssa.plot(out)

# Explicit tau-leap method
out <- ssa(x0,a,nu,parms,tf,method="ETL",simName,tau=0.01)
ssa.plot(out)

# Don't run: wrong results
# Binomial tau-leap method
# out <- ssa(x0,a,nu,parms,tf,method="BTL",simName)
# ssa.plot(out)

# Optimized tau-leap method
out <- ssa(x0,a,nu,parms,tf,method="OTL",simName)
ssa.plot(out)
