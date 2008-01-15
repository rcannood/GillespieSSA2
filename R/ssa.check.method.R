# $Id: ssa.check.method.R 214 2008-01-15 09:03:09Z pineda $

ssa.check.method <- function(x0,a,nu,method,tau,f) {

  # Check the consistency of the system dimensions, i.e. number of rows and 
  # columns in the state-change matrix and the number of elements in the initial 
  # state vector and the vector of propensity functions  
  if ((length(a)/dim(nu)[2]) != (length(x0)/dim(nu)[1])) 
    stop("inconsistent system dimensions (unequal 'nu' tessallation)")
  if (((length(a)%%dim(nu)[2])>0) || ((length(x0)%%dim(nu)[1])>0)) 
  stop("inconsistent system dimensions (fractional tessallation)")

  # For the ETL method tau>0
  if ((method=="ETL") & (!(tau>0))) stop("ETL method requires tau>0") 

  # Check that f (used in the BTL method) is >1 
  if (method=="BTL" & f<=1) stop("f has to be >1") 
}