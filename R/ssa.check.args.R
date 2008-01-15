# $Id: ssa.check.args.R 214 2008-01-15 09:03:09Z pineda $

ssa.check.args <- function(x0,a,nu,tf,method,tau,f,epsilon,nc,hor,dtf,nd,ignoreNegativeState,consoleInterval,censusInterval,verbose) {

  # Do some basic check of the argument types
  if (!is.numeric(x0))              stop("'x0' is not numeric")
  if (!is.character(a))             stop("'a' is not of character type")
  if (!is.numeric(nu))              stop("'nu' is not numeric")
  if (!is.numeric(tf))              stop("'tf' is not numeric")
  if (!is.character(method))        stop("'method' is not of character type")
  if (!is.numeric(tau))             stop("'tau' is not numeric")
  if (!is.numeric(f))               stop("'f' is not numeric")
  if (!is.numeric(epsilon))         stop("'epsilon' is not numeric")
  if (!is.numeric(nc))              stop("'nc' is not numeric")
  if (!is.numeric(hor))             stop("'hor' is not numeric")
  if (!is.numeric(dtf))             stop("'dtf' is not numeric")
  if (!is.numeric(nd))              stop("'nd' is not numeric")
  if (!is.numeric(consoleInterval)) stop("'consoleInterval' is not numeric")
  if (!is.numeric(censusInterval))  stop("'censusInterval' is not numeric")
  if ((ignoreNegativeState != TRUE) & (ignoreNegativeState != FALSE)) 
    stop("'ignoreNegativeState' is not boolean")
  if ((verbose != TRUE) & (verbose != FALSE)) stop("'verbose' is not boolean")
  if (is.null(names(x0))) stop("'x0' is missing element names")
}
