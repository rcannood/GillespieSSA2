# $Id: ssa.R 155 2007-10-04 06:19:46Z pineda $

`ssa` <-
function(         x0 = stop("undefined 'x0'"), 
                   a = stop("undefined 'a'"),
                  nu = stop("undefined 'nu'"),
               parms = NULL, 
                  tf = stop("undefined 'tf'"), 
              method = "D",
             simName = "",
                 tau = 0.3,
                   f = 10,
             epsilon = 0.03,
                  nc = 10,
                 hor = NaN,
                 dtf = 10,
                  nd = 100,
 ignoreNegativeState = TRUE,
     consoleInterval = 0,
      censusInterval = 0,
             verbose = FALSE,
         maxWallTime = Inf
             ) 
{ # End of function arguments

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

  # Start the timer and make an announcement if running silent
  procTimeStart   <- proc.time()
  elapsedWallTime <- 0
  startWallTime   <- format(Sys.time())
  if (verbose) {
    cat("Running ",method,
        " method with console output every ",consoleInterval,
        " time step\n",sep="")
    cat("Start wall time: ",startWallTime,"...\n",sep="")  
    flush.console()
  }
  
  # Convert lower case method names to upper case (undocumented featurette)
  if (method=="d")   method <- "D"
  if (method=="etl") method <- "ETL"
  if (method=="btl") method <- "BTL"
  if (method=="otl") method <- "OTL"
  
  #############################################################################
  # More elaborated sanity checks...
  #############################################################################

  # Check the consistency of the system dimensions, i.e. number of rows and 
  # columns in the state-change matrix and the number of elements in the initial 
  # state vector and the vector of propensity functions  
  if ((length(a)/dim(nu)[2]) != (length(x0)/dim(nu)[1])) 
    stop("inconsistent system dimensions (unequal 'nu' tessallation)")
  if (((length(a)%%dim(nu)[2])>0) || ((length(x0)%%dim(nu)[1])>0)) 
    stop("inconsistent system dimensions (fractional tessallation)")
  
  # Is the system nu-tiled along the diagonal?
  if ((length(a)/dim(nu)[2]>1) && (length(x0)/dim(nu)[1])>1){
    if (method=="D")   method <- "D.diag" 
    if (method=="ETL") method <- "ETL.diag"
    if (method=="BTL") method <- "BTL.diag"
    if (method=="OTL") method <- "OTL.diag"
  }

  # For the ETL method tau>0
  if ((method=="ETL") & (!(tau>0))) stop("ETL method requires tau>0") 

  # Parse the propensity vector by recursive substitution of variables names 
  # with the corresponding state vector indexing
  x <- x0
  varNam <- names(x)
  if (is.null(varNam)) stop("element labels are missing in 'x'") 
  for (i in seq_len(length(varNam))) {
   pattern <- paste("{",varNam[i],"}",sep="")
   a <- gsub(pattern,paste("x[",i,"]",sep=""),a,fixed=TRUE)
  }
  parse_a <- parse(text=a)

  # Assign the parameters defined in the parms vector
  if (!is.null(parms)) {
    parmsNames <- names(parms)
    for (i in seq(length(parms))) {
      assign(parmsNames[i],parms[[i]])
    }
  }

  # Check that f (used in the BTL method) is >1 
  if (method=="BTL" & f<=1) stop("f has to be >1") 
  
  # Take a snapshot of all the options so they can be saved later
  args <- list(    x0 = x0, 
                    a = a, 
                   nu = nu, 
                parms = parms, 
                   tf = tf, 
               method = method, 
                  tau = tau, 
                    f = f, 
              epsilon = epsilon, 
                   nc = nc, 
                  hor = hor, 
                  dtf = dtf, 
                   nd = nd, 
  ignoreNegativeState = ignoreNegativeState, 
      consoleInterval = consoleInterval, 
       censusInterval = censusInterval,
              verbose = verbose,
              simName = simName)
  
  # Initialize miscellaneous counters 
  t <- 0 # Initialize the simulation time
  timeOfNextCensus <- t + censusInterval  # Time of first data census
  timeForConsole   <- t + consoleInterval # Time of first console output
  stepSize         <- NULL
  timeToTerminate  <- FALSE

  # Add the initial state of the system to the time series matrix and 'pre-grow' with NAs
  timeSeries <- c(t, x) # First data point
  numCols    <- length(timeSeries)
  timeSeries <- rbind(timeSeries, matrix(nrow=1000, ncol=(numCols)))
  currentRow <- 2 # First row contains (t0,x0)
  
  # Set up empty vectors for the evaluated propensity functions
  M <- length(a)
  eval_a <- rep(0,M)

  # Evaluate the initial transition rates by evaluating the propensity functions
  for (num in seq_len(length(parse_a))) eval_a[num] <- eval(parse_a[num])

  # If required (depends on the solver method) check if hor vector is defined 
  # as NA in which case the user did not submitt his/her own hor vector as an 
  # argument to ssa(). Fall back to the conservative default vaule of 2 for
  # each species. If the hor vector was defined by the user check it's length 
  # (should have the same number of elements as the state vector) and check that 
  # is only consists of 1, 2, or 22 (i.e. first-, second-order, or the homo-dimer 
  # reactions).
  if (method == "OTL") {
    if (any(is.na(hor))) hor <- rep(2,length(x0)) # Undefined hor - use default values
    else if (length(hor) != length(x0)) stop("length of hor vector is different from length of 'x0'")
    else if (any(hor!=1 & hor!=2 & hor!=22)) stop("wrong value(s) in hor vector (can only be 1, 2, or 22)")
  }

  #############################################################################
  # We are ready to roll, start the simulation loop...
  # Continue the simulation loop as long as we have not reached the end time, 
  # all of the populations have not gone extincs, as long as no population is 
  # negative (occasinal by-product of the ETL method), and as long as at 
  # least one propensity function is larger than zero
  #############################################################################

  # Display the first time step on the console (not necessary if 
  # consoleInterval=0 since it is displayed every time step anyway and hence is 
  # already taken care of below
  if ((verbose) & (consoleInterval>0)) {
    cat("t=",t," : ",sep="")
    cat(x,sep=",")
    cat("\n")
    flush.console()
  }

  suspendedTauLeapMethod <- FALSE
  nSuspendedTauLeaps <- 0
  while( (t<tf) & (any(x>0)) & 
         (all(x>=0)) & (any(eval_a>0)) & 
         (elapsedWallTime<=maxWallTime) ) {

    doCalc <- TRUE
    
    if ((verbose) & (timeForConsole<=t)) {
      cat("(",elapsedWallTime,"s) t=",t," : ",sep="")
      cat(x,sep=",")
      cat("\n")
      flush.console()
      timeForConsole <- timeForConsole + consoleInterval
    }
    
    switch( method,
            "D" = { out <- ssa.d(eval_a, nu) 
                    if (suspendedTauLeapMethod) { 
                      suspendedTauLeapMethod <- suspendedTauLeapMethod-1
                      nSuspendedTauLeaps <- nSuspendedTauLeaps + 1
                      if (!suspendedTauLeapMethod) method <- "OTL"
                    }
                  },
          "ETL" = { out <- ssa.etl(eval_a, nu, tau) }, 
          "BTL" = { out <- ssa.btl(eval_a, nu, x, f) },
          "OTL" = { out <- ssa.otl(x, eval_a, nu, hor, nc, epsilon, dtf, nd)
                    suspendedTauLeapMethod <- out$suspendedTauLeapMethod
                    if (suspendedTauLeapMethod) { 
                      method <- "D"
                      doCalc <- FALSE
                    }
                  },
       "D.diag" = { out <- ssa.d.diag(eval_a, nu) 
                    if (suspendedTauLeapMethod) { 
                      suspendedTauLeapMethod <- suspendedTauLeapMethod-1
                      nSuspendedTauLeaps <- nSuspendedTauLeaps + 1
                      if (!suspendedTauLeapMethod) method <- "OTL.diag"
                    }
                  },
     "ETL.diag" = { out <- ssa.etl.diag(eval_a, nu, tau) }, 
     "BTL.diag" = { out <- ssa.btl.diag(eval_a, nu, x, f) },
     "OTL.diag" = { out <- ssa.otl.diag(x, eval_a, nu, hor, nc, epsilon, dtf, nd)
                    suspendedTauLeapMethod <- out$suspendedTauLeapMethod
                    if (suspendedTauLeapMethod) { 
                      method <- "D.diag"
                      doCalc <- FALSE
                    }
                  },
                  stop("unknown SSA method")
    )

    if (doCalc) {
      t <- t + out$tau  # Update the time
      x <- x + out$nu_j # Update the state vector

      # Check that no states are negative (can occur in some tau-leaping methods)
      if ((any(x<0)) & (!ignoreNegativeState)) {
        cat("at least one population in 'x' is negative. Bailing to browser...\n")
        browser()
      }

      # We need to record the step size separatelly from the resultMatrix (below)
      # since the history may not be recorded at each step (depending on the value
      # of 'censusInterval')
      stepSize <- c(stepSize, out$tau)

      # If it's time..., record the current state of the system (t,x)
      if (timeOfNextCensus <= t) { 
        timeSeries[currentRow,] <- c(t, x)
        currentRow              <- currentRow + 1
        timeOfNextCensus        <- t + censusInterval

        # If necessary add empty rows to the time series matrix
        if (currentRow > dim(timeSeries)[1]) timeSeries <- rbind(timeSeries, matrix(nrow=1000, ncol=(numCols)))
      } # if()

      # Evaluate the transition rates for the next step by evaluating the 
      # propensity functions
      for (num in seq_len(length(parse_a))) eval_a[num] <- eval(parse_a[num])
      eval_a[is.na(eval_a)] <- 0 # Replace NA with zero (0/0 gives NA)
      if(any(eval_a<0)) warning("negative propensity function - coersing to zero\n")
      eval_a[eval_a<0] <- 0
    } # if (!suspendedTauLeapMethod)
    procTimeEnd <- proc.time()
    elapsedWallTime <- procTimeEnd[3]-procTimeStart[3]
 } # while()
  #############################################################################
  # ... we are done. End the simulation loop...
  #############################################################################

  # If applicable, display the last time step on the console
  if (verbose) {
    cat("t=",t," : ",sep="")
    cat(x,sep=",")
    cat("\n")
    flush.console()
  }

  # Record the final state of the system
  timeSeries[currentRow,] <- c(t, x)
  endWallTime <- format(Sys.time())
 
  # Figure out all the reasons why the simulation terminated
  terminationStatus <- NULL 
  if (t>=tf)          terminationStatus <- c(terminationStatus, "finalTime")
  if (all(x==0))      terminationStatus <- c(terminationStatus, "extinction")
  if (any(x<0))       terminationStatus <- c(terminationStatus, "negativeState")
  if (all(eval_a==0)) terminationStatus <- c(terminationStatus, "zeroProp")
  if (elapsedWallTime>=maxWallTime) terminationStatus <- c(terminationStatus, "maxWallTime") 

  # Calculate some stats for the used method
  stats <- list(startWallime       = startWallTime,
                endWallTime        = endWallTime,
                elapsedWallTime    = elapsedWallTime,
                terminationStatus  = terminationStatus,
                nSteps             = length(stepSize),
                meanStepSize       = mean(stepSize),
                sdStepSize         = sd(stepSize),
                nSuspendedTauLeaps = nSuspendedTauLeaps)

  # Figure out why the simulation terminated and print some info/stats
  if (verbose) { 
    cat("--------------------\n")
    cat("tf: ",t,"\n","TerminationStatus: ",sep="")
    cat(stats$terminationStatus,sep=",")
    cat("\nDuration: ",stats$elapsedWallTime," seconds\n",
        "Method: ",method,"\n",
        "Nr of steps: ",stats$nSteps,"\n",
        "Mean step size: ",stats$meanStepSize,"+/-",stats$sdStepSize,"\n",sep="")
    if (method=="OTL") {
      cat("Nr suspended tau leaps: ",stats$nSuspendedTauLeaps,
          "(",100*(round(stats$nSuspendedTauLeaps/stats$nSteps)),"%)\n",sep="")
    }
    cat("End wall time: ",endWallTime,"\n",sep="")
    }

  # Return simulation results ('chopping' off any rows in the timeSeries matrix that have no values (NA))    
  return(list(data  = timeSeries[!is.na(timeSeries[,1]),], 
              stats = stats, 
              args  = args))
}
