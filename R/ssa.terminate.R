 ssa.terminate <- function(args,out.rxn,tf,method,maxWallTime,verbose) {

  # Get the final time and state vector
  t <- out.rxn$timeSeries[dim(out.rxn$timeSeries)[1],1]
  x <- out.rxn$timeSeries[dim(out.rxn$timeSeries)[1],2:dim(out.rxn$timeSeries)[2]]

  # Figure out all the reasons why the simulation terminated
  terminationStatus <- NULL 
  if (t>=tf)          terminationStatus <- c(terminationStatus, "finalTime")
  if (all(x==0))      terminationStatus <- c(terminationStatus, "extinction")
  if (any(x<0))       terminationStatus <- c(terminationStatus, "negativeState")
  if (all(out.rxn$eval_a==0)) terminationStatus <- c(terminationStatus, "zeroProp")
  if (out.rxn$elapsedWallTime>=maxWallTime) terminationStatus <- c(terminationStatus, "maxWallTime") 

  # Calculate some stats for the used method
  stats <- list(startWallime       = out.rxn$startWallTime,
                endWallTime        = out.rxn$endWallTime,
                elapsedWallTime    = out.rxn$elapsedWallTime,
                terminationStatus  = terminationStatus,
                nSteps             = length(out.rxn$stepSize),
                meanStepSize       = mean(out.rxn$stepSize),
                sdStepSize         = sd(out.rxn$stepSize),
                nSuspendedTauLeaps = out.rxn$nSuspendedTauLeaps)

  # Figure out why the simulation terminated and print some info/stats
  if (verbose) {
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
    cat("End wall time: ",stats$endWallTime,"\n",sep="")
    cat("--------------------\n")
    }

  # Return simulation results ('chopping' off any rows in the timeSeries matrix that have no values (NA))    
  return(list(data  = out.rxn$timeSeries, 
              stats = stats, 
              args  = args))
}