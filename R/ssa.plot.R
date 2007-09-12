# $Id: ssa.plot.R 146 2007-09-12 22:04:21Z pineda $

ssa.plot <- function(out = stop("requires simulation output object"), 
             plot.device = "x11", 
                    file = "ssaplot",
                      by = 1,
               plot.from = 2,
                 plot.to = dim(out$data)[2],
                 plot.by = 1){ # number: increment of the sequence.

  if ((plot.from == 1) || (plot.from > dim(out$data)[2])|| (plot.from > plot.to)) stop("error in plot.from/plot.to arguments")
                    
  # Set the correct graphics device
  if (plot.device == "x11") x11()
  else if (plot.device=="pdf")  pdf(file=paste(file,".pdf",sep=""))
  else if (plot.device=="png")  png(file=paste(file,".png",sep=""))
  else if (plot.device=="jpeg") png(file=paste(file,".jpeg",sep=""))
  else if (plot.device=="bmp")  png(file=paste(file,".bmp",sep=""))
  else stop("unrecognized graphics device")
  
  # Render the plot(s)
  colorVector <- rainbow(dim(out$data)[2]-1)
  mask <- seq(1,dim(out$data)[1],by)
  matplot(out$data[mask,1],
          out$data[mask, seq(plot.from,plot.to,plot.by)],
          pch=19,
          cex=0.1,
          col=colorVector,
          bty="n",
          xlab="Time",
          ylab="Frequency",...)
  title(out$args$simName)
  legendTxt <- names(out$arg$x0)

  # If there are more states than 20 the legend starts to look crazy
  if (length(legendTxt) < 20) legend("topright",legend=legendTxt,bty="y",pch=19,col=colorVector) 

  if (by==1) stepShowStr <- paste(" (showing all steps)")
  else stepShowStr <- paste(" (points plotted every ",by," steps)",sep="")

  textStr <- paste("Method: ", out$args$method,", Elapsed wall time: ",round(out$stats$elapsedWallTime,2)," sec, ",out$stats$nSteps," steps",stepShowStr,sep="") 
  if (out$arg$method=="OTL") textStr <- paste(textStr,", ",out$stats$nSuspendedTauLeaps," (",round((out$stats$nSuspendedTauLeaps/out$stats$nSteps),2),"%) susp. tau-leaps",sep="")
  mtext(textStr,line=0)
  
  # Wrap up
  if (plot.device == "x11") return()
  dev.off()
  cat("ssa.plot saved as ",file,"\n",sep="")
}
