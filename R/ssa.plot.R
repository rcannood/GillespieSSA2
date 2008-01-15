# $Id: ssa.plot.R 207 2008-01-09 15:19:24Z pineda $

ssa.plot <- function(out = stop("requires simulation output object"), 
                    file = "ssaplot",
                      by = 1,
               plot.from = 2,
                 plot.to = dim(out$data)[2],
                 plot.by = 1, # number: increment of the sequence.
              show.title = TRUE,
             show.legend = TRUE){
  if ((plot.from == 1) || (plot.from > dim(out$data)[2])|| (plot.from > plot.to)) stop("error in plot.from/plot.to arguments")
                      
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
          ylab="Frequency")
  if (show.title) title(out$args$simName)
  legendTxt <- names(out$arg$x0)

  # If there are more states than 20 the legend starts to look crazy, so we don't show it...
  if (length(legendTxt) < 20 & show.legend) legend("topright",legend=legendTxt,bty="y",pch=19,col=colorVector) 

  stepShowStr <- paste("(",by," steps/point)",sep="")
  textStr <- paste(out$args$method,", ",round(out$stats$elapsedWallTime,2)," sec, ",out$stats$nSteps," steps ",stepShowStr,sep="") 
  mtext(textStr,line=0,cex=0.75)
}
