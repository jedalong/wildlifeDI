# ---- roxygen documentation ----
#
#' @title Contact distance plot
#'
#' @description
#' This function is an exploratory tool to examine the pairwise distances between individuals within a large telemetry dataset.
#'
#' @details
#' The dcPlot function can be used to study the frequency distribution of pairwise distances between individual in a large telemetry dataset. It can be applied to a single group (if \code{mtraj2} is ignored) or two-groups of individuals. The code attempts to find natural breaks (local minima) in the frequency histogram using an approach based on the \code{peaks} function attributed to B. Ripley (see \url{https://stackoverflow.com/questions/6324354/add-a-curve-that-fits-the-peaks-from-a-plot-in-r} ). This tool is meant to be used for exploratory data analysis.

#' @param traj an object of the class \code{move2} which contains the time-stamped movement fixes of at least two individuals. For more information on objects of this type see \code{help(mt_as_move2)}.
#' @param traj2 (optional) same as traj, but for the second group of individuals. See \code{checkTO}.
#' @param tc time threshold for determining simultaneous fixes -- see function: \code{GetSimultaneous}.
#' @param histplot (logical) whether to output a histogram, along with a list of the natural breaks in the histogram (\code{histplot = TRUE}) or the dataframe of all paired distances used to construct the histogram (\code{histplot=FALSE}) to be used for further analysis.
#' @param dmax (optional) distance value to 'cut-off' the distance histogram.
#'
#' @return
#' If \code{histplot = TRUE} a list of the natural breaks (local minima) identified from the frequency histogram and a plot of the frequency histogram.
#' If \code{histplot = FALSE} a dataframe containing all the pairwise and simultaneous distances between all individuals in the trajectory dataset.
#'
# @references
#'
#' @seealso GetSimultaneous, conProcess, Prox, Don, IAB
#' @examples
#' 
#' \dontrun{
#' data(does)
#' dcPlot(does,tc=15*60,dmax=1000)
#' }
#' 
#' @export
#
# ---- End of roxygen documentation ----
dcPlot <- function(traj,traj2,tc=0,histplot=TRUE,dmax){
  
  #Combine trajectories and identify overlap pairs
  if (missing(traj2)){
    pairs <- checkTO(traj)
    pairs <- pairs[pairs$TO==TRUE,]
    mtraj <- traj
  } else {
    pairs <- checkTO(traj,traj2)
    pairs <- pairs[pairs$TO==TRUE,]
    if (st_crs(traj2) != st_crs(traj)){
      traj2 <- st_transform(traj2,crs=st_crs(traj))
    }
    mtraj <- mt_stack(traj,traj2,track_combine='check_unique')
  }
  
  
  dfr <- NULL
  for (i in 1:nrow(pairs)){
    traj1 <- mtraj[mt_track_id(traj)==pairs$ID1[i],]
    traj2 <- mtraj[mt_track_id(traj)==pairs$ID2[i],]
    
    trajs <- GetSimultaneous(traj1,traj2,tc)
    
    tr1 <- trajs[mt_track_id(trajs)==pairs$ID1[i],]
    tr2 <- trajs[mt_track_id(trajs)==pairs$ID2[i],]
    
    if (nrow(tr1) > 0 & nrow(tr2) > 0){
      prox <- data.frame(id1=pairs$ID1[i],id2=pairs$ID2[i],prox=as.numeric(st_distance(tr1,tr2,by_element=TRUE)))
      if (is.null(dfr)) {dfr <- prox} else {dfr <- rbind(dfr,prox)}
    }
  }
  
  #produce plot or return dataframe
  
  if (histplot){
    ###make dc plot here
    if (missing(dmax)) dmax <- max(dfr$prox,na.rm=TRUE)

    dfr <- subset(dfr,prox <= dmax)
    breaks = seq(0, max(dfr$prox), length.out=100) 
    
    cut = cut(dfr$prox, breaks, right=FALSE) 
    freq = table(cut)
    cumfreq0 = c(0, cumsum(freq))
    cumfreq0 <- cumfreq0 / max(cumfreq0)
    
    par(mar=c(4,4,1,4))
    
    hist <- hist(dfr$prox, breaks=breaks,col='grey',border='grey',axes=T, ann=T,prob=F,
                 xlim=c(0,dmax),xlab='Distance',ylab='Frequency',main='')
    par(new=T)
    plot(breaks,cumfreq0,type='l',axes=F,ann=F)
    axis(4, at=c(0, 0.5, 1))
    mtext('Cumulative Frequency',4,line=3)
    box()
    
    peaks<-function(series,span){
      z <- embed(series, span)
      s <- span%/%2
      v<- max.col(z) == 1 + s
      result <- c(rep(FALSE,s),v)
      result <- result[1:(length(result)-s)]
      result
    } 
    
    y <- -1*hist$counts
    pks <- hist$mids[which(peaks(y, span=11))] #span should be odd numbered
    pks. <- round(pks)
    abline(v=pks,lty=3,col='red')
    pky <- seq(0.1,0.7,length.out=length(pks))
    text(x=pks,pky, labels=pks.,pos=4,offset=0.1)
    return(pks)
  } else {
    return(dfr)
  }
}
  