# ---- roxygen documentation ----
#
#' @title Mapping wildlife contacts
#'
#' @description
#' The function \code{contacts} is a simple function for mapping where wildlife contacts occur on the landscape with wildilfe telemetry data. 
#' 
#' @details
#' The function \code{contacts} can be used to map where contacts occur on the lansdcape, contacts being defined spatially based on a distance threshold \code{dc} and temporally based on the time threshold \code{tc} -- see the function \code{getsimultaneous}. The location of the contact is calculated as the midpoint between the two fixes that are determined to be a "contact" based on \code{dc} and \code{tc}.  
#'
#' @param traj1 an object of the class \code{ltraj} which contains the time-stamped movement fixes of the first object. Note this object must be a \code{type II ltraj} object. For more information on objects of this type see \code{help(ltraj)}.
#' @param traj2 same as \code{traj1}.
#' @param tc time threshold for determining simultaneous fixes -- see function: \code{GetSimultaneous}.
#' @param dc distance tolerance limit (in appropriate units) for defining when two fixes are spatially together.
#' @param type whether to generate contacts as a \code{SpatialPointsDataFrame} of \code{SpatialLinesDataFrame}, points are the default, but lines can be useful for plotting and exploratory analysis.
#' @param proj4string a string object containing the projection information to be passed included in the output \code{SpatialPolygonsDataFrame} object. For more information see the \code{CRS-class} in the packages \code{sp} and \code{rgdal}. Default is \code{NA}.
#' @param tau (only used when \code{type = 'l'}), the time interval for interpolating between fixes, defaults to 1/10 the median fix interval.
#'
#' @return
#' A \code{SpatialPointsDataFrame} or \code{SpatialLinesDataFrame} containing the locations/paths of the contacts. The time of the contact is stored in the attributes of the \code{SpatialPointsDataFrame} object, along with the actual distance between fixes. The \code{SpatialLinesDataFrame} contains attributes of the time of contact, and the min, max, and mean distance apart along a line segment.
#'
#' @keywords mapping
#' @seealso GetSimultaneous, Prox
#' @examples
#' data(deer)
#' deer37 <- deer[1]
#' deer38 <- deer[2]
#' #tc = 7.5 minutes, dc = 50 meters
#' spts <- contacts(deer37, deer38, tc=7.5*60, dc=50)
#' slin <- contacts(deer37, deer38, tc=7.5*60, dc=50, type='l')
#' #plot(spts)
#' #plot(slin,add=TRUE)
#' 
#' @export
#

contacts <- function(traj1,traj2,dc=0,tc=0,type='p',proj4string=CRS(as.character(NA)),tau=NA){
  if (type=='p'){
    #convert ltraj objects to dataframes
    trajs <- GetSimultaneous(traj1, traj2, tc)
    #convert ltraj objects to dataframes
    tr1 <- ld(trajs[1])
    tr2 <- ld(trajs[2])
    trdf <- data.frame(x=(0.5*(tr1$x+tr2$x)),y=(0.5*(tr1$y+tr2$y)),date=tr1$date,prox=sqrt(((tr1$x - tr2$x)^2) + ((tr1$y - tr2$y)^2)))
    ind <- which(trdf$prox <= dc)
    spts <- SpatialPointsDataFrame(trdf[ind,1:2],data=trdf[ind,],proj4string=proj4string)
    return(spts)
  } else if (type=='l'){
    #Internal linear interpolation function
    linear <- function(traj,tau){
      
      tr <- ld(traj)
      #Function gets index of fix before time point.
      #-----------------------------------------------
      get.anchor <- function(traj,tt){
        traj.df <- ld(traj)
        indF <- max(which(difftime(traj.df$date,tt) <= 0))
        return(indF)
      }
      #-----------------------------------------------
      
      if (length(tau) > 1){
        #interpolate at user-specified times
        times <- tau
      } else {
        #otherwise interpolate by consistent interval tau
        times <- seq(min(tr$date), max(tr$date), by=tau)
      }
      
      #set-up the output dataframe.
      n <- length(times)
      out.df <- data.frame(x=rep(0,n),y=rep(0,n),date=times)
      
      for (i in 1:n){
        t <- times[i]
        ind <- get.anchor(traj,t)
        dt <- tr$date[ind+1] - tr$date[ind]
        t.slice <- t - tr$date[ind]
        rat <- as.numeric(t.slice)/as.numeric(dt)
        bx <- tr$x[ind] + rat*(tr$x[ind+1]-tr$x[ind])
        by <- tr$y[ind] + rat*(tr$y[ind+1]-tr$y[ind])
        out.df[i,1:2] <- c(bx,by) 
      }
      out.traj <- as.ltraj(out.df[,1:2],out.df$date, id=attr(traj[[1]],'id'))
      return(out.traj)
    }
    
    #Get overlapping paths
    trajs <- GetTemporalOverlap(traj1,traj2,tc)
    
    #convert to dataframes
    traj1 <- trajs[1]
    traj2 <- trajs[2]
    tr1 <- ld(traj1)
    tr2 <- ld(traj2)
    
    #If tau is not user defined default to 10 points between median fix interval.
    if (is.na(tau)) tau <- round(median(c(tr1$dt,tr2$dt),na.rm=TRUE)/10)
    
    #identify the times at which to interpolate 
    t.min <- max(c(min(tr1$date),min(tr2$date)))
    t.max <- min(c(max(tr1$date),max(tr2$date)))
    times <- seq(t.min, t.max, by=tau)
    
    #linear interpolate each trajectory at the same time points
    traj1 <- linear(traj1, times)
    traj2 <- linear(traj2, times)
    tr1 <- ld(traj1)
    tr2 <- ld(traj2)
    
    trdf <- data.frame(x=(0.5*(tr1$x+tr2$x)),y=(0.5*(tr1$y+tr2$y)),date=tr1$date,prox=sqrt(((tr1$x - tr2$x)^2) + ((tr1$y - tr2$y)^2)))
    trdf$bin <- 0
    ind <- which(trdf$prox <= dc)
    trdf$bin[ind] <- 1
    
    a <- rle(trdf$bin)
    
    val <- 0
    j <- 0
    n <- length(which(a$values == 1))
    lines.list <- list()
    ## time attributes being forced to numeric...
    out.df <- data.frame(t1=rep(0,n),t2=rep(0,n),min.prox=rep(0,n),max.prox=rep(0,n),mean.prox=rep(0,n))
    
    for (i in 1:length(a$lengths)){
      #if it is a 1 create a line object
      if (a$values[i] == 1){
        ind  <- val + 1:(a$lengths[i])
        #get coordinates of line
        xy <- trdf[ind,1:2]
        #if it is only a single point, duplicate it
        if (dim(xy)[1] == 1){xy[2,] <- xy[1,]}
        #append the line opect to the lines list
        lines.list <- c(lines.list, Lines(list(Line(xy)),as.character(j)))    #this could be slow/cause memory issues
        #grab some attribute data
        j <- j + 1
        out.df[j,] <- c(min(trdf$date[ind]),max(trdf$date[ind]),min(trdf$prox[ind]),max(trdf$prox[ind]),mean(trdf$prox[ind]))
      }
      val <- val + a$lengths[i]
    }
    
    sl <- SpatialLines(lines.list,proj4string=proj4string)
    sldf <- SpatialLinesDataFrame(sl,data=out.df,match.ID=FALSE)
    
    return(sldf)
  }

}