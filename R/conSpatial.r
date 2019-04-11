# ---- roxygen documentation ----
#
#' @title Mapping wildlife contacts
#'
#' @description
#' The function \code{contacts} is a simple function for mapping where wildlife contacts occur on the landscape with wildilfe telemetry data. 
#' 
#' @details
#' The function \code{conSpatial} can be used to map where contacts occur on the lansdcape, contacts being defined spatially based on a distance threshold \code{dc} and temporally based on the time threshold \code{tc} -- see the function \code{getsimultaneous}. The location of the contact can be calculated in a number of ways, and represented as points for each contact, or as line grouped by the contact phases. Which contacts to map can be defined in a number of ways using the \code{def} parameter: \cr
#'  \cr i) \code{def = 'all'} (the default) all fixes where column \code{contacts = 1} are returned in the Spatial* object;
#'  \cr ii) \code{def = 'phase'} all fixes which are part of a phase are returned, note the number of points when \code{def = 'phase'} should be greater than or equal to that when \code{def = 'all'} because of how phases are defined;
#'  \cr iii) \code{def = 'first'} the first location fix of each phase is returned;
#'  \cr iv) \code{def = 'last'} the last location fix of each phase is returned;
#'  \cr v) \code{def = 'minDist'} the location fix of each phase which has the minimal contact distance is returned;
#'  \cr vi) \code{def = 'minTime'} the location fix of each phase with the minimal time difference with contact fixes is returned;
#'
#' @param ltraj an object of the class \code{ltraj} which should be output from the function \code{conPhase}.
#' @param type one of ('p' - the default or 'l'). Whether to generate contacts as a \code{SpatialPointsDataFrame} or phases as a \code{SpatialLinesDataFrame}, points are the default, but lines can be useful for plotting and exploratory analysis.
#' @param def if type = 'p' one of ('all','phase','first','last','minDist','minTime') which defines how contacts are to be mapped using all or part of a contact phase. (see Details) 
#' 
#' @return
#' A \code{SpatialPointsDataFrame} or \code{SpatialLinesDataFrame} containing the locations/paths of the contacts. The time of the contact is stored in the attributes of the \code{SpatialPointsDataFrame} object, along with the actual distance between fixes. The \code{SpatialLinesDataFrame} contains attributes of the time of contact, and the min, max, and mean distance apart along a line segment.
#'
#' @keywords Contact Analysis
#' @seealso conProcess, conPhase
#' 
#' @examples 
#' 
#' #' \dontrun{
#' data(does)
#' doecons <- conProcess(does,tc=15*60,dc=50)
#' doephas <- conPhase(doecons,pc=60*60)
#' 
#' pts <- conSpatial(doephas)
#' lns <- conSpatial(doephas,type='l')
#' plot(pts)
#' plot(lns, add=TRUE)
#' }
#' 
#' @export
#

conSpatial <- function(ltraj,type='p',def='all'){
  
  proj4string <- attr(ltraj,'proj4string')
  #Function to extract minTime and minDist from phases.
  funPhase <- function(phase, df, def){
    ind <- which(df$contact_pha == phase)
    if (def=='first'){
      i1 <- ind[which.min(df$date[ind])]    
    } else if (def=='last'){
      i1 <- ind[which.max(df$date[ind])]
    } else if (def=='minTime'){
      sub <- df[df$contact_pha == phase,]
      sub$id <- as.character(sub$id)
      sub$burst <- as.character(sub$burst)
      sub <- dl(sub)
      dfpairs <- conPairs(sub)
      i1 <- ind[dfpairs$contact_orig_rowid[which.min(dfpairs$contact_dt)]]
    } else if (def=='minDist'){
      sub <- df[df$contact_pha == phase,]
      sub$id <- as.character(sub$id)
      sub$burst <- as.character(sub$burst)
      sub <- dl(sub)
      dfpairs <- conPairs(sub)
      i1 <- ind[dfpairs$contact_orig_rowid[which.min(dfpairs$contact_d)]]
    } 
    return(i1)
  }
  
  df <- ld(ltraj)
  if (type ==  'p'){
    #convert ltraj object to dataframe
    if (def =='all'){
      dfs <- df[df$contacts==1,]
    } else if (def == 'phase'){
      dfs <- df[!is.na(df$contact_pha),]
    } else {
      phaid <- unique(df$contact_pha[!is.na(df$contact_pha)])
      #Get the contacts, and optionally the BefAft Phases
      indo <- sapply(phaid,funPhase,df=df,def=def)
      dfs <- df[indo,]
    }
    spo <- SpatialPointsDataFrame(dfs[,1:2],data=dfs,proj4string=proj4string)
  } else if (type=='l'){
    df2 <- df[!is.na(df$contact_pha),]
    df2$id<- as.character(df2$id)
    df2$burst <- as.character(df2$burst)
    df2$oid <- df2$id
    df2$id <- df2$contact_pha
    df2$burst <- df2$contact_pha
    
    t2 <- dl(df2,proj4string=proj4string)
    spo <- ltraj2sldf(t2)
    names(spo@data) <- c('id','contact_pha')
  }
  return(spo)
}



# # OLD LINEAR INTERPOLATION CODE?
# #Internal linear interpolation function
# linear <- function(traj,tau){
#   #Function gets index of fix before time point.
#   #-----------------------------------------------
#   get.anchor <- function(traj,tt){
#     traj.df <- ld(traj)
#     indF <- max(which(difftime(traj.df$date,tt) <= 0))
#     return(indF)
#   }
#   #-----------------------------------------------
#   
#   if (length(tau) > 1){
#     #interpolate at user-specified times
#     times <- tau
#   } else {
#     #otherwise interpolate by consistent interval tau
#     times <- seq(min(tr$date), max(tr$date), by=tau)
#   }
#   
#   #set-up the output dataframe.
#   n <- length(times)
#   out.df <- data.frame(x=rep(0,n),y=rep(0,n),date=times)
#   
#   for (i in 1:n){
#     t <- times[i]
#     ind <- get.anchor(traj,t)
#     dt <- tr$date[ind+1] - tr$date[ind]
#     t.slice <- t - tr$date[ind]
#     rat <- as.numeric(t.slice)/as.numeric(dt)
#     bx <- tr$x[ind] + rat*(tr$x[ind+1]-tr$x[ind])
#     by <- tr$y[ind] + rat*(tr$y[ind+1]-tr$y[ind])
#     out.df[i,1:2] <- c(bx,by) 
#   }
#   out.traj <- as.ltraj(out.df[,1:2],out.df$date, id=attr(traj[[1]],'id'))
#   return(out.traj)
# }
# 
# #Get overlapping paths
# trajs <- GetTemporalOverlap(traj1,traj2,tc)
# 
# #convert to dataframes
# traj1 <- trajs[1]
# traj2 <- trajs[2]
# tr1 <- ld(traj1)
# tr2 <- ld(traj2)
# 
# #If tau is not user defined default to 10 points between median fix interval.
# if (is.na(tau)) tau <- round(median(c(tr1$dt,tr2$dt),na.rm=TRUE)/10)
# 
# #identify the times at which to interpolate 
# t.min <- max(c(min(tr1$date),min(tr2$date)))
# t.max <- min(c(max(tr1$date),max(tr2$date)))
# times <- seq(t.min, t.max, by=tau)
# 
# #linear interpolate each trajectory at the same time points
# traj1 <- linear(traj1, times)
# traj2 <- linear(traj2, times)
# tr1 <- ld(traj1)
# tr2 <- ld(traj2)
# 
# trdf <- data.frame(x=(0.5*(tr1$x+tr2$x)),y=(0.5*(tr1$y+tr2$y)),date=tr1$date,prox=sqrt(((tr1$x - tr2$x)^2) + ((tr1$y - tr2$y)^2)))
# trdf$bin <- 0
# ind <- which(trdf$prox <= dc)
# trdf$bin[ind] <- 1
# 
# a <- rle(trdf$bin)
# 
# val <- 0
# j <- 0
# n <- length(which(a$values == 1))
# lines.list <- list()
# ## time attributes being forced to numeric...
# out.df <- data.frame(t1=rep(0,n),t2=rep(0,n),min.prox=rep(0,n),max.prox=rep(0,n),mean.prox=rep(0,n))
# 
# for (i in 1:length(a$lengths)){
#   #if it is a 1 create a line object
#   if (a$values[i] == 1){
#     ind  <- val + 1:(a$lengths[i])
#     #get coordinates of line
#     xy <- trdf[ind,1:2]
#     #if it is only a single point, duplicate it
#     if (dim(xy)[1] == 1){xy[2,] <- xy[1,]}
#     #append the line opect to the lines list
#     lines.list <- c(lines.list, Lines(list(Line(xy)),as.character(j)))    #this could be slow/cause memory issues
#     #grab some attribute data
#     j <- j + 1
#     out.df[j,] <- c(min(trdf$date[ind]),max(trdf$date[ind]),min(trdf$prox[ind]),max(trdf$prox[ind]),mean(trdf$prox[ind]))
#   }
#   val <- val + a$lengths[i]
# }
# 
# sl <- SpatialLines(lines.list,proj4string=proj4string)
# sldf <- SpatialLinesDataFrame(sl,data=out.df,match.ID=FALSE)
# 
# return(sldf)