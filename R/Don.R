# ---- roxygen documentation ----
#
#' @title Doncaster's measure of dynamic interaction
#'
#' @description
#' The function \code{Don} measures the dynamic interaction between two moving objects following
#' the methods outlined by Doncaster (1990).
#'
#' @details
#' This function can be used to compute the Doncaster (1990) methods for measuring
#' dynamic interaction between two objects. The Doncaster method tests the proportion
#' of simultaneous fixes that are below \code{dc} against that which would be
#' expected based on the distribution of distances between all fixes.
#'
#' @param traj an object of the class \code{move2} which contains the time-stamped movement fixes of at least two individuals. For more information on objects of this type see \code{help(mt_as_move2)}.
#' @param traj2 (optional) same as traj, but for the second group of individuals. See \code{checkTO}
#' @param tc time threshold for determining simultaneous fixes -- see function: \code{GetSimultaneous}.
#' @param dc distance tolerance limit (in appropriate units) for defining when two fixes are spatially together.
#' @param plot logical, whether or not to plot the Doncaster plot. Default = TRUE.
#'
#' @return
#' A data.frame is returned that contains the values for the contingency table of simultaneous fixes (paired) and non-paired fixes below and above \code{dc}, along with the associated \emph{p}-value from the Chi-squared test.This function can optionally return a plot, for distance values ranging from 0 to the maximum distance separating two fixes, of the observed proportion of simultaneous fixes below each distance value (for each pair). The expected values based on all fixes are also included. 
#'
#' @references
#' Doncaster, C.P. (1992) Non-parametric estimates of interaction from radio-tracking
#' data. \emph{Journal of Theoretical Biology}, \bold{143}: 431-443.
#'
#' @keywords indices
#' @seealso GetSimultaneous
#' @examples
#' data(deer)
#' #tc = 7.5 minutes, dc = 50 meters
#' Don(deer, tc = 7.5*60, dc = 50)
#' @export
#
# ---- End of roxygen documentation ----
Don <- function(traj,traj2,tc=0,dc=0,plot=TRUE){
  
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
  
  n.pairs <- nrow(pairs)
  pairs$paired_lt_dc <- NA
  pairs$paired_gt_dc <- NA
  pairs$unpaired_lt_dc <- NA
  pairs$unpaired_gt_dc <- NA
  pairs$p_value <- NA
  
  plot.df <- NULL
  
  for (i in 1:n.pairs){
    traj1 <- mtraj[mt_track_id(traj)==pairs$ID1[i],]
    traj2 <- mtraj[mt_track_id(traj)==pairs$ID2[i],]
    A <- nrow(traj1)
    B <- nrow(traj2)
    
    trajs <- GetSimultaneous(traj1,traj2,tc)
    
    traj1 <- trajs[mt_track_id(trajs)==pairs$ID1[i],]
    traj2 <- trajs[mt_track_id(trajs)==pairs$ID2[i],]
    n <- nrow(traj1)
    #calculate the observed and expected distances
    De <- st_distance(traj1,traj2)
    Do <- diag(De)
    
    diag(De) <- NA
    
    De <- as.numeric(De)
    Do <- as.numeric(Do)
    
    pairs$paired_lt_dc[i] <- length(which(Do < dc))
    pairs$paired_gt_dc[i] <- length(which(Do >= dc))
    pairs$unpaired_lt_dc[i] <- length(which(De < dc))
    pairs$unpaired_gt_dc[i] <- length(which(De >= dc))
    con <- matrix(as.numeric(pairs[i,6:9]),nrow=2,byrow=T)
    pairs$p_value[i] <- chisq.test(con)$p.value
    
    if (plot){
      #compute the cumulative frequency plot
      pd <- seq(0,max(c(Do,De),na.rm=T),length.out=50)
      Co <- rep(0,length(pd))
      Ce <- Co
      for (j in 1:length(pd)){
        Co[j] <- length(which(Do < pd[j]))
        Ce[j] <- length(which(De < pd[j])) + Co[j]
      }
      Co <- Co / n
      Ce <- Ce / (n^2)
      
      temp.df <- data.frame(id1=pairs$ID1[i],id2=pairs$ID2[i],pd=pd,Ce=Ce,Co=Co)
      
      if (is.null(plot.df)){
        plot.df <- temp.df
      } else {
        plot.df <- rbind(plot.df,temp.df)
      }
    }

  }
  
  if (plot){
    pfun <- function(pdf) {
      mlab <- paste0(pdf$id1[1],' - ',pdf$id2[1])
      plot(pdf$pd,pdf$Ce,type="l",col='grey',ylab="probability",xlab="Distance <= (m)",main=mlab)
      points(pdf$pd,pdf$Co,pch=20)
    }
    
    x.mar <- round(sqrt(n.pairs))
    y.mar <- ceiling(sqrt(n.pairs))
    plot.df$pairID <- factor(paste0(plot.df$id1,'_',plot.df$id2))
    
    #hackjob facet wrap with base R
    par(mfcol=c(x.mar,y.mar))
    for (id in unique(plot.df$pairID)){
      pdf <- plot.df[plot.df$pairID == id,]
      pfun(pdf)
    }
  }

  
  #print(paste("A critical distance of ",dc," was used, resulting in a p-value = ", chi.p,".",sep=""))
  return(pairs)
}
#======================== End of Doncaster Function ============================