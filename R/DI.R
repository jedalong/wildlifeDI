# ---- roxygen documentation ----
#
#' @title Dynamic interaction index
#'
#' @description
#' The function \code{DI} measures dynamic interaction between two moving objects. It
#' calculates the local level di statistic for movement displacement, direction, 
#' and overall. DI can compute time- and/or distance-based weighting schemes
#' following Long and Nelson (2013).
#' 
#' @details
#' This function can be used for calculating the dynamic interaction (DI) statistic 
#' as described in Long and Nelson (2013). The DI statistic can be used to measure 
#' the local level of dynamic interaction between two moving objects. Specifically, 
#' it measures dynamic interaction in movement direction and displacement. 
#'
#' @param traj an object of the class \code{move2} which contains the time-stamped movement fixes of at least two individuals. For more information on objects of this type see \code{help(mt_as_move2)}.
#' @param traj2 (optional) same as traj, but for the second group of individuals. See \code{checkTO}
#' @param tc time threshold for determining simultaneous fixes -- see function: \code{GetSimultaneous}.
#' @param local logical value indicating whether a dataframe (\code{local = TRUE}) containing the DI index for each simultaneous fix should be returned (with a local permutation test), or (if \code{local = FALSE} - the default) the global index along with associated global permutation test.
#' @param rand number of permutations to use in the local permutation test.
#' @param alpha value for the \eqn{\alpha} parameter in the formula for di\eqn{_d} (default = 1).
#'
#' @return
#' If \code{local=FALSE} (the default) DI returns the numeric value of the DI index (along with DI\eqn{_{theta}}{_theta} and DI\eqn{_d}), and the associated p-value from a permutation test (see \code{IAB}).
#' If \code{local=TRUE} DI returns a large dataframe that contains the localized \code{di} values as a column (see Long and Nelson 2013). The columns for \code{di}, \code{di.theta}, and \code{di.d} represent dynamic interaction overall, in direction (azimuth), and in displacement, respectively for each segment. A localized p-value for a one sided test for positive interaction (and z-score) is computed based on \code{rand} permutations of the segments. The row.name columns can be used to match the simultaneous segments to the original trajectory (see \code{IAB}). 
#'
#' @references
#' Long, J.A., Nelson, T.A. 2013. Measuring dynamic interaction in movement data. \emph{Transactions in GIS}. 17(1): 62-77.
#'
#' @keywords indices
#' @seealso GetSimultaneous, Cr, IAB
#' @examples
#' \dontrun{
#' data(deer)
#' #tc = 7.5 minutes
#' DI(deer, tc = 7.5*60)
#' df <- DI(deer, tc = 7.5*60, local = TRUE)
#' }
#' 
#' @export
#
# ---- End of roxygen documentation ----
DI <- function(traj,traj2,tc=0,local=FALSE,rand=0,alpha=1){
  
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
  
  #----- INTERNAL FUNCTIONS FOR DI ----------
  #Interaction in azimuth function
  f.theta <- function(a,b){
    di.t <- cos(a-b)
    if (is.na(di.t) == T){
      if (is.na(a)==T & is.na(b)==T){
        di.t <- 1
      } else {
        di.t <- 0
      }
    }
    return(di.t)
  }
  #interaction in displacement function
  f.disp <- function(a,b,alpha){
    di.d <- 1 - (abs(a - b)/(a + b))^alpha
    if (is.na(di.d)==TRUE){
      di.d <- 1
    } 
    return(di.d)
  }
  
  #function to compute permutation p-value
  perm.p <- function(di,di.){
    k <- length(di.)
    ng <- length(which(di. > di))
    nb <- length(which(di. < di))
    p.above <- (ng+1)/k
    p.below <- (nb+1)/k
    return(p.above)
  }
  #function ot compute permutation z-score
  perm.z <- function(di,di.){
    u <- mean(di.)
    s <- sd(di.)
    z <- (di - u)/s
    return(z)
  }
  #-------------------------------------
  
  #========================
  #GLOBAL DI ANALYSIS
  #========================
  pairs$DI <- NA
  pairs$DI.theta <- NA
  pairs$DI.d <- NA
  if (rand > 0){
    pairs$P.positive <- NA
    pairs$P.negative <- NA
  }
  
  #output dataframe for local analysis 
  if (local) { return.df <- NULL }
  #Compute DI for all pairs
  for (i in 1:n.pairs){
    
    traj1 <- mtraj[mt_track_id(traj)==pairs$ID1[i],]
    traj2 <- mtraj[mt_track_id(traj)==pairs$ID2[i],]
    
    trajs <- GetSimultaneous(traj1, traj2, tc)
    
    traj1 <- trajs[mt_track_id(trajs)==pairs$ID1[i],]
    traj2 <- trajs[mt_track_id(trajs)==pairs$ID2[i],]
    
    #use as.numeric to speed up cosine function
    traj1$abs.angle <- as.numeric(mt_azimuth(traj1))
    traj2$abs.angle <- as.numeric(mt_azimuth(traj2))
    
    traj1$dist <- as.numeric(mt_distance(traj1))
    traj2$dist <- as.numeric(mt_distance(traj2))
    
    nfix <- nrow(traj1)           #number of fixes
    tr1 <- traj1[1:(nfix-1),]
    tr2 <- traj2[1:(nfix-1),]
    theta <- mapply(f.theta, tr1$abs.angle,tr2$abs.angle)
    disp <- mapply(f.disp, tr1$dist,tr2$dist, alpha)
    
    #compute overall interaction
    di.total <-theta*disp
    
    DI.TOT <- mean(di.total,na.rm=TRUE)
    DI.theta <- mean(theta,na.rm=TRUE)
    DI.disp <- mean(disp,na.rm=TRUE)
    
    pairs$DI[i] <- DI.TOT
    pairs$DI.theta[i] <- DI.theta
    pairs$DI.d[i] <- DI.disp
    
    #Significance Testing
    #-------------------------
    if (rand > 0 & local == FALSE){
      theta. <- NULL
      disp. <- NULL
      n <- nrow(tr1)
      DI. <- rep(NA,n)
      tr2. <- data.frame(abs.angle=tr2$abs.angle,dist=tr2$dist)
      kk <- 1:(n-1) |>
        sample(size=rand)
      
      for (k in kk){
        tr2_perm <- rbind(tr2.[(k+1):n,],tr2.[1:k,])
        theta <- mapply(f.theta, tr1$abs.angle,tr2_perm$abs.angle)
        disp <- mapply(f.disp, tr1$dist,tr2_perm$dist,alpha)
        di <- theta*disp    
        #theta. <- c(theta.,mean(theta,na.rm=TRUE))
        #disp. <- c(disp., mean(disp,na.rm=TRUE))
        DI.[k] <- mean(di,na.rm=TRUE)
      }
      ng <- length(which(DI. > DI.TOT))
      nb <- length(which(DI. < DI.TOT))
      P.positive <- (ng + 1)/n
      P.negative <- (nb + 1)/n
      #-------------------------
      
      pairs$P.positive[i] <- P.positive
      pairs$P.negative[i] <- P.negative
    }
    
    if (local) {
      outdf <- data.frame(id1 = pairs$ID1[i], id2 = pairs$ID2[i], date = mt_time(tr1), di.theta = theta, di.d = disp,di = di.total)
      
      # Local Significance Testing
      if (rand > 0){
        #Compute the permutations
        m <- nrow(tr1)            #n-1 segments
        df.rand <-expand.grid(m,m)
        names(df.rand) <- c('i','j')
        df.rand <- df.rand[which(df.rand$i != df.rand$j),] #remove the segments with no shift!
        
        #check here to make sure requested number of permutations is not greater than actual number available
        rr <- dim(df.rand)[1]   #should be (n-1)^2 - (n-1)
        if (rand > rr){
          print(paste(rand, ' permutations were requested, but only ',rr,' are possible; rand changed to ',rr,sep=''))
          rand <- rr
        }
        df.rand <- df.rand[sample(1:rr,rand),]
        perm1 <- tr1[df.rand$i,]
        perm2 <- tr2[df.rand$j,]
        
        #These are the distributions for testing based on the 'rand' number of permutations
        theta. <- mapply(f.theta, perm1$abs.angle,perm2$abs.angle)
        disp. <- mapply(f.disp, perm1$dist,perm2$dist,alpha)
        di. <- theta. * disp.
        #get the p-value and z-score (could do for both theta and disp components as well but overly complicated output)
        di.p <- sapply(di.total,perm.p,di.)
        di.z <- sapply(di.total,perm.z,di.)
        
        outdf <- cbind(outdf, data.frame(di.p = di.p, di.z = di.z))
        
      }
      ##ADD ROW NAMES
      outdf <- cbind(outdf,data.frame(row.name1 = row.names(tr1),row.name2 = row.names(tr2)))
      
      if (is.null(return.df)){
        return.df <- outdf
      } else{
        return.df <- rbind(return.df,outdf)
      }
    }
  }
  
  if (local) {
    return(return.df)
  } else {
    return(pairs)
  }
  
}

#==================== End of di Function =======================================
