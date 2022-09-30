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
#' @param traj1 an object of the class \code{ltraj} which contains the time-stamped 
#' movement fixes of the first object. Note this object must be a \code{type II 
#' traj} object. For more information on objects of this type see \code{help(ltraj)}.
#' @param traj2 same as \code{traj1}.
#' @param tc time threshold for determining simultaneous fixes -- see function: \code{GetSimultaneous}.
#' @param local logical value indicating whether a dataframe (\code{local = TRUE}) containing the IAB 
#'    index for each simultaneous fix should be returned (with a local permutation test), or (if \code{local = FALSE} - the default) 
#'    the global index along with associated global permutation test.
#' @param rand number of permutations to use in the local permutation test.
#' @param alpha value for the \eqn{\alpha} parameter in the formula for di\eqn{_d} (default = 1).
#'
#' @return
#' If \code{local=FALSE} (the default) DI returns the numeric value of the DI index (along with DI\eqn{_{theta}}{_theta} and DI\eqn{_d}), and the associated p-value from a permutation test (see \code{IAB}).
#' If \code{local=TRUE} DI returns a dataframe that contains the localized \code{di} values (see Long and Nelson 2013). The columns for \code{di}, \code{di.theta}, and \code{di.d} represent dynamic interaction overall, in direction (azimuth), and in displacement, respectively for each segment. A localized p-value for a one sided test for positive interaction (and z-score) is computed based on \code{rand} permutations of the segments. The pkey columns can be used to match the simultaneous segments to the original trajectory (see \code{IAB}). 
#'
#' @references
#' Long, J.A., Nelson, T.A. 2013. Measuring dynamic interaction in movement data. \emph{Transactions in GIS}. 17(1): 62-77.
#'
#' @keywords indices
#' @seealso GetSimultaneous, Cr, IAB
#' @examples
#' data(deer)
#' deer37 <- deer[1]
#' deer38 <- deer[2]
#' #tc = 7.5 minutes
#' DI(deer37, deer38, tc = 7.5*60)
#' \dontrun{
#' df <- DI(deer37, deer38, tc = 7.5*60, local = TRUE)
#' }
#' 
#' @export
#
# ---- End of roxygen documentation ----
DI <- function(traj1,traj2,tc=0,local=FALSE,rand=99,alpha=1){
  #convert ltraj objects to dataframes
  trajs <- GetSimultaneous(traj1, traj2, tc)
  #convert ltraj objects to dataframes
  tr1 <- ld(trajs[1])
  tr2 <- ld(trajs[2])

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
  
  n <- nrow(tr1)           #number of fixes
  tr1 <- tr1[1:(n-1),]
  tr2 <- tr2[1:(n-1),]
  theta <- mapply(f.theta, tr1$abs.angle,tr2$abs.angle)
  disp <- mapply(f.disp, tr1$dist,tr2$dist, alpha)
  
  #compute overall interaction
  di.total <-theta*disp
  
  if (local == FALSE){
    DI.TOT <- mean(di.total)
    DI.theta <- mean(theta)
    DI.disp <- mean(disp)
    
    #Significance Testing
    #-------------------------
    theta. <- NULL
    disp. <- NULL
    DI. <- NULL
    n <- nrow(tr1)
    for (k in 1:(n-1)){
      tr2. <- rbind(tr2[(k+1):n,],tr2[1:k,])
      theta <- mapply(f.theta, tr1$abs.angle,tr2.$abs.angle)
      disp <- mapply(f.disp, tr1$dist,tr2.$dist,alpha)
      di <- theta*disp    
      #theta. <- c(theta.,mean(theta,na.rm=TRUE))
      #disp. <- c(disp., mean(disp,na.rm=TRUE))
      DI. <- c(DI.,mean(di,na.rm=TRUE))  
    }
    ng <- length(which(DI. > DI.TOT))
    nb <- length(which(DI. < DI.TOT))
    P.positive <- (ng + 1)/n
    P.negative <- (nb + 1)/n
    #-------------------------
    return(list(DI=DI.TOT,DI.theta=DI.theta,DI.d = DI.disp,P.positive=P.positive,P.negative=P.negative))
    
  } else {
    #Compute the permutations
    i <- 1:(n-1)             #n-1 segments
    df.rand <-expand.grid(i,i)
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
    
    #get the p-value and z-score (could do for both theta and disp components as well but overly complicated output)
    di.p <- sapply(di.total,perm.p,di.)
    di.z <- sapply(di.total,perm.z,di.)
    #get the z-scores
    
    outdf <- data.frame(date = tr1$date, di.theta = theta, di.d = disp, di = di.total,  di.p = di.p, di.z = di.z)
    #need to add an extra row as it is only n-1 rows
    outdf <- rbind(outdf,rep(NA,6))
    outdf[n,1] <- ld(trajs[1])$date[n]
    outdf$pkey1 <- ld(trajs[1])$pkey
    outdf$pkey2 <- ld(trajs[2])$pkey
    return(outdf)
  }
}
#==================== End of di Function =======================================
