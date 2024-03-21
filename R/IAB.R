# ---- roxygen documentation ----
#
#' @title Benhamou's IAB Index
#' 
#' @description
#' The function \code{IAB} computes the IAB index following the methods described in the paper by 
#' Benhamou et al. (2014). It facilitates global analysis, with the significance testing procedure
#' described in the paper, but also a local level output, to explore the IAB statistic through time.
#'
#' @details
#' The function \code{IAB} can be used to test for direct interaction in wildlife telemetry data and affords a novel significance testing procedure that takes into account the serially correlated structure of telemetry data. Specifically, it computes an index analogous to the Bhattacharyya coefficient between the potential influence domains of two animals. Like the other indices, IAB is dependent on the selection of an appropriate value for \code{dc} (which is termed \eqn{\Delta}{Delta} in the article). The \code{dc} parameter here is not a threshold distance, but rather the distance at which the function shows maximum slope (see Benhamou et al. 2014). 
#' 
#' The significance testing procedure uses a wrapped shifting method in order to maintain the serially correlated structure of the data. At each shift, a sample value of IAB (termed MAB) is computed in order to generate a distribution of values to test against (for more information see Benhamou et al. 2014). Here a local version of this statistical testing procedure is implemented by taking \code{rand} samples of the \eqn{(n^2 - n} permutations of unpaired fixes. The p-values are computed following Benhamou et al. (2014), z-scores are calculated based on the mean and standard deviation of this hypothetical distribution.
#'
#' @param traj an object of the class \code{move2} which contains the time-stamped movement fixes of at least two individuals. For more information on objects of this type see \code{help(mt_as_move2)}.
#' @param traj2 (optional) same as traj, but for the second group of individuals. See \code{checkTO}
#' @param tc time threshold for determining simultaneous fixes -- see function: \code{GetSimultaneous}.
#' @param dc critical distance where the IAB function will show maximum slope -- see Benhamou et al. (2014) for more advice on selecting this parameter.
#' @param local logical value indicating whether a dataframe (\code{local = TRUE}) containing the IAB index for each simultaneous fix should be returned (with a local permutation test), or (if \code{local = FALSE} - the default) the global index along with associated global permutation test.
#' @param rand number of permutations to use in the local permutation test.
#'
#' @return
#' If \code{local=FALSE} (the default) IAB returns a dataframe with the values of the IAB index and the associated p-values for one-sided tests for attraction or avoidance. If \code{local=TRUE} IAB returns a dataframe (containing the date/times of \emph{all} simultaneous fixes (NOTE: times are associated with traj1), along with the distance between fixes at each time , and the IAB index value for each simultaneous fix. A localized p-value (.pa signifies the test for attraction and pb the test for avoidance) and z-score is computed based on \code{rand} permutations of the fixes. The row.name columns can be used to match the simultaneous fixes to the original trajectory.
#'
#' @references
#' Benhamou, S., Valeix, M., Chamaille-Jammes, S., Macdonald, D., Loveridge, A.J. (2014) Movement-based analysis of interactions in African lions. \emph{Animal Behaviour}, \bold{90}: 171-180.
#'
#' @keywords indices
#' @seealso GetSimultaneous, DI, Prox
#' @examples
#' data(deer)
#' #tc = 7.5 minutes, dc = 50 meters
#' IAB(deer, tc=7.5*60, dc=50)
#' df <- IAB(deer, tc=7.5*60, dc=50, local=TRUE)
#' 
#' @export
# ---- End of roxygen documentation ----


### Simon Benhamou's IAB method for ST interaction
IAB <- function(traj,traj2,tc=0,dc=0,local=FALSE,rand=99){
  
  
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
  
  ## INTERNAL FUNCTIONS
  #-----------------------------------
  MAB <- function(k,n,diags,dc){
    De <- c(diags[[k]],diags[[k+n]])
    Ie <- exp((-1/2)*(De/dc)^2)
    Ie
  }
  #function to compute permutation p-value (one-sided attract)
  perm.pa <- function(di,di.){
    k <- length(di.)
    ng <- length(which(di. > di))
    p.above <- (ng+1)/k
    return(p.above)
  }
  #function to compute permutation p-value (one-sided avoid)
  perm.pb <- function(di,di.){
    k <- length(di.)
    nb <- length(which(di. < di))
    p.below <- (nb+1)/k
    return(p.below)
  }
  #function ot compute permutation z-score
  perm.z <- function(di,di.){
    u <- mean(di.)
    s <- sd(di.)
    z <- (di - u)/s
    return(z)
  }
  #---------------------------------------------------
  
  
  #======================================
  # Global Analysis
  #======================================
  if (local==FALSE){
    pairs$IAB.obs <- NA
    pairs$IAB.exp<- NA
    pairs$p.attract <- NA
    pairs$p.avoid <- NA
    #Compute DI for all pairs
    for (i in 1:n.pairs){
      
      traj1 <- mtraj[mt_track_id(traj)==pairs$ID1[i],]
      traj2 <- mtraj[mt_track_id(traj)==pairs$ID2[i],]
      
      trajs <- GetSimultaneous(traj1, traj2, tc)
      
      traj1 <- trajs[mt_track_id(trajs)==pairs$ID1[i],]
      traj2 <- trajs[mt_track_id(trajs)==pairs$ID2[i],]
      
      n <- nrow(traj1)
      
      dM <- st_distance(traj1,traj2)
      units(dM) <- NULL
      
      IAB.df <- data.frame(date=traj1$date,Dab=diag(dM))
      IAB.df$Iab <- exp((-1/2)*(IAB.df$Dab/dc)^2)
      
      #Significance Testing Stuff
      #-------------------------
      dI = row(dM) - col(dM)
      diags = split(dM,dI)
      kk <- 1:(n-1)
      EM <-  as.matrix(sapply(kk,MAB,n,diags,dc),nrow=n,by.row=FALSE)
      IAB. <- mean(IAB.df$Iab,na.rm=TRUE)
      EM. <- colMeans(EM)
      ng <- length(which(EM. > IAB.))
      nb <- length(which(EM. < IAB.))
      P.attract <- (ng + 1)/n
      P.avoid <- (nb + 1)/n
      #-------------------------
      
      pairs$IAB.obs[i] <- IAB.
      pairs$IAB.exp[i] <- mean(EM.)
      pairs$p.attract[i] <- P.attract
      pairs$p.avoid[i] <- P.avoid
    }
    
    return(pairs)
      
  #======================================
  # Local Analysis
  #======================================
  } else {
    return.df <- NULL
    for (i in 1:n.pairs){
      
      traj1 <- traj[mt_track_id(traj)==pairs$ID1[i],]
      traj2 <- traj[mt_track_id(traj)==pairs$ID2[i],]
      
      trajs <- GetSimultaneous(traj1, traj2, tc)
      
      traj1 <- trajs[mt_track_id(trajs)==pairs$ID1[i],]
      traj2 <- trajs[mt_track_id(trajs)==pairs$ID2[i],]
      
      n <- nrow(traj1)
      
      dM = st_distance(traj1,traj2)
      units(dM) <- NULL
      
      IAB.df <- data.frame(id1=pairs$ID1[i],id2=pairs$ID2[i],date=traj1$date,Dab=diag(dM))
      IAB.df$Iab <- exp((-1/2)*(IAB.df$Dab/dc)^2)
      
      
      #Significance Testing Stuff
      #-------------------------
      dI = row(dM) - col(dM)
      diags = split(dM,dI)
      kk <- 1:(n-1)
      EM <-  as.matrix(sapply(kk,MAB,n,diags,dc),nrow=n,by.row=FALSE)
      #check here to make sure requested number of permutations is not greater than actual number available
      rr <- length(EM)  #should be n^2 - n
      if (rand > rr){
        print(paste(rand, ' permutations were requested, but only ',rr,' are possible; rand changed to ',rr,sep=''))
        rand <- rr
      }
      
      Iab. = sample(EM,rand)
      IAB.df$Iab.pa <- sapply(IAB.df$Iab,perm.pa,Iab.)
      IAB.df$Iab.pb <- sapply(IAB.df$Iab,perm.pb,Iab.)
      IAB.df$Iab.z <- sapply(IAB.df$Iab,perm.z,Iab.)
      IAB.df$row.name1 <- row.names(traj1)
      IAB.df$row.name2 <- row.names(traj2)
      if (is.null(return.df)){
        return.df <- IAB.df
      } else{
        return.df <- rbind(return.df,IAB.df)
      }
    }
    return(return.df)
  }
  
}
