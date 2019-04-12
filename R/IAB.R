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
#' @param traj1 an object of the class \code{ltraj} which contains the time-stamped
#'    movement fixes of the first object. Note this object must be a \code{type II
#'    ltraj} object. For more information on objects of this type see \code{help(ltraj)}.
#' @param traj2 same as \code{traj1}.
#' @param tc time threshold for determining simultaneous fixes -- see function: \code{GetSimultaneous}.
#' @param dc critical distance where the IAB function will show maximum slope -- see Benhamou et al. (2014) 
#'    for more advice on selecting this parameter.
#' @param local logical value indicating whether a dataframe (\code{local = TRUE}) containing the IAB 
#'    index for each simultaneous fix should be returned (with a local permutation test), or (if \code{local = FALSE} - the default) 
#'    the global index along with associated global permutation test.
#' @param rand number of permutations to use in the local permutation test.
#'
#' @return
#' If \code{local=FALSE} (the default) IAB returns the numeric value of the IAB index and the associated p-values for one-sided tests for attraction or avoidance.
#' If \code{local=TRUE} IAB returns a dataframe (containing the date/times of \emph{all} simultaneous fixes (NOTE: times are associated with traj1), along with the distance between fixes at each time , and the IAB index value for each simultaneous fix. A localized p-value (.pa signifies the test for attraction and pb the test for avoidance) and z-score is computed based on \code{rand} permutations of the fixes. The pkey columns can be used to match the simultaneous fixes to the original trajectory.
#'
#' @references
#' Benhamou, S., Valeix, M., Chamaille-Jammes, S., Macdonald, D., Loveridge, A.J. (2014) Movement-based analysis of interactions in African lions. \emph{Animal Behaviour}, \bold{90}: 171-180.
#'
#' @keywords indices
#' @seealso GetSimultaneous, DI, Prox
#' @examples
#' data(deer)
#' deer37 <- deer[1]
#' deer38 <- deer[2]
#' #tc = 7.5 minutes, dc = 50 meters
#' IAB(deer37, deer38, tc=7.5*60, dc=50)
#' df <- IAB(deer37, deer38, tc=7.5*60, dc=50, local=TRUE)
#' 
#' @export
# ---- End of roxygen documentation ----


### Simon Benhamou's IAB method for ST interaction
IAB <- function(traj1,traj2,tc=0,dc=50,local=FALSE,rand=99){
  trajs <- GetSimultaneous(traj1, traj2, tc)
  #convert ltraj objects to dataframes
  tr1 <- ld(trajs[1])
  tr2 <- ld(trajs[2])
  n <- nrow(tr1)
  #Calculate the observed distances
  fIAB <- function(tr1,tr2,dc){
    df <- data.frame(date=tr1$date,Dab=sqrt((tr1$x - tr2$x)^2 + (tr1$y - tr2$y)^2))
    df$Iab <- exp((-1/2)*(df$Dab/dc)^2)
    return(df)
  }
  IAB.df <- fIAB(tr1,tr2,dc=dc)
  
  
  if (local == FALSE){
    IAB. <- mean(IAB.df$Iab,na.rm=TRUE)
    
    #Significance Testing
    #-------------------------
    MAB <- NULL
    for (k in 1:(n-1)){
      tr2. <- rbind(tr2[(k+1):n,],tr2[1:k,])
      MAB.df <- data.frame(Dab=sqrt((tr1$x - tr2.$x)^2 + (tr1$y - tr2.$y)^2))
      MAB.df$Mab <- exp((-1/2)*(MAB.df$Dab/dc)^2)
      MAB <- c(MAB,mean(MAB.df$Mab,na.rm=TRUE))
    }
    ng <- length(which(MAB > IAB.))
    nb <- length(which(MAB < IAB.))
    P.attract <- (ng + 1)/n
    P.avoid <- (nb + 1)/n
    #-------------------------
    
    return(list(IAB.obs=IAB.,IAB.exp=mean(MAB), P.attract=P.attract,P.avoid=P.avoid))
    
  } else {
    #Compute the permutations
    i <- 1:n             #n fixes
    df.rand <-expand.grid(i,i)
    names(df.rand) <- c('i','j')
    df.rand <- df.rand[which(df.rand$i != df.rand$j),] #remove the fixes with no shift!

    #check here to make sure requested number of permutations is not greater than actual number available
    rr <- dim(df.rand)[1]   #should be n^2 - n
    if (rand > rr){
      print(paste(rand, ' permutations were requested, but only ',rr,' are possible; rand changed to ',rr,sep=''))
      rand <- rr
    }
    df.rand <- df.rand[sample(1:rr,rand),]
    perm1 <- tr1[df.rand$i,]
    perm2 <- tr2[df.rand$j,]
    
    #These are the distributions for testing based on the 'rand' number of permutations
    IAB. <- fIAB(perm1,perm2,dc=dc)$Iab
    
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
    
    IAB.df$Iab.pa <- sapply(IAB.df$Iab,perm.pa,IAB.)
    IAB.df$Iab.pb <- sapply(IAB.df$Iab,perm.pb,IAB.)
    IAB.df$Iab.z <- sapply(IAB.df$Iab,perm.z,IAB.)
    IAB.df$pkey1 <- ld(trajs[1])$pkey
    IAB.df$pkey2 <- ld(trajs[2])$pkey
    return(IAB.df)
  }
}



