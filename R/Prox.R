# ---- roxygen documentation ----
#
#' @title Proximity Index
#'
#' @description
#' The function \code{Prox} simply computes the proportion of (simultaneous) fixes that are proximal, based on some spatial
#' threshold -- \code{dc} (Bertrand et al. 1996). It also facilitates local-level proximity analysis
#'
#' @details
#' The function \code{Prox} can be used to test for the presence of attraction (via proximity) in wildlife telemetry data. Prox is simply the proportion of simultaneous fixes within the threshold distance -- \code{dc}. The local output (dataframe) can be useful for examining variation in proximity through time.
#'
#' @param traj1 an object of the class \code{ltraj} which contains the time-stamped
#'    movement fixes of the first object. Note this object must be a \code{type II
#'    ltraj} object. For more information on objects of this type see \code{help(ltraj)}.
#' @param traj2 same as \code{traj1}.
#' @param tc time threshold for determining simultaneous fixes -- see function: \code{GetSimultaneous}.
#' @param dc distance tolerance limit (in appropriate units) for defining when 
#'    two fixes are spatially together.
#' @param local logical value indicating whether or not a dataframe, containing the
#'    distance between each simultaneous fix, should be returned.
#' @param GetSimultaneous logical value indicating whether proximity analysis is based on simultaneous fixes (if \code{TRUE} the default) -- see function \code{GetSimultaneous} or (if \code{FALSE}) a one-way mapping from \code{traj1} to {traj2} is used. 
#'
#' @return
#' If \code{local=FALSE} (the default) Prox returns the numeric value of the Prox index.
#' If \code{local=TRUE} Prox returns a dataframe containing the date/times of \emph{all} simultaneous fixes from \code{traj1}, 
#' and in the case of \code{GetSimultaneous = FALSE} the time of the fixes that were deemed simultaneous in \code{traj2}.
#' If \code{GetSimultaneous = TRUE} (the default) the Prox considers only the simultaneous fixes, as defined in \code{GetSimultaneous}. If \code{FALSE} Prox considers all the fixes in \code{traj1} relative to \code{traj2}. The latter functionality is useful when the time between fixes for one trajectory (\code{traj1}) is  much shorter than the second trajectory. 
#'
#' @references
#' Bertrand, M.R., DeNicola, A.J., Beissinger, S.R, Swihart, R.K. (1996) Effects of parturition
#' on home ranges and social affiliations of female white-tailed deer.
#' \emph{Journal of Wildlife Management}, \bold{60}: 899-909.
#'
#' @keywords indices
#' @seealso GetSimultaneous, contacts
#' @examples
#' data(deer)
#' deer37 <- deer[1]
#' deer38 <- deer[2]
#' #tc = 7.5 minutes, dc = 50 meters
#' Prox(deer37, deer38, tc=7.5*60, dc=50)
#' df <- Prox(deer37, deer38, tc=7.5*60, dc=50, local=TRUE)
#' 
#' @export
#
# ---- End of roxygen documentation ----
#================= New Prox function ============================================
Prox <- function(traj1,traj2,tc=0,dc=50,local=FALSE,GetSimultaneous=TRUE){
  
  if (GetSimultaneous==TRUE){
    
    trajs <- GetSimultaneous(traj1, traj2, tc)
    #convert ltraj objects to sf
    tr1 <- ltraj2sf(trajs[1])
    tr2 <- ltraj2sf(trajs[2])

    #Calculate the observed distances
    prox.df <- data.frame(date1=tr1$date,row1=row.names(tr1),date2=tr2$date,row2=row.names(tr2))
    prox.df$dt <- abs(difftime(prox.df$date2,prox.df$date1,units="secs"))
    prox.df$prox <- diag(st_distance(tr1,tr2))
  } else {
    ## This is just a subset of the code from GetSimultaneous.
    #store as dataframes
    #convert ltraj objects to sf
    tr1 <- ltraj2sf(traj1)
    tr2 <- ltraj2sf(traj2)
    n1 <- nrow(tr1)
    
    #Under this new system the proximity is directional so it is from traj1 to traj2
    # i.e., the order the traj's are input matters.
    # the output dataframe will have the same number of records as traj1
    prox.df <- data.frame(date1=tr1$date,row1=as.numeric(row.names(tr1)),date2=tr1$date,row2=NA,dt=NA,prox=NA)
    for (i in 1:n1){
      matched <- which.min(abs(difftime(tr1$date[i],tr2$date,units="secs")))
      prox.df$date2[i] <- tr2$date[matched]
      prox.df$row2[i] <- as.numeric(row.names(tr2)[matched])
      prox.df$dt[i] <- abs(difftime(prox.df$date2[i],prox.df$date1[i],units="secs"))
      prox.df$prox[i] <- st_distance(tr1[i,],tr2[matched,])
    }
  }
  #compute the Prox index
  n <- nrow(prox.df)
  nprox <- length(which(prox.df$prox < dc))
  val <- nprox/n
  #Return a list of the Prox index value, and a dataframe for plotting if set to true
  if (local == TRUE){
    return(prox.df)
  }
  else {return(val)}
}
#================= end of Prox2 function ========================================