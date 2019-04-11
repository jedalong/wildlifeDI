# ---- roxygen documentation ----
#
#' @title Get period where two tracks overlap
#'
#' @description
#'   The function \code{GetTemporalOverlap} identifies and extracts parts of a trajectory that overlap in time with another trajectory.
#'
#' @details
#'   This function is used to determine the fixes that overlap in time between two trajectories.
#'
#' @param traj1 an object of the class \code{ltraj} which contains the time-stamped
#'    movement fixes of the first object. Note this object must be a \code{type II
#'    ltraj} object. For more information on objects of this type see \code{help(ltraj)}.
#' @param traj2 same as \code{traj1}.
#' @param tc time threshold for considering if fixes are in the overlap period.
#'
#' @return
#' A single ltraj object containing two bursts, representing the two original \code{ltraj} 
#' objects, each containing only those fixes that are overlap in time.
#'
# @references
#'
#' @keywords processing
#' @seealso checkTO
#' @examples
#' data(deer)
#' deer37 <- deer[1]
#' deer38 <- deer[2]
#' trajs <- GetTemporalOverlap(deer37, deer38)
#' deer37 <- trajs[1]
#' deer38 <- trajs[2]
#' 
#' @export
#
# ---- End of roxygen documentation ----

GetTO <- function(traj1,traj2,tc=0){
  #convert to dataframe
  tr1 <- ld(traj1)
  tr2 <- ld(traj2)
  
  #identify the temporal 
  t.min <- max(c(min(tr1$date),min(tr2$date))) - tc
  t.max <- min(c(max(tr1$date),max(tr2$date))) + tc
  
  #check to see if there is an overlap period
  if (t.min > t.max){stop('There is no temporal overlap.')}
  
  #find the fixes that are in the temporal overlap
  ind1 <- which(tr1$date >= t.min & tr1$date <= t.max)
  ind2 <- which(tr2$date >= t.min & tr2$date <= t.max)
  
  #convert to ltraj objects
  tr1.to <- dl(tr1[ind1,])
  tr2.to <- dl(tr2[ind2,])
  
  #Return the two ltraj objects
  return(c(tr1.to,tr2.to))
}

