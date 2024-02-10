# ---- roxygen documentation ----
#
#' @title Coefficient of Association
#'
#' @description
#'  This function measures the dynamic interaction between two moving objects following
#'  the methods first described by Cole (1949), and more recently employed by Bauman (1998).
#'
#' @details
#'  This function can be used to calculate the Cole (1949) measure of dynamic
#'  interaction between two animals. Termed a coefficient of association, the Ca
#'  statistic tests the number of fixes the animals are observed together against the
#'  total number of fixes following:
#'  \deqn{Ca = \frac{2AB}{A+B}}{2AB/(A+B)}
#'  where \eqn{A} (respectively \eqn{B}) is the number of times animal 1 (resp. 2) are
#'  observed, and \eqn{AB} is the number of times the two animals are observed together.
#'  Several works, including Bauman (1998) have suggested that Ca > 0.5 indicates
#'  affiliation or fidelity, while Ca < 0.5 indicates no association between the
#'  two animals. Note that this function calls \code{GetSimultaneous} to identify the temporal
#'  component of identifying when fixes together.
#'
#' @param traj an object of the class \code{move2} which contains the time-stamped movement fixes of at least two individuals. For more information on objects of this type see \code{help(mt_as_move2)}.
#' @param traj2 (optional) same as traj, but for the second group of individuals. See \code{checkTO}
#' @param tc temporal tolerance limit (in seconds) for defining when two fixes
#'         are simultaneous or together. Parameter passed to function \code{GetSimultaneous}.
#' @param dc distance tolerance limit (in appropriate units) for defining when 
#'         two fixes are spatially together.
#'
#' @return
#'  This function returns a numeric result of the Ca statistic for each pair in the dataset.
#'
#' @references
#'  Bauman, P.J. (1998) The Wind Cave National Park elk herd: home ranges, seasonal movements, and alternative control methods.
#'    M.S. Thesis. South Dakota State University, Brookings, South Dakota, USA. \cr\cr
#'  Cole, L.C. (1949) The measurement of interspecific association. \emph{Ecology}. \bold{30}, 411--424.
#'
#' @keywords indices
#' @seealso GetSimultaneous, Prox, HAI
#' @examples
#' data(deer)
#' #tc = 7.5 minutes, dc = 50 meters
#' Ca(deer, tc = 7.5*60, dc = 50)
#' 
#' @export
#
# ---- End of roxygen documentation ----

Ca <- function(traj,traj2,tc=0,dc=0){
  
  #Time Units set to seconds
  units(tc) <- as_units("s")
  
  if (missing(traj2)){
    pairs <- checkTO(traj)
    pairs <- pairs[pairs$TO==TRUE,]
  } else {
    pairs <- checkTO(traj,traj2)
    pairs <- pairs[pairs$TO==TRUE,]
    traj <- rbind(traj,traj2)
  }
  
  n.pairs <- nrow(pairs)
  pairs$Ca <- NA
  
  for (i in 1:n.pairs){
    traj1 <- traj[mt_track_id(traj)==pairs$ID1[i],]
    traj2 <- traj[mt_track_id(traj)==pairs$ID2[i],]
    A <- nrow(traj1)
    B <- nrow(traj2)
    
    trajs <- GetSimultaneous(traj1,traj2,tc)
    
    traj1 <- trajs[mt_track_id(trajs)==pairs$ID1[i],]
    traj2 <- trajs[mt_track_id(trajs)==pairs$ID2[i],]
    
    trDist = st_distance(traj1,traj2,by_element=TRUE)
    
    #Unit control
    units(dc) <- units(trDist)

    AB <- length(which(trDist <= dc))
    
    #Compute coefficent of association
    pairs$Ca[i] <- 2*AB/(A+B)
  }
  
  return(pairs)
}
#=============== End of Ca Function =======================================
