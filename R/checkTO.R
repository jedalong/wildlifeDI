# ---- roxygen documentation ----
#
#' @title Check for temporal overlap
#'
#' @description
#' The function \code{checkTO} is a simple function for identifying if, and for how long, two telemetry datasets overlap (temporally) with each other. The function returns a dataframe with5 columns of information: the ids of the first an second individuals in a dyad, a logical variable indicating if the two trajectories overlap temporally, and timings of the beginning and end of the overlap period. If only a single move2 object is provided it considers all pairwise dyads within that move2 object. If two move2 objects are passed in it considers only the dyad pairs from traj against traj2. This can be used to test only the pairwise dyads between two groups (e.g., inter-species).
#' 
#' @details
#' The function \code{checkTO} can be used to identify if, when, and for how long the tracking data of two individuals overlap temporally.   
#'
#' @param traj an object of the class \code{move2} which contains the time-stamped movement fixes. For more information on objects of this type see \code{help(mt_as_move2)}.
#' @param traj2 (optional) an object of the class \code{move2} which contains the time-stamped movement fixes. For more information on objects of this type see \code{help(mt_as_move2)}.
#'
#' @return
#' A \code{data.frame} of with five columns, ID1, ID2, TO (logical indicating if the two tracking dataset overlap temporally), the beginning (\code{tmin}), and end (\code{tmax}) of the overlap period, stored as \code{POSIX} objects.  
#'
#' @keywords processing
#' @seealso GetSimultaneous, GetTO
#' @examples
#' data(does)
#' dyads <- checkTO(does)
#' 
#' @export
#

checkTO <- function(traj,traj2){
  
  
  if (missing(traj2)){
    id.list <- unique(mt_track_id(traj))
    if (length(id.list) <= 1) {
      print("There was only 1 unique individual in the input dataset. Please try again with a dataset containing at least 2 individuals.")
      return(NULL)
    }
    pairs <- expand.grid(id.list,id.list,stringsAsFactors=F)
  } else {
    id.list1 <- unique(mt_track_id(traj))
    id.list2 <- unique(mt_track_id(traj2))
    pairs <- expand.grid(id.list1,id.list2,stringsAsFactors=F)
    traj <- mt_stack(traj,traj2)
  }
  
  #Get all the unique combinations between one group
  names(pairs) <- c('ID1','ID2')
  pairs <- pairs[which(pairs$ID1 != pairs$ID2),]
  pairs <- pairs[order(pairs$ID1),]
  pairs <- pairs[!duplicated(t(apply(pairs, 1, sort))),]   #Do we always want to get rid of duplicated pairs, i think so...
  n.pairs <- nrow(pairs)
  
  if (n.pairs == 0){
    print('No pairs with temporal overlap found.')
    return(NULL)
  }
  
  pairs$TO <- FALSE
  pairs$t.min <- NA
  pairs$t.max <- NA
  dfr <- NULL
  
  for (i in 1:n.pairs){
    traj1 <- traj[mt_track_id(traj)==pairs$ID1[i],]
    traj2 <- traj[mt_track_id(traj)==pairs$ID2[i],]
    
    #Get nearest fixes that are within tc from one another.
    t1 <- mt_time(traj1)
    t2 <- mt_time(traj2)
    
    #identify the temporal overlap period
    t.min <- max(c(min(t1),min(t2)))
    t.max <- min(c(max(t1),max(t2)))
    
    #check to see if there is an overlap period
    if (t.min > t.max){
      TO <- FALSE
      t.min <- t.max <- NA
    } else {
      TO <- TRUE
    }
    pairs$TO[i] <- TO
    pairs$t.min[i] <- t.min
    pairs$t.max[i] <- t.max
    
  }
  
  #origin of numerical dates
  origin = "1970-01-01"
  
  if (is.numeric.POSIXt(pairs$t.min)) {
    #do nothing
  } else {
    pairs$t.min <- as.POSIXct(pairs$t.min, origin = origin)
  }
  if (is.numeric.POSIXt(pairs$t.max)) {
    #do nothing
  } else {
    pairs$t.max <- as.POSIXct(pairs$t.max, origin = origin)
  }
  
  return(pairs)
}
