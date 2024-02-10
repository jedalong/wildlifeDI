# ---- roxygen documentation ----
#
#' @title Get period where two tracks overlap
#'
#' @description
#'   The function \code{GetTO} identifies and extracts fixes of a tracking dataset that overlap in time with all other trajectories.
#'
#' @details
#'   This function is used to determine the fixes that overlap in time between two trajectories.
#'
#' @param traj an object of the class \code{move2} which contains the time-stamped movement fixes of >= 1 individual. For more information on objects of this type see \code{help(mt_as_move2)}.
#' @param tb (optional) time threshold (i.e., time buffer) for considering if fixes are in the overlap period (in seconds).
#'
#' @return
#' A single \code{move2} object containing the fixes from traj that temporally overlap. If more than 2 individuals it will take the maximum of the earliest start-time from all individuals and minimum of the latest end-time of all individuals.
#'
# @references
#'
#' @keywords processing
#' @seealso checkTO
#' @examples
#' data(deer)
#' deer_to <- GetTO(deer)
#' 
#' @export
#
# ---- End of roxygen documentation ----

GetTO <- function(traj,tb=0){
  
  #global variables in group_by hack
  id <- NULL
  time <- NULL
  
  #units(tb) <- 's'
  #Get nearest fixes that are within tc from one another.
  df <- data.frame(id=mt_track_id(traj),time=mt_time(traj))
  
  #identify the temporal overlap period
  tdf <- df |>
    dplyr::group_by(id) |>
    dplyr::summarise(start=min(time),
                     end=max(time))
  
  t.min <- max(tdf$start) - tb
  t.max <- min(tdf$end) + tb
  
  #check to see if there is an overlap period
  if (t.min > t.max){stop('There is not a defined period of temporal overlap.')}
  
  #find the fixes that are in the temporal overlap
  traj <- traj[mt_time(traj) >= t.min & mt_time(traj) <= t.max, ]

  #Return the move2 object in overlap period
  return(traj)
}

