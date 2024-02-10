# ---- roxygen documentation ----
#
#' @title Convert move2 to ltraj object
#'
#' @description
#' The function \code{move2_ltraj} is a simple function for quickly converting \code{move2} to \code{ltraj} objects. 
#' 
#' @details
#' The function \code{move2_ltraj} can be used to convert a \code{move2} object into a \code{ltraj} object.
#' 
#' @param traj an object of the class \code{move2} which contains the time-stamped movement fixes of the object. For more information on objects of this type see \code{help(move2)}.
#'
#' @return
#' A \code{ltraj} object.  
#'
#' @keywords processing
#' 
#' @seealso ltraj_move2
#' @examples
#' data(deer)
#' deer_ltraj <- move2_ltraj(deer)
#' 
#' @export
#

move2_ltraj <- function(traj){
  xy <- st_coordinates(traj)
  date <- mt_time(traj)
  id <- mt_track_id(traj)
  infolocs <- traj |>
    st_drop_geometry() |>
    data.frame()
 
  if (is.na(st_crs(traj))){final_crs = NULL} else {final_crs <- sp::CRS(st_crs(traj)$proj4string)}

  l2 <- adehabitatLT::as.ltraj(xy,date,id,infolocs=infolocs,proj4string = final_crs)
  return(l2)
}