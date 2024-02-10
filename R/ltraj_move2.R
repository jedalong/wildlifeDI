# ---- roxygen documentation ----
#
#' @title Convert ltraj to move2 object
#'
#' @description
#' The function \code{ltraj_move2} is a simple function for quickly converting \code{ltraj} to \code{move2} objects. 
#' 
#' @details
#' The function \code{ltraj_move2} can be used to convert an \code{ltraj} object into an \code{move2} object.
#' 
#' @param ltraj an object of the class \code{ltraj} which contains the time-stamped movement fixes of the object. For more information on objects of this type see \code{help(ltraj)}.
#'
#' @return
#' A \code{move2} object.  
#'
#' @keywords processing
#' 
#' @seealso move2_ltraj
#' @examples
#' data(deer)
#' deer_ltraj <- move2_ltraj(deer)
#' deer_move <- ltraj_move2(deer_ltraj)
#' 
#' @export
#

ltraj_move2 <- function(ltraj){
  prj4string <- attr(ltraj,'proj4string')
  if (is.null(prj4string)){prj4string = NA}
  m2 <- mt_as_move2(adehabitatLT::ld(ltraj), coords=c('x','y'),
                    time_column='date',
                    track_id_column='id',
                    na.fail=FALSE) |> st_set_crs(prj4string)
  return(m2)
}