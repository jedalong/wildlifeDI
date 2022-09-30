# ---- roxygen documentation ----
#
#' @title Convert ltraj to sf spatial object
#'
#' @description
#' The function \code{ltraj2sf} is a simple function for converting ltraj to sf objects. 
#' 
#' @details
#' The function \code{ltraj2sf} can be used to convert an \code{ltraj} object into an \code{sf} spatial object (either as points or lines).
#' 
#' @param traj an object of the class \code{ltraj} which contains the time-stamped movement fixes of the object. For more information on objects of this type see \code{help(ltraj)}.
#' @param type One of \code{"point"} (the default) or \code{"line"}.
#'
#' @return
#' A \code{sf} object either points or lines.  
#'
#' @keywords processing
#' @seealso conSpatial,sf2ltraj
#' @examples
#' data(deer)
#' #points
#' deer_pt <- ltraj2sf(deer)
#' plot(deer_pt['id'])
#' 
#' #lines
#' deer_ln <- ltraj2sf(deer,type='line')
#' plot(deer_ln['id'])
#' @export
#

ltraj2sf <- function(traj, type='point'){
  tr1 <- ld(traj)
  
  prj4string <- attr(traj,'proj4string')
  if (is.null(prj4string)){prj4string = NA}
  spo = st_as_sf(tr1,coords=c('x','y'),crs=prj4string)
  
  id = NULL #fix global variable issue in package
  if (type=='line'){
    g = dplyr::group_by(spo,id)
    s = dplyr::summarise(g,
                         traj_begin = min(date),
                         traj_end = max(date),
                         npt = dplyr::n(),
                         do_union=FALSE)
    s2 = s[s$npt > 1,]
    spl = st_cast(s2,"LINESTRING")
    return(spl)
  } else {
    return(spo)
  }
}


