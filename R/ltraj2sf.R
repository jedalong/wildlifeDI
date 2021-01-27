# ---- roxygen documentation ----
#
#' @title Convert ltraj to sf spatial object
#'
#' @description
#' The function \code{ltraj2sf} is a simplee function for converting ltraj to sf objects. 
#' 
#' @details
#' The function \code{ltraj2sf} can be used to convert an \code{ltraj} object into an \{sf} spatial ojbect (either as points or lines).
#' 
#' @param traj an object of the class \code{ltraj} which contains the time-stamped movement fixes of the object. For more information on objects of this type see \code{help(ltraj)}.
#' @param type One of \code{"point"} (the default) or \code{"line"}.
#'
#' @return
#' A \code{sf} object either points or lines.  
#'
#' @keywords processing
#' @seealso conSpatial
#' @examples
#' data(deer)
#' deer37 <- deer[1]
#' 
#' #Temporal Filter
#' t1 <- as.POSIXct(strptime('2005-03-09 00:00:00', format= '%Y-%m-%d %H:%M:%S'))
#' t2 <- as.POSIXct(strptime('2005-03-11 00:00:00', format= '%Y-%m-%d %H:%M:%S'))
#' twin <- c(t1,t2)
#' x <- FilterTraj(deer37,type='temporal',filter=twin)
#' 
#' #tod Filter
#' tod <- c('06:00:00', '10:00:00')
#' x <- FilterTraj(deer37,type='tod',filter=tod)
#' 
#' #attribute Filter
#' q <- 'dist > 100'
#' x <- FilterTraj(deer37,type='attribute',filter=q)
#'
#' q <- 'dist > 100 & rel.angle < 1'
#' x <- FilterTraj(deer37,type='attribute',filter=q)
#' 
#' @export
#

ltraj2sf <- function(traj, type='point'){
  tr1 <- ld(traj)
  
  prj4string <- attr(traj,'proj4string')
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


