# ---- roxygen documentation ----
#
#' @title Convert ltraj to sp object
#'
#' @description
#' The function \code{traj2sp} is a simple function to convert from a \code{ltraj} object to a \code{sp} object for plotting purposes. Note that this function is extremely basic.  
#' 
#' @details
#' If \code{type = 'p'} a \code{SpatialPointsDataFrame} is created from the \code{ltraj} object directly (see \code{?ld} in \code{adehabitatLT} for the columns that are generated). If \code{type = 'l'} a \code{SpatialLines} object is created, which allows simple plotting of the movement trajectory using a connect-the-dots method. At the moment, multiple individuals or 'bursts' are not considered.    
#'
#' @param traj an object of the class \code{ltraj}.
#' @param type 'p' for \code{SpatialPoints*}, 'l' for \code{SpatialLines}.
#' @param proj4string projection information in the form of a CRS object see \code{?CRS}.
#'
#' @return
#' A \code{SpatialPoints*} or \code{SpatialLines} object.  
#'
#' @keywords mapping
#' @seealso contacts
#' @examples
#' data(deer)
#' deer37 <- deer[1]
#' prj <- CRS("+proj=utm +zone=14 +ellps=WGS84 +datum=NAD83 +units=m +no_defs")
#' sp <- traj2sp(deer37,type='p',proj4string=prj)
#' sl <- traj2sp(deer37,type='l',proj4string=prj)
#' 
#' plot(sp,axes=T)
#' plot(sl,add=T)
#' 
#' @export
#

traj2sp <- function(traj,type='p',proj4string=NA){
  tr <- ld(traj)
  #how to deal with bursts?
  if (type == 'p'){
    sp <- SpatialPointsDataFrame(tr[,1:2],data=tr,proj4string=proj4string)
  } else if(type == 'l'){
    l1 <- Line(tr[,1:2])
    s1 <- Lines(list(l1),ID='a')
    sp <- SpatialLines(list(s1))
  }
  slot(sp,'proj4string',check=TRUE) <- proj4string
  return(sp)
}


