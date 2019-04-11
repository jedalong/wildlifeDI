# ---- roxygen documentation ----
#
#' @title Filter trajectory based on conditions
#'
#' @description
#' The function \code{FilterTraj} is a function for extracting portions of a trajectory based on some filter criteria. 
#' 
#' @details
#' The function \code{FilterTraj} can be used to extract a portion of a trajectory based on a spatial area (i.e., a \code{SpatialPolygons} object) or based on temporal criteria (such as time-of-day or a temporal window) or based on attributes of the trajectory. In the case of attribute filtering, the criteria can be any of the default attributes of an \code{ltraj} object, or based on additional attributes stored alongside the trajectory. The \code{filter} parameter is a flexible way to define the criteria for filtering. The \code{filter} object can be one of: \cr
#' - \code{SpatialPolygons*} objects for spatial filtering;\cr
#' - \code{POSIX} class objects (list of length=2) for filtering based on some temporal window (see examples);\cr
#' - \code{character} object (list of length=2) for filtering based on time-of-day (see examples);\cr
#' - \code{character} object containing the attribute filter criteria as a logical expression (see examples and \code{?subset}).\cr
#' \cr
#'NOTE: When using \code{FilterTraj} be very careful about using the output \code{ltraj} object in further analysis. When a new \code{ltraj} object is created, all the movement parameters (e.g., dist, dt) are recalculated and thus not necessarily valid. Therefore, subsequent analysis should ideally focus only on the raw fix information (i.e., x, y, date).
#'
#' @param traj an object of the class \code{ltraj} which contains the time-stamped movement fixes of the object. Note this object must be a \code{type II ltraj} object. For more information on objects of this type see \code{help(ltraj)}.
#' @param type The type of filter to apply (one of 'attribute','spatial','temporal','tod').
#' @param filter The filter criteria (see details).
#'
#' @return
#' A \code{ltraj} object with only those fixes satisfying the \code{filter} criteria.  
#'
#' @keywords processing
#' @seealso GetTemporalOverlap
#' @examples
#' data(deer)
#' deer37 <- deer[1]
#' deer38 <- deer[2]
#' 
#' #Spatial Filter 
#' #----NOT RUN----
#' \dontrun{
#' library(adhabitatHR)
#' oz <- mcp(SpatialPoints(cbind(ld(deer38)$x,ld(deer38)$y)))
#' x <- FilterTraj(deer37,type='spatial',filter=oz)
#' }
#' #
#' #---------------
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

FilterTraj <- function(traj, type='attribute',filter=NA){
  tr1 <- ld(traj)
  if (type == 'spatial'){
    #Check territory.sp
    if (class(filter) %in% c('SpatialPolygonsDataFrame','SpatialPolygons')){check <- TRUE} else {stop('Parameter Error: filter is not a SpatialPolygons* object.')}
    tr1.sp <- SpatialPointsDataFrame(SpatialPoints(tr1[,1:2]),data=tr1)
    tr1.int <- tr1.sp[gWithin(tr1.sp,filter,byid=TRUE)[1,],]  # Better way to do this?
    tr1 <- tr1.int@data
  } else if (type == 'temporal'){
    if (class(filter)[1] %in% c('POSIXct','POSIXlt')){check <- TRUE} else {stop('Parameter Error: filter is not a POSIX object.')}
    if (length(filter) != 2){stop('Parameter Error: filter is not a POSIX list length=2.')}
    ind <- which(tr1$date > filter[1] & tr1$date <= filter[2])
    tr1 <- tr1[ind,]
  } else if(type == 'tod'){
    if (class(filter) != 'character'){stop('Parameter Error: filter is not a character object.')}
    if (length(filter) != 2){stop('Parameter Error: filter is not a character list length=2.')}
    tr1.time <- strftime(tr1$date, format = '%H:%M:%S')
    ind <- which(tr1.time > filter[1] & tr1.time <= filter[2])
    tr1 <- tr1[ind,]
  } else if (type == 'attribute'){
    if (class(filter) != 'character'){stop('Parameter Error: filter is not a character object.')}
    tr1 <- subset(tr1,eval(parse(text=filter)))
  }
  #convert to ltraj object
  if (dim(tr1)[1] == 0){
    return(NA)
  } else {
    out.traj <- dl(tr1)
    #Return the ltraj object
    return(out.traj)
  }
}


