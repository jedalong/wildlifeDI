# ---- roxygen documentation ----
#
#' @title Convert sf spatial object to ltraj
#'
#' @description
#' The function \code{sf2ltraj} is a simple function for converting sf objects to ltraj objects. 
#' 
#' @details
#' The function \code{sf2ltraj} can be used to convert an \code{sf} object (stored as points) into an \code{ltraj} object.
#' 
#' @param sfp an object of the class \code{sf} which contains the time-stamped movement fixes (as points) of the object(s).
#' @param date the name of the column containing the date-time information of each fix. Must be converted to a POSIXct type.
#' @param id  either a character string indicating the identity of the animal or the name of the column with the ID of the individuals.
#'
#' @return
#' A \code{ltraj} object. For more information on objects of this type see \code{help(ltraj)} and \code{?as.ltraj}.
#'
#' @keywords processing
#' @seealso conSpatial, ltraj2sf
#' @examples
#' data(deer)
#' #points
#' deer_pt <- ltraj2sf(deer)
#' deer_tr <- sf2ltraj(deer_pt, date='date', id='id')
#'
#' @export
#

sf2ltraj <- function(sfp, date, id){
  
  #get coords
  xy = data.frame(st_coordinates(sfp))
  names(xy) = c('x','y')
  
  #get date
  df = st_drop_geometry(sfp)
  datetime = df[,date]
  
  #get id
  if (id %in% names(df)){
    id = df[,id]
  } else {
    id = id
  }
  
  
  #create ltraj object
  crs <- st_crs(sfp)
  if (is.na(crs)){
    trj = as.ltraj(xy,date=datetime,id=id)
  } else {
    prj4string <- crs$proj4string
    trj = as.ltraj(xy,date=datetime,id=id,proj4string=CRS(prj4string))
  }
  
  
  return(trj)
}