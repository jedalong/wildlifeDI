# ---- roxygen documentation ----
#
#' @title Half-weight Association Index
#'
#' @description
#' This function computes the Half-weight Association Index for examining the presence of dynamic interaction in wildlife telemetry studies. This implementation follows that outlined in the paper Atwood and Weeks (2003).
#'
#' @details
#' This function can be used to test for the presence of dynamic interaction within the shared area (often termed the overlap zone) of the two animals home ranges. Specifically, HAI is calculated in identical fashion to that for \code{Ca}, but considers only those fixes in the shared area. Typically, the overlap zone (OZ) is easily obtained by taking the spatial intersection of two polygon home ranges.
#'
#' @param traj an object of the class \code{move2} which contains the time-stamped movement fixes of at least two individuals. For more information on objects of this type see \code{help(mt_as_move2)}.
#' @param traj2 (optional) same as traj, but for the second group of individuals. See \code{checkTO}
#' @param hr (optional)spatial polygon \code{sf} object associated with the home range (or some other form of) spatial range estimate for each individual in \code{traj}. The hr polygon should have a corresponding ID column with the same column name as in \code{traj}. If NULL (the default) the MCP home range estimate will be used for each individual.
#' @param tc time threshold for determining simultaneous fixes -- see function: \code{GetSimultaneous}.
#' @param dc distance tolerance limit (in appropriate units) for defining when two fixes are spatially together.
#'
#' @return
#' This function returns the numeric value of the HAI statistic. Values near 1 indicate attraction within the shared home range area, while values near 0 indicate avoidance within this shared area.
#'
#' @references
#' Atwood, T.C. and Weeks Jr., H.P. (2003) Spatial home-range overlap and temporal
#' interaction in eastern coyotes: The influence of pair types and fragmentation.
#' \emph{Canadian Journal of Zoology}, \bold{81}: 1589-1597.\cr\cr
#'
#' @keywords indices
#' @seealso GetSimultaneous, Ca
#' @examples
#' \dontrun{
#' data(deer)
#' 
#' #uses as a default minimum convex polygon for home range...
#' #tc = 7.5 minutes, dc = 50 meters
#' HAI(deer, tc=7.5*60, dc=50)
#' }
#' @export
#
# ---- End of roxygen documentation ----
HAI <- function(traj, traj2, hr=NULL, tc = 0,dc = 50){
  
  if (missing(traj2)){
    pairs <- checkTO(traj)
    pairs <- pairs[pairs$TO==TRUE,]
    mtraj <- traj
  } else {
    pairs <- checkTO(traj,traj2)
    pairs <- pairs[pairs$TO==TRUE,]
    if (st_crs(traj2) != st_crs(traj)){
      traj2 <- st_transform(traj2,crs=st_crs(traj))
    }
    mtraj <- data.frame(id = c(mt_track_id(traj),mt_track_id(traj2)),
                        time = c(mt_time(traj),mt_time(traj2)),
                        geometry = c(traj[[attr(traj,'sf_column')]],traj2[[attr(traj2,'sf_column')]])) |>
      st_as_sf(sf_column_name = "geometry", crs=st_crs(traj)) |>
      mt_as_move2(time_column='time',track_id_column='id')
  }
  
  n.pairs <- nrow(pairs)
  pairs$hai <- NA
  
  #get column name of ID column
  idcol <- mt_track_id_column(traj)
  
  #IF no hr is specified
  if (is.null(hr)){
    hr <- traj |>
      dplyr::group_by_at(idcol) |>
      dplyr::summarise() |>
      st_convex_hull()
  } else {
    #Check if the polygon has a column with the correct name
    if (any(idcol %in% names(hr))==FALSE){
      print(paste0("hr object does not have correctly named id column: ",idcol))
      return(NULL)
    }
  }
  
  for (i in 1:n.pairs){
    traj1 <- mtraj[mt_track_id(traj)==pairs$ID1[i],]
    traj2 <- mtraj[mt_track_id(traj)==pairs$ID2[i],]
    
    trajs <- GetSimultaneous(traj1,traj2,tc)
    
    traj1 <- trajs[mt_track_id(trajs)==pairs$ID1[i],]
    traj2 <- trajs[mt_track_id(trajs)==pairs$ID2[i],]
  
    n1 <- nrow(traj1)
    n2 <- nrow(traj2)
    
    hr1 <- hr[ hr[[idcol]] == pairs$ID1[i], ]
    hr2 <- hr[ hr[[idcol]] == pairs$ID2[i], ]
    
    oz <- suppressWarnings(st_intersection(hr1,hr2))
    
    #identify only those points in the OZ
    suppressWarnings(tr1x <- st_intersection(traj1,oz))
    suppressWarnings(tr2x <- st_intersection(traj2,oz))
  
    #Get only those fixes that are simultaneously in OZ
    tr1x <- tr1x[mt_time(tr1x) %in% mt_time(tr2x),]
    tr2x <- tr2x[mt_time(tr2x) %in% mt_time(tr1x),]
  
    #calculate the count of distances below threshold
    Do <- diag(st_distance(tr1x,tr2x))
    units(dc) <- units(Do)
    noz <- length(which(Do < dc))
    #Compute count of all other fixes that are not proximal
    x <- n1 - noz
    y <- n2 - noz
  
    pairs$hai[i] <- noz/(noz + 0.5*(x + y))
  }
  #print(paste("A critical distance of ",dc," was used, resulting in HAI = ", hai,".",sep=""))
  return(pairs)
}

#======================= End of HAI function ======================================
