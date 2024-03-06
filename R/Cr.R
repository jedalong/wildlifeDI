# ---- roxygen documentation ----
#
#' @title Movement Correlation Coefficient
#'
#' @description
#'   The function \code{Cr} computes the correlation statistic for movement data as presented 
#'   in the paper by Shirabe (2006). The statistic is essentially a Pearson product-moment 
#'   correlation statistic formulated for use with movement data.
#'   
#' @details
#'   The function \code{Cr} can be used to measure the level of dynamic interaction (termed correlation)
#'   between a pair of simultaneously moving objects. The statistic is sensitive to
#'   interaction in both movement direction (azimuth) and displacement, but is unable to
#'   disentangle the effects of these components. 
#'   NOTE: This function is only appropriate with projected coordinates.
#'
#' @param traj an object of the class \code{move2} which contains the time-stamped movement fixes of at least two individuals. For more information on objects of this type see \code{help(mt_as_move2)}.
#' @param traj2 (optional) same as traj, but for the second group of individuals. See \code{checkTO}
#' @param tc time threshold for determining simultaneous fixes -- see function: \code{GetSimultaneous}.
#'
#' @return
#'   This function returns the Shirabe (2006) correlation statistic for two moving objects.
#'
#' @references
#' Shirabe, T. 2006. Correlation analysis of discrete motions. In: Raubal, M.,
#' Miller, HJ, Frank, AU, and Goodchild, M. eds. GIScience 2006, LNCS 4197. Berlin: Springer-Verlag;
#' 370-382.
#'
#' @keywords indices
#' @seealso GetSimultaneous, DI
#' @examples
#' data(deer)
#' #tc = 7.5 minutes
#' Cr(deer, tc = 7.5*60)
#' 
#' @export
#
# ---- End of roxygen documentation ----
Cr <- function(traj, traj2, tc = 0){
  
  
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
  pairs$Cr <- NA
  
  for (i in 1:n.pairs){
    traj1 <- mtraj[mt_track_id(traj)==pairs$ID1[i],]
    traj2 <- mtraj[mt_track_id(traj)==pairs$ID2[i],]
    
    trajs <- GetSimultaneous(traj1, traj2, tc)
    
    traj1 <- trajs[mt_track_id(trajs)==pairs$ID1[i],]
    traj2 <- trajs[mt_track_id(trajs)==pairs$ID2[i],]
    
    tr1 <- st_coordinates(traj1)
    tr2 <- st_coordinates(traj2)
    n <- nrow(tr1)
    
    #compute vectors for each movement segment
    for (j in n:2)
    {
      tr1[j,] <- tr1[j,] - tr1[j-1,]
      tr2[j,] <- tr2[j,] - tr2[j-1,]
    }
    #Remove the beginning fix
    tr1 <- tr1[2:n,]
    tr2 <- tr2[2:n,]
    
    #compute path means for x,y coords
    tr1.bar <- apply(tr1,2,mean)
    tr2.bar <- apply(tr2,2,mean)
    #compute mean subtracted data matrix
    tr1.m <- t(t(as.matrix(tr1)) - tr1.bar)
    tr2.m <- t(t(as.matrix(tr2)) - tr2.bar)
    #function for numerator
    numer <- function(a,b)
    {
      numer <- sum(a*b)
      return(numer)
    }
    #compute numerator
    R.numer <- rep(0,(n-1))
    for (k in 1:(n-1))
    {
      R.numer[k] <- numer(tr1.m[k,],tr2.m[k,])
    }
    #compute denominator
    len <- function(v)
    {
      l <- sqrt(sum(v*v))
      return(l)
    }
    tr1.denom <- sqrt(sum(apply(tr1.m,1,len)^2))
    tr2.denom <- sqrt(sum(apply(tr2.m,1,len)^2))
    R.denom <- tr1.denom*tr2.denom
    #compute and return the statistic
    pairs$Cr[i] <- sum(R.numer) / R.denom
  }
  return(pairs)
}


#======================== End of shirabe Function ==============================