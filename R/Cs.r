# ---- roxygen documentation ----
#
#' @title Coefficient of Sociality
#'
#' @description
#'    The function \code{Cs} computes the coefficient of sociality between two moving
#'    objects following the methods outlined by Kenward et al. (1993). It also uses a
#'    signed Wilcoxon-rank test to test for significance.
#'
#' @details
#'  This function can be used to calculate the Kenward et al. (1993) coefficient of sociality (Cs)
#'  between two animals. The Cs statistic tests the observed mean
#'  distance between simultaneous fixes against that expected by the overall
#'  distribution of distances between all fixes.
#'  \deqn{Cs=\frac{D_E-D_O}{D_O+D_E}}{Cs=(De-Do)/(De+Do)}
#'  Where \eqn{D_O}{Do} is the mean observed distance between simultaneous fixes, and \eqn{D_E}{De}
#'  is the mean expected distance between all fixes. Kenward et al. (1993) propose Cs
#'  as a useful metric for exploring attraction or avoidance behaviour.
#'  Values for Cs closer to 1 indicate
#'  attraction, while values for Cs closer to -1 indicate avoidance. Values of Cs
#'  near 0 indicate that the two animals' movements have no influence on one another.
#'  \cr \cr
#'  Further, the difference between the observed and expected distances are compared
#'  using a paired signed-rank test (both one-sided tests, indicative of attraction
#'  or avoidance). See the function \code{GetSimultaneous} for details on how
#'  simultaneous fixes are determined from two trajectories.
#'
#' @param traj an object of the class \code{move2} which contains the time-stamped movement fixes of at least two individuals. For more information on objects of this type see \code{help(mt_as_move2)}.
#' @param traj2 (optional) same as traj, but for the second group of individuals. See \code{checkTO}
#' @param tc time threshold for determining simultaneous fixes -- see function: \code{GetSimultaneous}.
#'
#' @return
#'  This function returns a list of objects representing the calculated values from the
#'  Cs statistic and associated \emph{p}-values from the signed rank test.
#'  \itemize{
#'    \item Do -- The mean distance of simultaneous fixes.
#'    \item De -- The mean expected distance, from all fixes.
#'    \item Cs -- The coefficient of sociality, see \bold{Details}.
#'    \item p.Attract -- One sided \emph{p}-value from signed rank test, testing for attraction.
#'    \item p.Avoid -- One sided \emph{p}-value from signed rank test, testing for avoidance.
#'    }
#'
#' @references
#' Kenward, R.E., Marcstrom, V. and Karlbom, M. (1993) Post-nestling behaviour in
#'  goshawks, \emph{Accipiter gentilis: II}. Sex differences in sociality and nest-switching.
#'  \emph{Animal Behaviour}. \bold{46}, 371--378.
#'
#'
#' @keywords indices
#' @seealso GetSimultaneous
#' @examples
#' data(deer)
#' #tc = 7.5 minutes
#' Cs(deer, tc = 7.5*60) 
#' 
#' @export
#
# ---- End of roxygen documentation ----
Cs <- function(traj, traj2, tc=0){
  
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
  pairs$Do <- NA
  pairs$De <- NA
  pairs$Cs <- NA
  pairs$p.attract <- NA
  pairs$p.avoid <- NA
  
  for (i in 1:n.pairs){
    traj1 <- mtraj[mt_track_id(traj)==pairs$ID1[i],]
    traj2 <- mtraj[mt_track_id(traj)==pairs$ID2[i],]
  
    trajs <- GetSimultaneous(traj1, traj2, tc)
  
    traj1 <- trajs[mt_track_id(trajs)==pairs$ID1[i],]
    traj2 <- trajs[mt_track_id(trajs)==pairs$ID2[i],]
    n <- nrow(traj1)

    #calculate the observed distances
    dM = st_distance(traj1,traj2)
    
    Do <- as.numeric(diag(dM))

    #calculate the expected distances
    De = rowSums(dM)/n
  
    #calculate the significance of differences b/w Do and De using a Wilcoxon signed rank test
    p.less <- wilcox.test(Do,De,paired=T,alternative="less")$p.value
    p.great <- wilcox.test(Do,De,paired=T,alternative="greater")$p.value

    #Compute the Coefficient of sociality (Kenward et al. 1993)
    DO <- sum(Do)/n
    DE <- sum(De)/n
    Cs <- (DE - DO) / (DO + DE)

    #return output
    pairs$Do[i] <- DO
    pairs$De[i] <- DE
    pairs$Cs[i] <- Cs
    pairs$p.attract[i] <- p.less
    pairs$p.avoid[i] <- p.great
    
  }
  return(pairs)
}


