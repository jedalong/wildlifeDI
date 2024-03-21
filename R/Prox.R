# ---- roxygen documentation ----
#
#' @title Proximity Index
#'
#' @description
#' The function \code{Prox} simply computes the proportion of (simultaneous) fixes that are proximal, based on some spatial
#' threshold -- \code{dc} (Bertrand et al. 1996). It also facilitates local-level proximity analysis
#'
#' @details
#' The function \code{Prox} can be used to test for the presence of attraction (via proximity) in wildlife telemetry data. Prox is simply the proportion of simultaneous fixes within the threshold distance -- \code{dc}. The local output (dataframe) can be useful for examining variation in proximity through time.
#'
#' @param traj an object of the class \code{move2} which contains the time-stamped movement fixes of at least two individuals. For more information on objects of this type see \code{help(mt_as_move2)}.
#' @param traj2 (optional) same as traj, but for the second group of individuals. See \code{checkTO}
#' @param tc time threshold for determining simultaneous fixes -- see function: \code{GetSimultaneous}.
#' @param dc distance tolerance limit (in appropriate units) for defining when two fixes are spatially together.
#' @param local logical value indicating. When local = FALSE (the default) prox returns a data.frame with the global proximity ratio (proportion of all fixes below dc and tc) for each pair of individuals. When local = TRUE, prox returns the input move2 object with the distance to the most proximal fix, and the number of fixes that are considered proximal for each fix in the dataset.
#'
#' @return
#' If \code{local=FALSE} (the default) Prox returns the numeric value of the Prox index for each pair of individuals.
#' If \code{local=TRUE} Prox returns a \code{move2} containing the original trajectory (or both trajectories) with three additional columns prox (the distance to the nearest proximal fix), prox.id (the id of the nearest proximal fix) and prox.n (the number of individuals with proximal fixes)
#'
#' @references
#' Bertrand, M.R., DeNicola, A.J., Beissinger, S.R, Swihart, R.K. (1996) Effects of parturition
#' on home ranges and social affiliations of female white-tailed deer.
#' \emph{Journal of Wildlife Management}, \bold{60}: 899-909.
#'
#' @keywords indices
#' @seealso GetSimultaneous, contacts
#' @examples
#' data(deer)
#' #tc = 7.5 minutes, dc = 50 meters
#' Prox(deer, tc=7.5*60, dc=50)
#' deer <- Prox(deer, tc=7.5*60, dc=50, local=TRUE)
#' 
#' @export
#
# ---- End of roxygen documentation ----

#================= New Prox function ============================================
Prox <- function(traj,traj2,tc=0,dc=50,local=FALSE){
  
  #Combine trajectories and identify overlap pairs
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
    mtraj <- mt_stack(traj,traj2,track_combine='check_unique')
  }
  
  n.pairs <- nrow(pairs)
  pairs$prox <- NA
  
  if (local) {
    local.df <- NULL
    traj$prox <- NA
    traj$prox.id <- NA
    traj$prox.n <- 0
  }
  
  for (i in 1:n.pairs){
    traj1 <- mtraj[mt_track_id(traj)==pairs$ID1[i],]
    traj2 <- mtraj[mt_track_id(traj)==pairs$ID2[i],]
    
    trajs <- GetSimultaneous(traj1,traj2,tc)
    
    traj1 <- trajs[mt_track_id(trajs)==pairs$ID1[i],]
    traj2 <- trajs[mt_track_id(trajs)==pairs$ID2[i],]
    n <- nrow(traj1)
    
    prox.dist <- as.numeric(st_distance(traj1,traj2,by_element=TRUE))
    
    #compute the Prox index
    nprox <- length(which(prox.dist < dc))
    pairs$prox[i] <- nprox/n
    
    if (local){
      #add proximity values for each individual
      df1 <- data.frame(id1=pairs$ID1[i],id2=pairs$ID2[i],prox=prox.dist,row.name = row.names(traj1))
      df2 <- data.frame(id1=pairs$ID2[i],id2=pairs$ID1[i],prox=prox.dist,row.name = row.names(traj2))
      df3 <- rbind(df1,df2)
      df3 <- df3[df3$prox < dc,]
      
      if (is.null(local.df)){
        local.df <- df3
      } else {
        local.df <- rbind(local.df,df3)
      }
    }
  }
  
  #organize outputs
  if (local){
    for (j in 1:nrow(traj)){
      pdf <- local.df[local.df$row.name == row.names(traj)[j],]
      if (nrow(pdf) > 0){
        traj$prox[j] <- min(pdf$prox)
        traj$prox.id[j] <- pdf$id2[which.min(pdf$prox)]
        traj$prox.n[j] <- nrow(pdf)
      } 
    }
    return(traj)
  } else {
    return(pairs)
  }
}
  
#================= end of Prox function ========================================