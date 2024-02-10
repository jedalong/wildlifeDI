# ---- roxygen documentation ----
#' @title Process contacts
#'
#' @description
#' This function performs basic contact analysis between individuals in a group of tracked animals, or between two different groups of tracked animals.
#'
#' @details
#' This function can be used to identify all fixes defined as contacts in space and time between individuals in one or two groups.

#' @param traj an object of the class \code{move2} which contains the time-stamped movement fixes of at least two individuals. For more information on objects of this type see \code{help(mt_as_move2)}. If traj2 is specified traj may have only one individual.
#' @param traj2 (optional) same as traj, but for the second group of individuals. See \code{checkTO}
#' @param tc time threshold for determining simultaneous fixes -- see function: GetSimultaneous.
#' @param dc distance tolerance limit (in appropriate units) for defining when two fixes are spatially together.
#' @param GetSim (logical) whether or not to use GetSimultaneous to time match fixes between pairs of individuals. Default = TRUE.
#' @param return What to return (one of 'move2' (default) or 'contacts'). See Return below.
#' 
#' @return If return = 'move2' (the default) this function returns the input traj move2 object with additional columns: contact - (binary) whtether or not a fix is a contact, contact_id - the id of the individual with which a contact occurs, contact_d - the proximity distance of the contact, contact_dt - the difference in time between the two fixes in the contact, contact_n - the number of contacts at that time. In the event that there is more than one contact for a given fix, the contact_id, contact_d, and contact_dt values are all associated with the most proximal (in geographical space) contact. If return = 'contacts' this function returns a data.frame with the columns: (id1,id2) the id's of the individuals involved in a contact, (row1,row2) the rownames from the original data associated with each of the fixes involved in a contact, (dist) the distance between the two fixes associated with the contact, and (difftime) the difference in time between the two fixes involved in the contact. 
#'
#' @references
#'  Long, JA, Webb, SL, Harju, SM, Gee, KL (2022) Analyzing Contacts and Behavior from High Frequency 
#'  Tracking Data Using the wildlifeDI R Package. \emph{Geographical Analysis}. \bold{54}, 648--663.
#' 
#'
#' @keywords contacts
#' @seealso GetSimultaneous, dcPlot, conPhase
#' 
#' @examples 
#' \dontrun{
#' data(does)
#' doecons <- conProcess(does,tc=15*60,dc=50)
#' }
#' 
#' @export
#
# ---- End of roxygen documentation ----

conProcess <- function(traj,traj2,dc=0,tc=0,GetSim=TRUE,return='move2'){
  
  #global variables in group_by hack
  id1 <- NULL
  dist <- NULL
  row1 <- NULL
  
  #Unit control
  units(tc) <- as_units("s")
  

  if (missing(traj2)){
    pairs <- checkTO(traj)
    pairs <- pairs[pairs$TO==TRUE,]
    mtraj <- traj
  } else {
    pairs <- checkTO(traj,traj2)
    pairs <- pairs[pairs$TO==TRUE,]
    mtraj <- rbind(traj,traj2)
  }
  
  n.pairs <- nrow(pairs)
  condf <- NULL
  for (i in 1:n.pairs){
    traja <- mtraj[mt_track_id(mtraj)==pairs$ID1[i],]
    trajb <- mtraj[mt_track_id(mtraj)==pairs$ID2[i],]
    
    if (GetSim){
      trajs <- GetSimultaneous(traja,trajb,tc)
      tr1 <- trajs[mt_track_id(trajs)==pairs$ID1[i],]
      tr2 <- trajs[mt_track_id(trajs)==pairs$ID2[i],]
      proxdf <- data.frame(id1=pairs$ID1[i],id2=pairs$ID2[i],
                           row1=row.names(tr1),row2=row.names(tr2),
                           dist=st_distance(tr1,tr2,by_element=TRUE),
                           difftime=as.numeric(abs(mt_time(tr1)-mt_time(tr2))))
      #Unit control
      units(proxdf$difftime) <- as_units('s')
      units(dc) <- units(proxdf$dist)
      proxdf <- proxdf[ proxdf$dist < dc & proxdf$difftime < tc, ]
      
    } else { #If we don' use GetSim We can have multiple contacts from different sources.
      dM <- st_distance(traja,trajb)
      tM <- abs(-outer(as.numeric(mt_time(traja)),as.numeric(mt_time(trajb)),'-'))
      
      units(tM) <- as_units('s')
      units(dc) <- units(dM)
      
      rownames(tM) <- rownames(traja)
      colnames(tM) <- rownames(trajb)
      
      ind <- which(dM < dc & tM < tc, arr.ind=TRUE)
      
      if (length(ind) > 0){
        rnm1 <- rownames(tM)[which(dM < dc & tM < tc, arr.ind = TRUE)[, 1]]
        rnm2 <- colnames(tM)[which(dM < dc & tM < tc, arr.ind = TRUE)[, 2]]
        proxdf <- data.frame(id1=pairs$ID1[i],id2=pairs$ID2[i],row1=rnm1,row2=rnm2,dist=dM[ind],difftime=tM[ind])
      } else {
        proxdf <- NULL
      }
      
    }
    if (is.null(condf)) {condf <- proxdf} else {condf <- rbind(condf,proxdf)}
  }
  
  
  
  #arrange return object
  if (return=='contact'){
    return(condf)
  } else {
    #Create contact list for both pair directions in case of 1 group
    if (missing(traj2)){
      dfr <- rbind(condf, 
                   data.frame(id1=condf$id2,id2=condf$id1,row1=condf$row2,row2=condf$row1,dist=condf$dist,difftime=condf$difftime))
    } else {
      dfr <- condf
    }
    
    dfr.d <- dfr |>
      dplyr::group_by(id1,row1) |>
      dplyr::slice_min(dist)
    
    dfr.n <- dfr |>
      dplyr::group_by(id1,row1) |>
      dplyr::summarise(ncon=dplyr::n(),
                       .groups='drop')
    
    traj$contact <- 0
    traj$contact_id <- NA
    traj$contact_d <- NA
    traj$contact_dt <- NA
    traj$contact_n <- NA
    
    ind.con <- match(dfr.d$row1, row.names(traj))
    
    traj$contact[ind.con] <- 1
    traj$contact_id[ind.con] <- dfr.d$id2
    traj$contact_d[ind.con] <- dfr.d$dist
    traj$contact_dt[ind.con] <- dfr.d$difftime
    traj$contact_n[ind.con] <- dfr.n$ncon
    
    return(traj)
  }

}























# 
# 
# 
# ### TO DELETE
# #Process contact analysis for all pairs of individuals in one group
# oneGroup <- function(traj,dc,tc,GetSim){
#   #Get Pairs
#   pairs <- checkTO(traj)
#   n.pairs <- nrow(pairs)
#   condf <- NULL
#   for (i in 1:n.pairs){
#     traj1 <- traj[mt_track_id(mtraj1)==pairs$ID1[i],]
#     traj2 <- traj[mt_track_id(mtraj1)==pairs$ID2[i],]
#     
#     if (GetSim){
#       trajs <- GetSimultaneous(traj1,traj2,tc)
#       tr1 <- trajs[mt_track_id(trajs)==pairs$ID1[i],]
#       tr2 <- trajs[mt_track_id(trajs)==pairs$ID2[i],]
#       proxdf <- data.frame(id1=pairs$ID1[i],id2=pairs$ID2[i],row1=row.names(tr1),row2=row.names(tr2),dist=st_distance(tr1,tr2,by_element=TRUE),difftime=as.numeric(abs(mt_time(tr1)-mt_time(tr2))))
#       #Unit control
#       units(proxdf$difftime) <- as_units('s')
#       units(dc) <- units(proxdf$dist)
#       proxdf <- subset(proxdf,dist < dc & difftime < tc)
#       
#       ##### NEED TO FIX WHEN GET SIM IS FALSE
#     } else {
#       dM <- st_distance(traj1,traj2)
#       tM <- abs(-outer(as.numeric(mt_time(traj1)),as.numeric(mt_time(traj2)),'-'))
#       
#       units(tM) <- as_units('s')
#       units(dc) <- units(dM)
#       
#       rownames(tM) <- rownames(traj1)
#       colnames(tM) <- rownames(traj2)
#       
#       ind <- which(dM < dc & tM < tc, arr.ind=TRUE)
#       
#       if (length(ind) > 0){
#         rnm1 <- rownames(tM)[which(dM < dc & tM < tc, arr.ind = TRUE)[, 1]]
#         rnm2 <- colnames(tM)[which(dM < dc & tM < tc, arr.ind = TRUE)[, 2]]
#         proxdf <- data.frame(id1=pairs$ID1[i],id2=pairs$ID2[i],row1=rnm1,row2=rnm2,dist=dM[ind],difftime=tM[ind])
#       } else {
#         proxdf <- NULL
#       }
#       
#     }
#     if (is.null(condf)) {condf <- proxdf} else {condf <- rbind(condf,proxdf)}
#   }
#   return(condf)
# }
# 
# #Process contact analysis for all pairs of individuals from group 1 to group 2
# twoGroup <- function(traj,traj2,dc,tc,GetSim)
#   #Get Pairs
#   pairs <- checkTO(traj,traj2)
# n.pairs <- nrow(pairs)
# condf <- NULL
# for (i in 1:n.pairs){
#   traj1 <- mtraj1[mt_track_id(traj)==pairs$ID1[i],]
#   traj2 <- mtraj2[mt_track_id(traj2)==pairs$ID2[i],]
#   
#   if (GetSim){
#     trajs <- GetSimultaneous(traj1,traj2,tc)
#     tr1 <- trajs[mt_track_id(trajs)==pairs$ID1[i],]
#     tr2 <- trajs[mt_track_id(trajs)==pairs$ID2[i],]
#     proxdf <- data.frame(id1=pairs$ID1[i],id2=pairs$ID2[i],row1=row.names(tr1),row2=row.names(tr2),dist=st_distance(tr1,tr2,by_element=TRUE),difftime=as.numeric(abs(mt_time(tr1)-mt_time(tr2))))
#     units(proxdf$difftime) <- as_units('s')
#     proxdf <- subset(proxdf,dist < dc & difftime < tc)
#   } else {
#     dM <- st_distance(traj1,traj2)
#     tM <- abs(-outer(as.numeric(mt_time(traj1)),as.numeric(mt_time(traj2)),'-'))
#     
#     units(tM) <- as_units('s')
#     units(dc) <- units(dM)
#     
#     rownames(tM) <- rownames(traj1)
#     colnames(tM) <- rownames(traj2)
#     
#     ind <- which(dM < dc & tM < tc, arr.ind=TRUE)
#     
#     if (length(ind) > 0){
#       rnm1 <- rownames(tM)[which(dM < dc & tM < tc, arr.ind = TRUE)[, 1]]
#       rnm2 <- colnames(tM)[which(dM < dc & tM < tc, arr.ind = TRUE)[, 2]]
#       proxdf <- data.frame(id1=pairs$ID1[i],id2=pairs$ID2[i],row1=rnm1,row2=rnm2,dist=dM[ind],difftime=tM[ind])
#     } else {
#       proxdf <- NULL
#     }
#   }
#   if (is.null(condf)) {condf <- proxdf} else {condf <- rbind(condf,proxdf)}
# }
# return(condf)
# }
