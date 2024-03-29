# ---- roxygen documentation ----
#
#' @title Calculate net displacement from contacts
#' @description
#' Calculate the net-displacement (distance) of fixes before and after a contact phase from the nearest contact pahse in time.
#'
#' @details
#' This function is used to compute the net displacement away from contacts by an animal before and after a contact phase. Net displacement represents an important variable related to the movement of the individual.  

#' @param traj an object of the class \code{move2} which should be output from the function \code{conPhase}.
#' @param def how to define the point-of-contact. The default is to define it as all fixes in a phase \code{type = 'all'}, alternatively contacts can be defined as a single point along the phase defined as one of: \code{'first','last','minDist','minTime'}, which corresponds to the first fix int he contact phase, the last fix in the contact phase, the fix with the minimum time difference and the fix with the closest contact distance.
#'
#' @return
#' An move2 object with a new 'contact_displacement' column indicating the straight-line distance to the nearest (in time) contact phase (defined using parameter def). If there are no contacts associated with an individual the contact displacement is NA.
#'
#' @references
#'  Long, JA, Webb, SL, Harju, SM, Gee, KL (2022) Analyzing Contacts and Behavior from High Frequency 
#'  Tracking Data Using the wildlifeDI R Package. \emph{Geographical Analysis}. \bold{54}, 648--663.
#'
#' @keywords contacts
#' @seealso conProcess, conPhase, conTimelag
#' 
#' @examples 
#' \dontrun{
#' data(does)
#' doecons <- conProcess(does,tc=15*60,dc=50)
#' doephas <- conPhase(doecons,pc=60*60)
#' disp_f <- conDisplacement(doephas,def='first')
#' disp_l <- conDisplacement(doephas,def='last')
#' }
#' 
#' @export
#
# ---- End of roxygen documentation ----

#ASSUMES PROJECTED COORDINATES
conDisplacement <- function(traj,def='all'){
  
  # Set up displacement analysis
  traj$contact_displacement <- NA
  
  #Peform Displacement individually for every Animal.
  anid <- unique(mt_track_id(traj))

  for (i in anid){
    trj <- traj[mt_track_id(traj)==i,]
    phas <- unique(trj$contact_pha)
    phas <- phas[!is.na(phas)]
    if (length(phas) == 0) { 
      #print(paste0('no contacts: ',i))
      next 
    }
    cid <- NULL
    for (p in phas){
      sub_rows <- which(trj$contact_pha == p)
      sub <- trj[sub_rows,]
      if (def == 'first'){
        cid <- c(cid,first(row.names(sub)))
      } else if (def == 'last'){
        cid <- c(cid,last(row.names(sub)))
      } else if (def == 'minTime'){
        cid <- c(cid,row.names(sub[which.min(sub$contact_dt),]))
      } else if (def == 'minDist') {
        cid <- c(cid,row.names(sub[which.min(sub$contact_d),]))
      } else {
        cid <- c(cid,row.names(sub))
      }
    }
    
    #CID are all the phase reference fixes row.names
    # Find the nearest phase reference fix for every fix in the animal
    refTimes <- mt_time(trj)[row.names(trj) %in% cid]
    j <-  sapply(mt_time(trj), function(x,refTimes,cid) {cid[which.min(abs(x-refTimes))]}, refTimes, cid)
    
    #compute distance from each fix to reference fix
    traj$contact_displacement[row.names(traj) %in% row.names(trj)] <- st_distance(trj,trj[j,],by_element=TRUE)
  }
  
  return(traj)
}