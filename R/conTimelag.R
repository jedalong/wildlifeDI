# ---- roxygen documentation ----
#
#' @title Compute time-lags from contact phases
#'
#' @description
#' Computes the time-lag from the nearest contact phase.
#'
#' @details
#' This function is used following the \code{conphase} function. One should choose how to define the contact point (i.e., the parameter \code{contact}) depending on the research question. 

#' @param traj an object of the class \code{move2} which should be output from the function \code{conPhase}.
#' @param def how to define the point-of-contact. The default is to define it as all fixes in a phase \code{def = 'all'}, alternatively contacts can be defined as a single point along the phase defined as one of: \code{'first','last','minDist','minTime'}, which corresponds to the first fix int he contact phase, the last fix in the contact phase, the fix with the minimum time difference and the fix with the closest contact distance.
#'
#' @return
#' A move2 object with an additional column contact_timelag with the time to the nearest (in time) contact phase. Negative values indicate times prior to the nearest contact phase and postive values indicate times after the nearest contact phase. If an individual has no contacts, the contact time-lag is NA.
#'
#' @references
#'  Long, JA, Webb, SL, Harju, SM, Gee, KL (2022) Analyzing Contacts and Behavior from High Frequency 
#'  Tracking Data Using the wildlifeDI R Package. \emph{Geographical Analysis}. \bold{54}, 648--663.
#'
#' @keywords contacts
#' @seealso conPhase, conProcess, conDisplacement
#' @examples
#' 
#' \dontrun{
#' data(does)
#' doecons <- conProcess(does,tc=15*60,dc=50)
#' doephas <- conPhase(doecons,pc=60*60)
#' conTL_first <- conTimelag(doephas,def='first')
#' conTL_all <- conTimelag(doephas,def='all')
#' }
#' 
#' @export
#
# ---- End of roxygen documentation ----

### JED STILL TO FIX
conTimelag <- function(traj,def='all'){
  
  # Set up contact timelag column
  traj$contact_timelag <- NA
  
  #Peform timelag to contact individually for every Animal.
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
    traj$contact_timelag[row.names(traj) %in% row.names(trj)] <- mt_time(trj) - mt_time(trj[j,])
  }
  return(traj)
  
}