# ---- roxygen documentation ----
#
#' @title conTemporal
#'
#' @description
#' Create a summary dataframe of the timing and and duration of contact phases.
#'
#' @details
#' This function is used to calculate the start and end times of contact phases, and their duration following use of the \code{conPhase} function.

#' @param mtraj an object of the class \code{ltraj} which is output from the funciton \code{conPhase}.
#' @param units units of duration e.g., \code{'mins'} (see \code{difftime}). 
#'
#' @return
#' A data frame, with the time and duration attributes associated with contact phases.
#'
# @references
#'
#' @keywords Contact Analysis
#' @seealso conPhase
#' @examples
#' @export
#
# ---- End of roxygen documentation ----

conTemporal <- function(traj,units='auto'){
  df <- ld(traj)
  df <- df[df$contacts>0,]
  phas <- unique(df$contact_pha)
  outdf <- data.frame(contact_pha = phas, start_time = df$date[1], end_time = df$date[1])

  for (i in 1:length(phas)){
    temp <- df[df$contact_pha==phas[i],]
    outdf$start_time[i] <- as.POSIXct(min(temp$date,na.rm=TRUE)) 
    outdf$end_time[i] <- as.POSIXct(max(temp$date,na.rm=TRUE)) 
    
  }
  outdf$duration <- difftime(outdf$end_time,outdf$start_time,units=units)
  return(outdf)
}
