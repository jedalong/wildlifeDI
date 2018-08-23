# ---- roxygen documentation ----
#
#' @title conDuration
#'
#' @description
#' Create a summary dataframe of the timing and and duration of contact segments.
#'
#' @details
#' This function is used to calculate the start and end times of contact segments, and their duration following use of the \code{conSegments} function.

#' @param mtraj an object of the class \code{ltraj} which is output from the funciton \code{conSegments}.
#' @param units units of duration e.g., \code{'mins'} (see \code{difftime}). 
#'
#' @return
#' A data frame, with the time and duration attributes associated with contact segments.
#'
# @references
#'
#' @keywords Contact Analysis
#' @seealso conSegment
#' @examples
#' @export
#
# ---- End of roxygen documentation ----

conDuration <- function(traj,units='auto'){
  df <- ld(traj)
  df <- df[df$contacts>0,]
  segs <- unique(df$contact_seg)
  outdf <- data.frame(contact_seg = segs, start_time = df$date[1], end_time = df$date[1])

  for (i in 1:length(segs)){
    temp <- df[df$contact_seg==segs[i],]
    outdf$start_time[i] <- as.POSIXct(min(temp$date,na.rm=TRUE)) 
    outdf$end_time[i] <- as.POSIXct(max(temp$date,na.rm=TRUE)) 
    
  }
  outdf$duration <- difftime(outdf$end_time,outdf$start_time,units=units)
  return(outdf)
}
