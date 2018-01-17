#new function for processing contacts from ContactAnalysis
# ---- roxygen documentation ----
#
#' @title Summarize contacts and segments
#' @description
#' Computes some basic summary statistics from a contact analysis.
#'
#' @details
#' This function is used following the \code{conSegment} function. It computes the following summary statistics from the contact analysis:
#' - total number of fixes in the dataset
#' - total number of fixes deemed a contact
#' - number of contact segments
#' - longest segment duration
#' - mean segment duration
#' - median segment duration
#' - no. of segments where the duration is only one fix (i.e., instantaneous contacts)

#' @param cs an object of the class \code{ltraj} which should be output from the function \code{conSegment}.
#'
#' @return
#' A dataframe that can be used to summarize contact segments.
#'
# @references
#'
#' @keywords Contact Analysis
#' @seealso contacts, conSegment
#' @examples
#' @export
#
# ---- End of roxygen documentation ----
conSummary <- function(cs){
  x <- ld(cs)
  outdf <- data.frame(stat = c('n Fixes','n Contacts','n Segments','Longest Segment (secs)','Mean Segment (secs)','Median Segment (secs)','no. one-fix segments'), result=NA)
  outdf$result[1] <- dim(x)[1]                                 
  outdf$result[2] <- length(which(x$contacts > 0))
  outdf$result[3] <- length(unique(x$contact_seg,na.rm=T))
  
  #Longest segment duration
  ind <- which(x$contact_seg == which.max(table(x$contact_seg)))
  outdf$result[4] <- difftime(max(x$date[ind]),min(x$date[ind]),units='secs')
  
  #Mean segment duration
  dt <- NULL
  segid <- unique(x$contact_seg)[!is.na(unique(x$contact_seg))]
  for (seg in segid){
    ind <- which(x$contact_seg == seg)
    dt <- c(dt,difftime(max(x$date[ind]),min(x$date[ind]),units='secs'))
  }
  outdf$result[5] <- round(mean(dt))
  outdf$result[6] <- median(dt)
  outdf$result[7] <- length(which(dt == 0))
  outdf
}