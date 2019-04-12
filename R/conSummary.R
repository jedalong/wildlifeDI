#new function for processing contacts from ContactAnalysis
# ---- roxygen documentation ----
#
#' @title Summarize contacts and phases
#' @description
#' Computes some basic summary statistics from a contact analysis.
#'
#' @details
#' This function is used following the \code{conPhase} function. It computes the following summary statistics from the contact analysis:
#' - total number of fixes in the dataset
#' - total number of fixes deemed a contact
#' - number of contact phases
#' - longest phase duration
#' - mean phase duration
#' - median phase duration
#' - no. of phase where the duration is only one fix (i.e., instantaneous contacts)

#' @param ltraj an object of the class \code{ltraj} which should be output from the function \code{conPhase}.
#'
#' @return
#' A dataframe that can be used to summarize contact phases.
#'
# @references
#'
#' @keywords contacts
#' @seealso conPhase
#' 
#' @examples 
#' \dontrun{
#' data(does)
#' doecons <- conProcess(does,tc=15*60,dc=50)
#' doephas <- conPhase(doecons,pc=60*60)
#' conSummary(doephas)
#' }
#' @export
#
# ---- End of roxygen documentation ----
conSummary <- function(ltraj){
  x <- ld(ltraj)
  outdf <- data.frame(stat = c('n Fixes','n Contacts','n Phases','Longest Phase (secs)','Mean Phase (secs)','Median Phase (secs)','no. one-fix Phases'), result=NA)
  outdf$result[1] <- dim(x)[1]                                 
  outdf$result[2] <- length(which(x$contacts > 0))
  outdf$result[3] <- length(unique(x$contact_pha,na.rm=T))
  
  #Longest segment duration
  #ind <- which(x$contact_pha == which.max(table(x$contact_pha)))
  #outdf$result[4] <- difftime(max(x$date[ind]),min(x$date[ind]),units='secs')
  
  #Mean segment duration
  dt <- NULL
  segid <- unique(x$contact_pha)[!is.na(unique(x$contact_pha))]
  for (seg in segid){
    ind <- which(x$contact_pha == seg)
    dt <- c(dt,difftime(max(x$date[ind]),min(x$date[ind]),units='secs'))
  }
  outdf$result[4] <- round(max(dt))
  outdf$result[5] <- round(mean(dt))
  outdf$result[6] <- median(dt)
  outdf$result[7] <- length(which(dt == 0))
  outdf
}