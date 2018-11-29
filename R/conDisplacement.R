# ---- roxygen documentation ----
#
#' @title Calculate net displacement from contacts
#' @description
#' Calcualte the net-displacement (distance) of fixes before and after a contact from that contact point.
#'
#' @details
#' This function is used to compute the net displacement away from contacts by an animal before and after a contact phase. Net displacement represents an important contextual variable, related to the mobility of the individual.  

#' @param ltraj an object of the class \code{ltraj} which should be output from the function \code{conPhase}.
#' @param contact how to define the point-of-contact. The default is to define it as all fixes in a phase \code{type = 'all'}, alternatively contacts can be defined as a single point along the phase defined as one of: \code{'first','last','minDist','minTime'}, which corresonds to the first fix int he contact phase, the last fix in the contact phase, the fix with the minimum time difference and the fix with the closest contact distance.
#' @param idcol column id associated with IDs of individuals, default is the 'burst'
#'
#' @return
#' An ltraj object with a new 'displacement' column in infolocs.
#'
# @references
#'
#' @keywords Contact Analysis
#' @seealso conPhase, conContext
#' @export
#
# ---- End of roxygen documentation ----
conDisplacement <- function(ltraj,contact='all',idcol='burst',nrand=0,nlag=0,lag=0,gap=0,phaid){
  df <- ld(ltraj)
  n <- dim(df)[1]
  df$displacement <- 0
  cid <- which(df$contact == 1)
  for (i in 1:n){
    j <- cid[which.min(difftime(df$date[i],df$date[cid],units='secs'))]
    df$displacement[i] <- sqrt( (df$x[i]-df$x[j])^2 + (df$y[i]-df$y[j])^2 )
  }
  
  outtraj <- dl(df,proj4string=attr(ltraj,'proj4string'))
  return(outtraj)
}