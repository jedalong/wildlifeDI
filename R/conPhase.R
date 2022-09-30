#new function for processing contacts from ContactAnalysis
# ---- roxygen documentation ----
#
#' @title Process contact phases
#'
#' @description
#' Computes phases where contacts occur based on a temporal tolerance.
#'
#' @details
#' This function is used following the \code{conProcess} function to arrange contacts into phases where continuous contact occurs (based on the user-defined time threshold \code{pc}. The idea is that we can consider a phase to be a continuous contact event (based on \code{dc} see \code{conProcess}) as long as the contact is only interrupted for no more than \code{pc} time units.

#' @param ltraj an object of the class \code{ltraj} which is output from the function \code{conProcess}.
#' @param pc time (in seconds) to allow for which to combine contact events (see details). 
#' @param idcol column used to identify individuals (default is the burst)
#'
#' @return
#' An \code{ltraj} object with new column \code{contact_pha}.
#'
#' @references
#'  Long, JA, Webb, SL, Harju, SM, Gee, KL (2022) Analyzing Contacts and Behavior from High Frequency 
#'  Tracking Data Using the wildlifeDI R Package. \emph{Geographical Analysis}. \bold{54}, 648--663.
#'
#' @keywords contacts
#' @seealso conProcess, conSpatial, conTemporal, conSummary
#' @examples
#' \dontrun{
#' data(does)
#' doecons <- conProcess(does,tc=15*60,dc=50)
#' doephas <- conPhase(doecons,pc=60*60)
#' }
#' 
#' @export
#
# ---- End of roxygen documentation ----
#function to group contacts into interaction 'phases'.
conPhase <- function(ltraj,pc=0,idcol='burst'){
  dfr <- ld(ltraj)
  dfr$contact_pha <- 0
  dfr$contact_pha[which(dfr$contacts > 0)] <- 1
  
  pha.id <- 1
  
  ids <- unique(dfr[,idcol])
  for (i in ids){
    ind <- which(dfr[,idcol]==i)
    ##identify clusters of 'contacts'
    x <- dfr[ind,]
    run <- rle(x$contact_pha)
    len <- run$lengths
    val <- run$values
    
    #Combine clusters that are within 'tol=sc' time units into phases
    if (pc > 0){
      i1 <- 1
      for (j in 1:length(len)){
        i2 <- i1 + len[j] - 1
        if (val[j] == 0){
          dt <- abs(as.numeric(difftime(x$date[i1],x$date[i2],units='secs')))
          if(dt < pc){
            x$contact_pha[i1:i2] <- 1
          }
        } 
        i1 <- i2+1
      }
    }
    
    #re-identify the clusters of 'contacts' with gaps removed.
    run <- rle(x$contact_pha)
    len <- run$lengths
    val <- run$values
    pha.id2 <- pha.id+length(which(val>0))-1
    val[which(val>0)] <- seq(pha.id,pha.id2)
    pha.id <- pha.id2 + 1
    dfr$contact_pha[ind] <- rep(val,len)
  }
  
  dfr$contact_pha[which(dfr$contact_pha == 0)] <- NA
  
  outtraj <- dl(dfr,proj4string=attr(ltraj,'proj4string'))
  return(outtraj)
}
