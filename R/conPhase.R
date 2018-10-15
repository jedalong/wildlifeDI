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

#' @param ltraj an object of the class \code{ltraj} which is output from the funciton \code{conProcess}.
#' @param pc time (in seconds) to allow for which to combine contact events (see details). 
#' @param idcol column used to identify indivduals (default is the burst)
#'
#' @return
#' An \code{ltraj} object with new column \code{contact_pha}.
#'
# @references
#'
#' @keywords Contact Analysis
#' @seealso conProcess, conSpatial, conTemporal, conSummary
#' @export
#
# ---- End of roxygen documentation ----
#function to group contacts into interaction 'phases'.
conPhase <- function(ltraj,pc=0,idcol='burst'){
  df <- ld(ltraj)
  df$contact_pha <- 0
  df$contact_pha[which(df$contacts > 0)] <- 1
  
  pha.id <- 1
  
  ids <- unique(df[,idcol])
  for (i in ids){
    ind <- which(df[,idcol]==i)
    ##identify clusters of 'contacts'
    x <- df[ind,]
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
    df$contact_pha[ind] <- rep(val,len)
  }
  
  df$contact_pha[which(df$contact_pha == 0)] <- NA
  return(dl(df))
}
