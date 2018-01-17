#new function for processing contacts from ContactAnalysis
# ---- roxygen documentation ----
#
#' @title Process contact segments
#'
#' @description
#' Computes segments where contacts occur based on a temporal tolerance.
#'
#' @details
#' This function is used following the \code{ContactAnalysis} function to arrange contacts into segments where continuous contact occurs (based on the user-defined time threshold \code{sc}. The idea is that we can consider a segment to be a continuous contact event (based on \code{dc} see \code{ContactAnalysis}) as long as the contact is only interrupted for no more than \col{sc} time units.

#' @param traj an object of the class \code{ltraj} which is output from the funciton \code{ContactAnalysis}.
#' @param sc time (in seconds) to allow for which to combine contact events (see details). 
#'
#' @return
#' An \code{ltraj} object with new column \code{contact_seg}.
#'
# @references
#'
#' @keywords Contact Analysis
#' @seealso dcPlot, contacts, conContext, conSummary
#' @examples
#' @export
#
# ---- End of roxygen documentation ----
#function to group contacts into interaction 'segments'.
conSegment <- function(ca,sc=0,idcol='burst'){
  df <- ld(ca)
  df$contact_seg <- 0
  df$contact_seg[which(df$contacts > 0)] <- 1
  
  seg.id <- 1
  
  ids <- unique(df[,idcol])
  for (i in ids){
    ind <- which(df[,idcol]==i)
    ##identify clusters of 'contacts'
    x <- df[ind,]
    run <- rle(x$contact_seg)
    len <- run$lengths
    val <- run$values
    
    #Combine clusters that are within 'tol=sc' time units into segments
    if (sc > 0){
      i1 <- 1
      for (j in 1:length(len)){
        i2 <- i1 + len[j] - 1
        if (val[j] == 0){
          dt <- abs(as.numeric(difftime(x$date[i1],x$date[i2],units='secs')))
          if(dt < sc){
            x$contact_seg[i1:i2] <- 1
          }
        } 
        i1 <- i2+1
      }
    }
    
    #re-identify the clusters of 'contacts' with gaps removed.
    run <- rle(x$contact_seg)
    len <- run$lengths
    val <- run$values
    seg.id2 <- seg.id+length(which(val>0))-1
    val[which(val>0)] <- seq(seg.id,seg.id2)
    seg.id <- seg.id2 + 1
    df$contact_seg[ind] <- rep(val,len)
  }
  
  df$contact_seg[which(df$contact_seg == 0)] <- NA
  return(dl(df))
}
