# ---- roxygen documentation ----
#
#' @title Extract contacts following contact analysis
#'
#' @description
#' Create a dataframe with only the attributes associated with contacts, and do some low-level processing.
#'
#' @details
#' This function is used to extract contacts from a trajectory following use of the \code{contacts} function. If \code{get='all'} then all contacts are returned. If \code{get='minDist' or 'minTime'}, then where multiple contacts occur for the same fix, only the contact that is nearest in space or time is returned.

#' @param mtraj an object of the class \code{ltraj} which is output from the funciton \code{contacts}.
#' @param get one of \code{c('all','minDist','minTime')} properties of the contacts to retreive (see Details). 
#'
#' @return
#' A data frame, with the fix attributes associated with contacts.
#'
# @references
#'
#' @keywords Contact Analysis
#' @seealso contacts
#' @examples
#' @export
#
# ---- End of roxygen documentation ----
conExtract <- function(ca,get = 'all',idcol='burst'){
  if (class(ca)[1]=='ltraj'){
    df <- ld(ca)
  }
  ind <- which(df$contacts > 0)
  outdf <- df[0,]
  for (i in ind){
    a <- df[i,]
    a <- a[rep(seq(nrow(a)),a$contacts), ] 
    a$contact_rowid1 <- i
    a$contact_id2 <- unlist(strsplit(as.character(df$contact_id[i]),','))
    a$contact_rowid2 <- unlist(strsplit(as.character(df$contact_rowid[i]),','))
    a$contact_d <- unlist(strsplit(as.character(df$contact_d[i]),','))
    a$contact_dt <- unlist(strsplit(as.character(df$contact_dt[i]),','))
    #If we only want the nearest contact in space
    if (get=='minDist'){
      a <- a[which.min(a$contact_d),]
    } else if (get == 'minTime'){
      a <- a[which.min(a$contact_dt),]
    }
    #save output dataframe
    outdf <- rbind(outdf,a)
  }
  cols <- c(idcol,'contact_rowid1','contact_id2','contact_rowid2','contact_d','contact_dt')
  outdf <- outdf[,cols]
  names(outdf) <- c('contact_id1','contact_rowid1','contact_id2','contact_rowid2','contact_d','contact_dt')
  return(outdf)
}