# ---- roxygen documentation ----
#
#' @title Calculate net displacement from contacts
#' @description
#' Calcualte the net-displacement (distance) of fixes before and after a contact from that contact point.
#'
#' @details
#' This function is used to compute the net displacement away from contacts by an animal before and after a contact phase. Net displacement represents an important contextual variable, related to the mobility of the individual.  

#' @param ltraj an object of the class \code{ltraj} which should be output from the function \code{conPhase}.
#' @param def how to define the point-of-contact. The default is to define it as all fixes in a phase \code{type = 'all'}, alternatively contacts can be defined as a single point along the phase defined as one of: \code{'first','last','minDist','minTime'}, which corresonds to the first fix int he contact phase, the last fix in the contact phase, the fix with the minimum time difference and the fix with the closest contact distance.
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
conDisplacement <- function(ltraj,def='all',idcol='burst'){
  df <- ld(ltraj)
  df$displacement <- 0
  
  #Compute displacement for each individual
  ind_displacement <- function(id,idcol, df, def){
    ind <- which(df[,idcol] == id)
    df <- df[ind,]
    n <- dim(df)[1]
    if (def == 'first') {
      cid <- which(!duplicated(df$contact_pha))
      cid <- cid[which(!is.na(df$contact_pha[cid]))] # remove NA
    } else if (def == 'last') {
      cid <- which(!duplicated(rev(df$contact_pha)))
      cid <- cid[which(!is.na(rev(df$contact_pha)[cid]))] # remove NA
      cid <- rev(1:n)[cid]
    } else if (def == 'minDist'){
      cid <- NULL
      cpdf <- conPairs(ltraj)
      for (i in 1:max(cpdf$contact_pha)){
        cpdf_sub <- subset(cpdf, contact_pha == i)
        cid <- c(cid, cpdf_sub$contact_orig_rowid[which.min(cpdf_sub$contact_d)])
      }
    } else if (def == 'minTime'){
      cid <- NULL
      cpdf <- conPairs(ltraj)
      for (i in 1:max(cpdf$contact_pha)){
        cpdf_sub <- subset(cpdf, contact_pha == i)
        cid <- c(cid, cpdf_sub$contact_orig_rowid[which.min(cpdf_sub$contact_dt)])
      }
    } else {
      cid <- which(df$contacts > 0)
    }
    
    disp <- rep(0,n)
    if (length(cid) > 0){
      for (i in 1:n){
        j <- cid[which.min(abs(difftime(df$date[cid],df$date[i],units='secs')))]
        disp[i] <- sqrt( (df$x[i]-df$x[j])^2 + (df$y[i]-df$y[j])^2 )
      }
    }
    return(disp)
  }
  
  ids <- levels(df[,idcol])
  for (id in ids){
    df$displacement[which(df[,idcol]==id)] <- ind_displacement(id,idcol,df,def)
  }  

  outtraj <- dl(df,proj4string=attr(ltraj,'proj4string'))
  return(outtraj)
}