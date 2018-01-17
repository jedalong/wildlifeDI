# ---- roxygen documentation ----
#
#' @title Calculate net displacement from contacts
#' @description
#' Calcualte the net-displacement (distance) of fixes before and after a contact from that contact point.
#'
#' @details
#' This function is used to compute the net displacement away from contacts by an animal before and after a contact segment. Net displacement represents an important contextual variable, related to the mobility of the individual, to investigate during contact analysis. 

#' @param cs an object of the class \code{ltraj} which should be output from the function \code{conSegment}.
#' @param contact how to define the point-of-contact. The default is to define it as all fixes in a segment \code{type = 'all'}, alternatively contacts can be defined as a single point along the segment defined as one of: \code{'first','last','minDist','minTime'}, which corresonds to the first fix int he contact segment, the last fix in the contact segment, the fix with the minimum time difference and the fix with the closest contact distance.
#' @param idcol column id associated with IDs of individuals, default is the 'burst'
#' @param nrand  number of random fixes to be selected (default = 100, set to 0 for strictly before and after analysis).
#' @param nlag number of lags to compute in the before and after phases of a contact. If lag = 0 then only contacts are used.
#' @param lag time (in seconds) for defining the lags in before and after periods of a conatact. 
#' @param gap time (in seconds) for exluding the lags in before and after periods of a contact.  
#' @param segid (optional) id(s) of the contact segment upon which to examine (default is all).
#'
#' @return
#' A dataframe that can be used to summarize contact segments.
#'
# @references
#'
#' @keywords Contact Analysis
#' @seealso contacts, conSegment, conContext
#' @examples
#' @export
#
# ---- End of roxygen documentation ----
conDisplacement <- function(cs,contact='all',idcol='burst',nrand=100,nlag=0,lag=0,gap=0,segid){
  df <- ld(cs)
  if (missing(segid)){ 
    segid <- unique(df$contact_seg) 
    segid <- segid[!is.na(segid)]
  }
  
  #Compute displacement from contacts
  fun.BefAft <- function(segid,df,var,nlag,lag,gap,idcol,contact){
    #subset data to only get individual associated with contact
    df.sub <- subset(df,get(idcol) == df[which(df$contact_seg == segid),idcol][1])
    
    #get only the segment
    ind <- which(df.sub$contact_seg == segid)
    seg <- df.sub[ind,]
    
    #Get start and end times of contact segment depending on how contacts are defined
    if (contact=='first'){
      i1 <- i2 <- ind[which.min(seg$date)]    
    } else if (contact=='last'){
      it <- i2 <- ind[which.max(seg$date)]
    } else if (contact=='minTime'){
      i1 <- i2 <- ind[which.min(seg$contact_dt)]
    } else if (contact=='minDist'){
      i1 <- i2 <- ind[which.min(seg$contact_d)]
    } else {                            #Use segments rather than contact points (default) 
      i1 <- ind[which.min(seg$date)]    #segment start index
      i2 <- ind[which.max(seg$date)]    #segment end index
    }
    
    cols <- c('date')
    #Get contact segment fixes 
    indseg <- seq(i1,i2)
    ccon <- data.frame(df.sub[indseg,cols])
    names(ccon) <- cols
    ccon$displacement <- 0
    ccon$dt_con <- 0
    ccon$phase <- 'Con'
    outdf <- ccon
    
    spcon <- SpatialPoints(df.sub[indseg,c('x','y')])
    
    #If nlag > 0 do the Before-After Analysis
    if (nlag > 0){
      for (l in 1:nlag){
        l1 <- gap+(l-1)*lag
        l2 <- gap+l*lag
        #Get before fixes
        ct1 <- difftime(df.sub$date[i1],df.sub$date)
        ibef <- which(ct1 > l1 & ct1 < l2)
        cbef <- data.frame(df.sub[ibef,cols])
        if (dim(cbef)[1] > 0){
          names(cbef) <- cols
          spbef <- SpatialPoints(df.sub[ibef,c('x','y')])
          cbef$displacement <- as.vector(apply(gDistance(spcon,spbef,byid=T),1,min))
          cbef$dt_con <- as.numeric(ct1[ibef])
          cbef$phase <- rep(paste('B',l,sep=''),length(ibef))
          outdf <- rbind(outdf,cbef)
        }
        
        #Get after fixes
        ct2 <- difftime(df.sub$date[i2],df.sub$date)
        iaft <- which(ct2 < -l1 & ct2 > -l2)
        caft <- data.frame(df.sub[iaft,cols])
        if (dim(caft)[1] > 0){
          names(caft) <- cols
          spaft <- SpatialPoints(df.sub[iaft,c('x','y')])
          caft$displacement <- as.vector(apply(gDistance(spcon,spaft,byid=T),1,min))
          caft$dt_con <- as.numeric(-ct2[iaft])
          caft$phase <- rep(paste('A',l,sep=''),length(iaft))
          outdf <- rbind(outdf,caft)
        }
      }
    }
    
    outdf$segid <- segid
    outdf$ID <- df.sub[1,idcol]
    return(outdf)
  }
  
  #Get the contacts, and optionally the BefAft Phases
  dff <- do.call(rbind, lapply(segid,fun.BefAft,df=df,var=var,nlag=nlag,lag=lag,gap=gap,idcol=idcol,contact=contact))
  #dff$phase <- factor(dff$phase,c('Bef','Con','Aft'))  #just for ordering later on
  
  #Get Random fixes if required
  if (nrand > 0){
    #Get Random pairs of fixes that are not 'contacts'
    ind2 <- which(df$contacts == 0)
    ind2 <- sample(ind2, nrand)
    rand <- data.frame(df[ind2,c('date','dist')])
    names(rand) <- c('date','displacement')
    rand$dt_con <- NA
    rand$phase <- 'Rnd'
    rand$segid <- NA
    rand$ID <- df[ind2,idcol]
    dff <- rbind(dff,rand)
  }
  
  #Organize factor levels
  lev <- NULL
  for (i in nlag:1){
    lev <- c(lev,paste('B',i,sep=''))
  }
  lev <- c(lev,'Con')
  for (i in 1:nlag){
    lev <- c(lev,paste('A',i,sep=''))
  }
  if (nrand > 0) {lev <- c(lev,'Rnd')}
  
  dff$phase <- factor(dff$phase,levels=lev)
  return(dff)
}