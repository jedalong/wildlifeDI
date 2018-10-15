#new function for processing contacts from ContactAnalysis
# ---- roxygen documentation ----
#
#' @title Examine context associated with contact phases
#'
#' @description
#' Extracts the variables associated with \code{var} before, during, and after contact phases, based on some specified time-lag.
#'
#' @details
#' This function is used following the \code{conphase} function. One should choose how to define the contact point (i.e., the parameter \code{contact}) depending on the research question. In most typical cases (with regular interval tracking data) the lag time should be set to the tracking interval and the gap should be set to 1/2 the tracking interval. 

#' @param ltraj an object of the class \code{ltraj} which should be output from the function \code{conPhase}.
#' @param var name(s) (as character) of columns (possibly from \code{infolocs}) to keep for contextual analysis.
#' @param contact how to define the point-of-contact. The default is to define it as all fixes in a phase \code{type = 'all'}, alternatively contacts can be defined as a single point along the phase defined as one of: \code{'first','last','minDist','minTime'}, which corresonds to the first fix int he contact phase, the last fix in the contact phase, the fix with the minimum time difference and the fix with the closest contact distance.
#' @param idcol column id associated with IDs of individuals, default is the 'burst'.
#' @param nrand  number of random fixes to be selected (default = 0).
#' @param nlag number of lags to compute in the before and after phases of a contact. If lag = 0 then only contacts are used.
#' @param lag time (in seconds) for defining the lags in before and after periods of a conatact. 
#' @param gap time (in seconds) for exluding the lags in before and after periods of a contact.  
#' @param phaid (optional) id(s) of the contact phase upon which to examine (default is all).
#'
#' @return
#' A dataframe that can be used to examine behaviour/context before, during, and after contact phases.
#'
# @references
#'
#' @keywords Contact Analysis
#' @seealso conPhase
# @examples
#' 
#' @export
#
# ---- End of roxygen documentation ----
conContext <- function(ltraj,var='dist',contact='all',idcol='burst',nrand=0,nlag=0,lag=0,gap=0,phaid){
  df <- ld(ltraj)
  if (missing(phaid)){ 
    phaid <- unique(df$contact_pha) 
    phaid <- phaid[!is.na(phaid)]
  }
  
  #Function to extract minTime and minDist from phases.
  funPhase <- function(phase, df, def){
    ind <- which(df$contact_pha == phase)
    if (def=='first'){
      i1 <- ind[which.min(df$date[ind])]    
    } else if (def=='last'){
      i1 <- ind[which.max(df$date[ind])]
    } else if (def=='minTime'){
      sub <- df[df$contact_pha == phase,]
      sub$id <- as.character(sub$id)
      sub$burst <- as.character(sub$burst)
      sub <- dl(sub)
      dfpairs <- conPairs(sub)
      i1 <- ind[dfpairs$contact_orig_rowid[which.min(dfpairs$contact_dt)]]
    } else if (def=='minDist'){
      sub <- df[df$contact_pha == phase,]
      sub$id <- as.character(sub$id)
      sub$burst <- as.character(sub$burst)
      sub <- dl(sub)
      dfpairs <- conPairs(sub)
      i1 <- ind[dfpairs$contact_orig_rowid[which.min(dfpairs$contact_d)]]
    } 
    return(i1)
  }
  
  #Before After analysis
  fun.BefAft <- function(phaid,df,var,nlag,lag,gap,idcol,contact){
    #subset data to only get individual associated with contact
    df.sub <- subset(df,get(idcol) == df[which(df$contact_pha == phaid),idcol][1])
    
    #get only the phase
    ind <- which(df.sub$contact_pha == phaid)
    pha <- df.sub[ind,]
    
    #Get start and end times of contact phase depending on how contacts are defined
    if (contact=='first'){
      i1 <- i2 <- ind[which.min(pha$date)]    
    } else if (contact=='last'){
      i1 <- i2 <- ind[which.max(pha$date)]
    } else if (contact=='minTime'){
      pha$id <- as.character(pha$id)
      pha$burst <- as.character(pha$burst)
      pha <- dl(pha)
      dfpairs <- conPairs(pha)
      i1 <- i2 <- ind[dfpairs$contact_orig_rowid[which.min(dfpairs$contact_dt)]]
    } else if (contact=='minDist'){
      pha$id <- as.character(pha$id)
      pha$burst <- as.character(pha$burst)
      pha <- dl(pha)
      dfpairs <- conPairs(pha)
      i1 <- i2 <- ind[dfpairs$contact_orig_rowid[which.min(dfpairs$contact_d)]]
    } else {                            #Use phases rather than contact points (default) 
      i1 <- ind[which.min(pha$date)]    #phase start index
      i2 <- ind[which.max(pha$date)]    #phase end index
    }
    
    cols <- c('date',var)
    #Get contact phase fixes 
    indpha <- seq(i1,i2)
    ccon <- data.frame(df.sub[indpha,cols])
    #names(ccon) <- var
    ccon$dt_con <- 0
    ccon$phase <- 'Con'
    outdf <- ccon
    
    
    #If nlag > 0 do the Before-After Analysis
    if (nlag > 0){
      for (l in 1:nlag){
        l1 <- gap+(l-1)*lag
        l2 <- gap+l*lag
        
        #Get before fixes
        ct1 <- difftime(df.sub$date[i1],df.sub$date)
        ibef <- which(ct1 > l1 & ct1 < l2)
        cbef <- data.frame(df.sub[ibef,cols]) 
        cbef$dt_con <- as.numeric(ct1[ibef])
        cbef$phase <- rep(paste('B',l,sep=''),length(ibef))
        
        #Get after fixes
        ct2 <- difftime(df.sub$date[i2],df.sub$date)
        iaft <- which(ct2 < -l1 & ct2 > -l2)
        caft <- data.frame(df.sub[iaft,cols])
        caft$dt_con <- as.numeric(-ct2[iaft])
        caft$phase <- rep(paste('A',l,sep=''),length(iaft))
        
        #append the rows
        outdf <- rbind(outdf,cbef,caft)
      }
    }
    
    outdf$phaid <- phaid
    outdf$ID <- df.sub[1,idcol]
    return(outdf)
  }
  
  #Get the contacts, and optionally the BefAft Phases
  dff <- do.call(rbind, lapply(phaid,fun.BefAft,df=df,var=var,nlag=nlag,lag=lag,gap=gap,idcol=idcol,contact=contact))
  #dff$phase <- factor(dff$phase,c('Bef','Con','Aft'))  #just for ordering later on
  
  #Get Random fixes if required
  if (nrand > 0){
    #Get Random fixes that are not 'contacts'
    ind2 <- which(df$contacts == 0)
    ind2 <- sample(ind2, nrand)
    rand <- data.frame(df[ind2,c('date',var)])
    rand$dt_con <- NA
    rand$phase <- 'Rnd'
    rand$phaid <- NA
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