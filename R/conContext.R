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
#' @param def how to define the point-of-contact. The default is to define it as all fixes in a phase \code{def = 'all'}, alternatively contacts can be defined as a single point along the phase defined as one of: \code{'first','last','minDist','minTime'}, which corresponds to the first fix int he contact phase, the last fix in the contact phase, the fix with the minimum time difference and the fix with the closest contact distance.
#' @param idcol column id associated with IDs of individuals, default is the 'burst'.
#' @param nrand  number of random fixes to be selected (default = 0).
#' @param nlag number of lags to compute in the before and after phases of a contact. If lag = 0 then only contacts are used.
#' @param lag time (in seconds) for defining the lags in before and after periods of a contact. 
#' @param gap time (in seconds) for excluding the lags in before and after periods of a contact.  
#' @param phaid (optional) id(s) of the contact phase upon which to examine (default is all).
#'
#' @return
#' A dataframe that can be used to examine behaviour/context before, during, and after contact phases.
#'
#' @references
#'  Long, JA, Webb, SL, Harju, SM, Gee, KL (2022) Analyzing Contacts and Behavior from High Frequency 
#'  Tracking Data Using the wildlifeDI R Package. \emph{Geographical Analysis}. \bold{54}, 648--663.
#'
#' @keywords contacts
#' @seealso conPhase
#' @examples
#' 
#' \dontrun{
#' data(does)
#' doecons <- conProcess(does,tc=15*60,dc=50)
#' doephas <- conPhase(doecons,pc=60*60)
#' cc <- conContext(var='dist',def='first',nlag=3,lag=30*60,gap=15*60)
#' head(cc)
#' }
#' 
#' @export
#
# ---- End of roxygen documentation ----
conContext <- function(ltraj,var='dist',def='all',idcol='burst',nrand=0,nlag=0,lag=0,gap=0,phaid){
  dfr <- ld(ltraj)
  if (missing(phaid)){ 
    phaid <- unique(dfr$contact_pha) 
    phaid <- phaid[!is.na(phaid)]
  }
  
  #Function to extract minTime and minDist from phases.
  funPhase <- function(phase, dfr, def){
    ind <- which(dfr$contact_pha == phase)
    if (def=='first'){
      i1 <- ind[which.min(dfr$date[ind])]    
    } else if (def=='last'){
      i1 <- ind[which.max(dfr$date[ind])]
    } else if (def=='minTime'){
      sub <- dfr[dfr$contact_pha == phase,]
      sub$id <- as.character(sub$id)
      sub$burst <- as.character(sub$burst)
      sub <- dl(sub)
      dfpairs <- conPairs(sub)
      i1 <- ind[dfpairs$contact_orig_rowid[which.min(dfpairs$contact_dt)]]
    } else if (def=='minDist'){
      sub <- dfr[dfr$contact_pha == phase,]
      sub$id <- as.character(sub$id)
      sub$burst <- as.character(sub$burst)
      sub <- dl(sub)
      dfpairs <- conPairs(sub)
      i1 <- ind[dfpairs$contact_orig_rowid[which.min(dfpairs$contact_d)]]
    } 
    return(i1)
  }
  
  #Before After analysis
  fun.BefAft <- function(phaid,dfr,var,nlag,lag,gap,idcol,def){
    #subset data to only get individual associated with contact
    df.sub <- subset(dfr,get(idcol) == dfr[which(dfr$contact_pha == phaid),idcol][1])
    
    #get only the phase
    ind <- which(df.sub$contact_pha == phaid)
    pha <- df.sub[ind,]
    
    #Get start and end times of contact phase depending on how contacts are defined
    if (def=='first'){
      i1 <- i2 <- ind[which.min(pha$date)]    
    } else if (def=='last'){
      i1 <- i2 <- ind[which.max(pha$date)]
    } else if (def=='minTime'){
      pha$id <- as.character(pha$id)
      pha$burst <- as.character(pha$burst)
      pha <- dl(pha)
      dfpairs <- conPairs(pha)
      i1 <- i2 <- ind[dfpairs$contact_orig_rowid[which.min(dfpairs$contact_dt)]]
    } else if (def=='minDist'){
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
    ccon$dt_lev <- 'Con'
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
        cbef$dt_con <- as.numeric(-ct1[ibef])
        cbef$dt_lev <- rep(paste('B',l,sep=''),length(ibef))
        
        #Get after fixes
        ct2 <- difftime(df.sub$date[i2],df.sub$date)
        iaft <- which(ct2 < -l1 & ct2 > -l2)
        caft <- data.frame(df.sub[iaft,cols])
        caft$dt_con <- as.numeric(-ct2[iaft])
        caft$dt_lev <- rep(paste('A',l,sep=''),length(iaft))
        
        #append the rows
        outdf <- rbind(outdf,cbef,caft)
      }
    }
    
    outdf$phaid <- phaid
    outdf$ID <- df.sub[1,idcol]
    return(outdf)
  }
  
  #Get the contacts, and optionally the BefAft Phases
  dff <- do.call(rbind, lapply(phaid,fun.BefAft,dfr=dfr,var=var,nlag=nlag,lag=lag,gap=gap,idcol=idcol,def=def))
  
  #Get Random fixes if required
  if (nrand > 0){
    #Get Random fixes that are not 'contacts'
    ind2 <- which(dfr$contacts == 0)
    ind2 <- sample(ind2, nrand)
    rand <- data.frame(dfr[ind2,c('date',var)])
    rand$dt_con <- NA
    rand$dt_lev <- 'Rnd'
    rand$phaid <- NA
    rand$ID <- dfr[ind2,idcol]
    dff <- rbind(dff,rand)
  }
  
  #Organize factor levels
  lev <- NULL
  if (nlag > 0){
    for (i in nlag:1){
      lev <- c(lev,paste('B',i,sep=''))
    }
  }
  lev <- c(lev,'Con')
  if (nlag > 0){
    for (i in 1:nlag){
      lev <- c(lev,paste('A',i,sep=''))
    } 
  }
  if (nrand > 0) {lev <- c(lev,'Rnd')}
  
  dff$dt_lev <- factor(dff$dt_lev,levels=lev)
  return(dff)
}