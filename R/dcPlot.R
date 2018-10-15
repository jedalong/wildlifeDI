# ---- roxygen documentation ----
#
#' @title Contact distance plot
#'
#' @description
#' This function is an exploratory tool to examine the pairwise distances between individuals within a large telemetry dataset.
#'
#' @details
#' The dcPlot function can be used to study the frequency distribution of pairwise distances bewteen individual in a large telemetry dataset. It can be applied to a single group (if \code{mtraj2} is ignored) or two-groups of indidividuals. The code attempts to find natural breaks (local minima) in the frequency histogram using an approach based on the \code{peaks} function attributed to B. Ripley (see \url{https://stackoverflow.com/questions/6324354/add-a-curve-that-fits-the-peaks-from-a-plot-in-r} ). This tool is meant to be used for exploratory data analysis.

#' @param mtraj1 an object of the class \code{ltraj} which contains the time-stamped movement fixes of the first group of individuals. Each individual should be stored with a unique 'id'. (see \code{?as.ltraj})
#' @param mtraj2 (optional) same as \code{mtraj1}, but for the second group of individuals.
#' @param tc time threshold for determining simultaneous fixes -- see function: \code{GetSimultaneous}.
#' @param idcol1 column id associated with IDs of the first group of individuals, default is the 'burst'.
#' @param idcol2 (optional) column id associated with IDs of the second group of individuals.
#' @param histplot (logical) whether to output a histogram, along with a list of the natural breaks in the histogram (\code{histplot = TRUE}) or the dataframe of all paired distances used to construct the histogram (\code{histplot=FALSE}) to be used for further analysis.
#' @param dmax (optional) distance value to 'cut-off' the distance histogram.
#'
#' @return
#' If \code{histplot = TRUE} a list of the natural breaks (local minima) identified from the frequency histogram and a plot of the frequency histogram.
#' If \code{histplot = FALSE} a dataframe containing all the pairwise and simulataneous distances between all individuals in the trajectory dataset.
#'
# @references
#'
#' @seealso GetSimultaneous, conProcess, Prox, Don, IAB
#' @export
#
# ---- End of roxygen documentation ----
dcPlot <- function(mtraj1,mtraj2,tc=0,idcol1='burst',idcol2,histplot=TRUE,dmax){
  #Process contact analysis for all pairs of individuals in one group
  oneGroup <- function(mtraj1,tc,idcol1){
    #no mtraj2 specified, get all unique pairs of individuals
    df <- ld(mtraj1)
    col1 <- which(names(df)==idcol1)
    id1 <- as.character(unique(df[,col1]))
    #Get all the unique combinations between one group
    pairs <- expand.grid(id1,id1,stringsAsFactors=F)
    pairs <- pairs[-which(pairs$Var1 == pairs$Var2),]
    pairs <- pairs[order(pairs$Var1),]
    pairs <- pairs[!duplicated(t(apply(pairs, 1, sort))),]
    
    outdf <- NULL
    for (i in 1:(dim(pairs)[1])){
      i1 <- pairs$Var1[i]
      i2 <- pairs$Var2[i]
      ind1 <- which(df[,col1] == i1)
      ind2 <- which(df[,col1] == i2)
      tr1 <- dl(df[ind1,1:10])
      tr2 <- dl(df[ind2,1:10])  
      #Only do contact analysis if they overlap temporally.
      if (checkTO(tr1,tr2)$TO){
        #Proximity analysis
        p1 <- Prox(tr1,tr2,tc=tc,dc=0,local=TRUE,GetSimultaneous=TRUE)
        if (is.null(outdf)) {outdf <- p1} else {outdf <- rbind(outdf,p1)}
      }
    }
    return(outdf)
  }
  
  #Process contact analysis for all pairs of indivudals from group 1 to group 2
  twoGroup <- function(mtraj1,mtraj2,tc,idcol1,idcol2){
    df <- ld(mtraj1)
    df2 <- ld(mtraj2)
    col1 <- which(names(df)==idcol1)
    col2 <- which(names(df)==idcol2)
    id1 <- as.character(unique(df[,col1]))
    id2 <- as.character(unique(df2[,col2]))
    #Get all the unique combinations between two groups
    pairs <- expand.grid(id1,id2,stringsAsFactors=F)
    
    outdf <- NULL
    for (i in 1:(dim(pairs)[1])){
      i1 <- pairs$Var1[i]
      i2 <- pairs$Var2[i]
      ind1 <- which(df[,col1] == i1)
      ind2 <- which(df2[,col2] == i2)
      tr1 <- dl(df[ind1,1:10]) 
      tr2 <- dl(df2[ind2,1:10])  
      #Only do contact analysis if they overlap temporally.
      if (checkTO(tr1,tr2)$TO){
        #Proximity analysis
        p1 <- Prox(tr1,tr2,tc=tc,dc=0,local=TRUE,GetSimultaneous=TRUE)
        if (is.null(outdf)){ outdf <- p1 } else {outdf <- rbind(outdf,p1) }
      }
    }
    return(outdf)
  }
  
  #no idcol2 specified, use idcol1 for both.
  if (missing(idcol2)){ idcol2 <- idcol1 }
  
  if (missing(mtraj2)){
    mtraj2 <- NULL
    df <- oneGroup(mtraj1,tc,idcol1)
  } else {
    df <- twoGroup(mtraj1,mtraj2,tc,idcol1,idcol2)
  }
  
  if (histplot==TRUE){
    ###make dc plot here
    if (missing(dmax)) dmax <- max(df$prox)
    ind <- which(df$prox <= dmax)  
    df <- df[ind,]
    
    breaks = seq(0, max(df$prox), length.out=100) 
    cut = cut(df$prox, breaks, right=FALSE) 
    freq = table(cut)
    cumfreq0 = c(0, cumsum(freq))
    cumfreq0 <- cumfreq0 / max(cumfreq0)
    
    par(mar=c(4,4,1,4))
    
    hist <- hist(df$prox, breaks=breaks,col='grey',border='grey',axes=T, ann=T,prob=F,
                 xlim=c(0,max(df$prox)),xlab='Distance',ylab='Frequency',main='')
    par(new=T)
    plot(breaks,cumfreq0,type='l',axes=F,ann=F)
    axis(4, at=c(0, 0.5, 1))
    mtext('Cumulative Frequency',4,line=3)
    box()
    
    peaks<-function(series,span){
      z <- embed(series, span)
      s <- span%/%2
      v<- max.col(z) == 1 + s
      result <- c(rep(FALSE,s),v)
      result <- result[1:(length(result)-s)]
      result
    } 
    
    y <- -1*hist$counts
    pks <- hist$mids[which(peaks(y, span=11))] #span should be odd numbered
    pks. <- round(pks)
    abline(v=pks,lty=3,col='red')
    pky <- seq(0.1,0.7,length.out=length(pks))
    text(x=pks,pky, labels=pks.,pos=4,offset=0.1)
    return(pks)
  } else {
    return(df)
  }
}
