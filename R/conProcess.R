# ---- roxygen documentation ----
#' @title Process contacts
#'
#' @description
#' This function performs basic contact analysis between individuals in a group of tracked animals, or between two different groups of tracked animals.
#'
#' @details
#' This function can be used to identify the nature of contacts in space and time between individuals in one or two groups.

#' @param mtraj1 an object of the class \code{ltraj} which contains the time-stamped movement fixes of the first group of individuals. Each individual should be stored with a unique 'id'. (see \code{?as.ltraj})
#' @param mtraj2 (optional) same as \code{mtraj1}, but for the second group of individuals.
#' @param tc time threshold for determining simultaneous fixes -- see function: \code{GetSimultaneous}.
#' @param dc distance tolerance limit (in appropriate units) for defining when two fixes are spatially together.
#' @param idcol1 column id associated with IDs of the first group of individuals, default is the 'burst'.
#' @param idcol2 (optional) column id associated with IDs of the second group of individuals.
#'
#' @return
#' This function returns the object \code{mtraj1} with three additional fields:
#' contact - the number of contacts associated with each given fix.
#' contact_id - the id(s) of the individual(s) associated with those contacts.
#' contact_d - the distance (in the same units as \code{mtraj1}) at which the contacts occur.
#' Note that if more than one contact occurs at a given time, the contact_id and contact_d fields will be a concatenated list of the contact IDs and distances.
#'
# @references
#'
#' @keywords Contact Analysis
#' @seealso GetSimultaneous, dcPlot, conPhase, conSummary
#' @export
#
# ---- End of roxygen documentation ----

conProcess <- function(mtraj1,mtraj2,dc=0,tc=0,idcol1='burst',idcol2){
  
  #Process contact analysis for all pairs of individuals in one group
  oneGroup <- function(mtraj1,dc,tc,idcol1){
    #no mtraj2 specified, get all unique pairs of individuals
    df <- ld(mtraj1)
    col1 <- which(names(df)==idcol1)
    id1 <- as.character(unique(df[,col1]))
    #Get all the unique combinations between one group
    pairs <- expand.grid(id1,id1,stringsAsFactors=F)
    pairs <- pairs[-which(pairs$Var1 == pairs$Var2),]
    pairs <- pairs[order(pairs$Var1),]
    #pairs <- pairs[!duplicated(t(apply(pairs, 1, sort))),]
    
    #make new columns
    df$contacts <- 0
    df$contact_id <- NA
    df$contact_rowid <- NA
    df$contact_d <- NA
    df$contact_dt <- NA
    
    for (i in 1:(dim(pairs)[1])){
      i1 <- pairs$Var1[i]
      i2 <- pairs$Var2[i]
      ind1 <- which(df[,col1] == i1)
      ind2 <- which(df[,col1] == i2)
      tr1 <- dl(df[ind1,1:10])
      tr2 <- dl(df[ind2,1:10])  
      #Only do contact analysis if they overlap temporally.
      if (checkTO(tr1,tr2)$TO){
        #Proximity analysis - Use GetSimultaneous = FALSE to apply only to Focal Individual 
        p1 <- Prox(tr1,tr2,tc=tc,dc=dc,local=TRUE,GetSimultaneous=FALSE)
        #Get which fixes are deemed an contact
        pind1 <- which(p1$prox <= dc & p1$dt <= tc)
        indc1 <- ind1[pind1]
        indc2 <- ind2[which(row.names(ld(tr2)) %in% p1$row2[pind1])]
        #update count of contacts at a time (count is >1 if an animal has a contact with multiple animals at the same time)
        #Uses as.character to access row.names slot
        df[indc1,'contacts'] <- df[indc1,'contacts'] + 1
        #update contact id lists
        df$contact_id[indc1] <- sapply(df$contact_id[indc1],con,i2)
        #update contact distance
        df$contact_d[indc1] <- mapply(con,df$contact_d[indc1],p1$prox[pind1])
        #add contact time difference
        df$contact_dt[indc1] <- mapply(con,df$contact_dt[indc1],p1$dt[pind1])
        #add contact rowid
        if (length(indc1) != length(indc2)) {print(c(i1,i2))}
        df$contact_rowid[indc1] <- mapply(con,df$contact_rowid[indc1],indc2)
      }
    }
    return(df)
  }
  
  #Process contact analysis for all pairs of indivudals from group 1 to group 2
  twoGroup <- function(mtraj1,mtraj2,dc,tc,idcol1,idcol2){
    df <- ld(mtraj1)
    df2 <- ld(mtraj2)
    col1 <- which(names(df)==idcol1)
    col2 <- which(names(df)==idcol2)
    id1 <- as.character(unique(df[,col1]))
    id2 <- as.character(unique(df2[,col2]))
    #Get all the unique combinations between two groups
    pairs <- expand.grid(id1,id2,stringsAsFactors=F)
    
    #make new columns
    df$contacts <- 0
    df$contact_id <- NA
    df$contact_rowid <- NA
    df$contact_d <- NA
    df$contact_dt <- NA
    
    for (i in 1:(dim(pairs)[1])){
      i1 <- pairs$Var1[i]
      i2 <- pairs$Var2[i]
      ind1 <- which(df[,col1] == i1)
      ind2 <- which(df2[,col2] == i2)
      tr1 <- dl(df[ind1,1:10]) 
      tr2 <- dl(df2[ind2,1:10])  
      #Only do contact analysis if they overlap temporally.
      if (checkTO(tr1,tr2)$TO){
        #Proximity analysis - both ways
        p1 <- Prox(tr1,tr2,tc=tc,dc=dc,local=TRUE,GetSimultaneous=FALSE)
        #p2 <- Prox(tr2,tr1,tc=tc,dc=dc,local=TRUE,GetSimultaneous=FALSE)
        #Get which fixes are deemed an contact
        pind1 <- which(p1$prox <= dc & p1$dt <= tc)
        indc1 <- ind1[pind1]
        indc2 <- ind2[which(row.names(ld(tr2)) %in% p1$row2[pind1])]
        #pind2 <- which(p2$prox <= dc & p2$dt <= tc)
        #indc2 <- ind2[pind2]
        #update count of contacts at a time (count is >1 if an animal has a contact with multiple animals at the same time)
        df$contacts[indc1] <- df$contacts[indc1] + 1
        #update contact id lists
        df$contact_id[indc1] <- sapply(df$contact_id[indc1],con,i2)
        #update contact distance
        df$contact_d[indc1] <- mapply(con,df$contact_d[indc1],p1$prox[pind1])
        #add contact time difference
        df$contact_dt[indc1] <- mapply(con,df$contact_dt[indc1],p1$dt[pind1])
        #add contact rowid
        df$contact_rowid[indc1] <- mapply(con,df$contact_rowid[indc1],indc2)
      }
    }
    return(df)
  }
  
  #concatenate data frame function
  con <- function(a,b){
    if (is.na(a)){
      o <- b
    } else {
      o <- paste(a,b,sep=',')
    }
    o
  }
  
  #no idcol2 specified, use idcol1 for both.
  if (missing(idcol2)){ idcol2 <- idcol1 }

  if (missing(mtraj2)){
    mtraj2 <- NULL
    df <- oneGroup(mtraj1,dc,tc,idcol1)
  } else {
    df <- twoGroup(mtraj1,mtraj2,dc,tc,idcol1,idcol2)
  }
  
  outtraj <- dl(df,proj4string = attr(mtraj1,'proj4string'))
  return(outtraj)
}

