##Temporal Overlap Function

GetTemporalOverlap <- function(traj1,traj2,tc=0){
  #convert to dataframe
  tr1 <- ld(traj1)
  tr2 <- ld(traj2)
  
  #identify the temporal 
  t.min <- max(c(min(tr1$date),min(tr2$date))) - tc
  t.max <- min(c(max(tr1$date),max(tr2$date))) + tc
  
  #check to see if there is an overlap period
  if (t.min > t.max){stop('There is no temporal overlap.')}
  
  #find the fixes that are in the temporal overlap
  ind1 <- which(tr1$date >= t.min & tr1$date <= t.max)
  ind2 <- which(tr2$date >= t.min & tr2$date <= t.max)
  
  #convert to ltraj objects
  tr1.to <- dl(tr1[ind1,])
  tr2.to <- dl(tr2[ind2,])
  
  #Return the two ltraj objects
  return(c(tr1.to,tr2.to))
}

