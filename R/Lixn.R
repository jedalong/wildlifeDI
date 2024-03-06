# ---- roxygen documentation ----
#
#' @title Minta's Spatial-temporal interaction statistics
#'
#' @description
#' The function \code{Lixn} measures dynamic interaction between two animals following
#' the methods outlined by Minta (1992).
#'
#' @details
#' The function \code{Lixn} can be used to calculate the Minta (1992) measures of dynamic
#' interaction between two animals. The Minta statistic tests how the two animals simultaneously utilize
#' an area shared between the two individuals. Three coefficients are produced \eqn{L_{AA}}, \eqn{L_{BB}}, 
#' and \eqn{L_{ixn}}. Each of these statistics are based on a contingency table that compares the observed 
#' frequency of those fixes that are simultaneous and within/outside the shared area to expectations based on 
#' area overlap proportions (if \code{method="spatial"}) or expectations derived from all fixes (if 
#' \code{method="frequency"}) -- see Minta (1992) for more details. A Chi-squared statistic can then
#' be used to examine the significance between the observed and expected use of the shared area.
#' \cr\cr
#' Minta (1992) suggests the following interpretations of the coefficients. When \eqn{L_{AA}}
#' is near 0, the first animal's use of the shared area is random (or as expected). When
#' \eqn{L_{AA} > 0} it signifies spatial attraction to the shared area, or greater than
#' expected use. When \eqn{L_{AA} < 0} it signifies spatial avoidance of the shared area, or
#' less than expected use. Interpretation of \eqn{L_{BB}} is the same as for \eqn{L_{AA}} with
#' respect to the second animal. \eqn{L_{ixn}} tells us far more about the nature of the
#' interaction between the two individuals. As \eqn{L_{ixn}} nears 0, both animals use the
#' shared area randomly, with regards to the other animal. If \eqn{L_{ixn} > 0} the animals
#' use the shared area more \emph{simultaneously}, whereas if \eqn{L_{ixn} < 0} it is an
#' indication of \emph{solitary} use, or avoidance. This is why \eqn{L_{ixn}} is termed the temporal
#' interaction coefficient. A Chi-squared test can be used to identify the significance
#' of the \eqn{L_{AA}}, \eqn{L_{BB}}, and \eqn{L_{ixn}} values.
#' \cr\cr
#' NOTEs: 
#' \cr
#' 1. With modern telemetry datasets, where home ranges are readily estimated, choosing \code{method = 'spatial'}
#' is most appropriate. If parmater \code{hr} is not specified, the code uses the minimum convex hull method to calculate individual 
#' home ranges.
#' \cr
#' 2. When the home ranges do not overlap the Lixn statistic is not defined and the function returns a 
#' string of NA's.
#' \cr
#' 3. When one home range completely encloses another the Lixn statistic is not defined and the function returns
#' a string of NA's and \code{'ContainsB'} (or \code{'ContainsB'}) under the p.IXN result.
#' \cr
#' 4. Further to points 2 and 3, the Lixn statistic is not appropriate in situations where the overlap area is 
#' either very large or very small relative to either home range (i.e., a situation with almost complete enclosure 
#' or virtually no overlap). The example data (deer) is an exampl of a near complete enclosure. Thus, it is advised 
#' that \code{Lixn} be used only in situations where there are suitable marginal areas for areaA, areaB, and areaAB 
#' -- see Minta (1992). 
#'
#' @param traj an object of the class \code{move2} which contains the time-stamped movement fixes of at least two individuals. For more information on objects of this type see \code{help(mt_as_move2)}.
#' @param traj2 (optional) same as traj, but for the second group of individuals. See \code{checkTO}
#' @param method method for computing the marginal distribution from which expected
#'    values are computed. If \code{method = "spatial"}, the marginal values are calculated based on areas
#'    of the shared and unshared portions of the home ranges. If \code{method = "frequency"}, the marginal 
#'    values are calculated based on the number of all fixes within the shared and unshared portions of 
#'    the home ranges -- see Details.
#' @param tc time threshold for determining simultaneous fixes -- see function: \code{GetSimultaneous}.
#' @param hr (optional) spatial polygon \code{sf} object associated with the home range (or some other form of) spatial range estimate for each individual in \code{traj}. The hr polygon should have a corresponding ID column with the same column name as in \code{traj}. If NULL (the default) the MCP home range estimate will be used for each individual.
#' @param OZ (-- required if method = 'frequency') A \code{sf} object representing the shared area polygon associated with spatial use overlap each pair of individuals in \code{traj}. Must be a \code{sf} polygon object and contain two columns id1 and id2 indicating the polygon associated with each pair.
#'
#' @return
#' This function returns a data.frame with values representing the calculated statistical values and associated \emph{p}-values from the Chi-squared test for each dyad.
#' \itemize{
#' \item pTable -- contingency table showing marginal probabilities of expected use 
#'    based on the selection of the \code{method} parameter.
#' \item nTable -- contingency table showing observed frequency of use of the shared 
#' area based on simultaneous fixes.
#' \item oTable -- the odds for each cell in the contingency table.
#' \item Laa -- the calculated value of the \eqn{L_{AA}} statistic
#' \item p.AA -- the associated \emph{p}-value
#' \item Lbb -- the calculated value of the \eqn{L_{BB}} statistic
#' \item p.BB -- the associated \emph{p}-value
#' \item Lixn -- the calculated value of the \eqn{L_{ixn}} statistic
#' \item p.IXN -- the associated \emph{p}-value
#' }
#'
#' @references
#' Minta, S.C. (1992) Tests of spatial and temporal interaction among animals.
#' \emph{Ecological Applications}, \bold{2}: 178-188
#'
#' @keywords indices
#' @seealso GetSimultaneous

#' @examples
#' \dontrun{
#' data(deer)
#' #tc = 7.5 minutes, dc = 50 meters
#' Lixn(deer,  method='spatial', tc=7.5*60)
#' 
#' #use internal buffer 500m of MCP for demonstration of frequency method
#' # NOTE: This is just an example, this is not an appropriate way to define overlap zone.
#' idcol <- mt_track_id_column(deer)
#' deercore <- deer |>
#'   st_union() |>
#'   st_convex_hull() |>
#'   st_buffer(-500)
#' Lixn(deer,  method='frequency', tc=7.5*60, OZ=deercore)
#' }
#' @export
#
# ---- End of roxygen documentation ----
Lixn <- function(traj,traj2,method="spatial",tc=0,hr=NULL,OZ=NULL){
  
  if (missing(traj2)){
    pairs <- checkTO(traj)
    pairs <- pairs[pairs$TO==TRUE,]
    mtraj <- traj
  } else {
    pairs <- checkTO(traj,traj2)
    pairs <- pairs[pairs$TO==TRUE,]
    if (st_crs(traj2) != st_crs(traj)){
      traj2 <- st_transform(traj2,crs=st_crs(traj))
    }
    mtraj <- data.frame(id = c(mt_track_id(traj),mt_track_id(traj2)),
                        time = c(mt_time(traj),mt_time(traj2)),
                        geometry = c(traj[[attr(traj,'sf_column')]],traj2[[attr(traj2,'sf_column')]])) |>
      st_as_sf(sf_column_name = "geometry", crs=st_crs(traj)) |>
      mt_as_move2(time_column='time',track_id_column='id')
  }

    n.pairs <- nrow(pairs)
  
  pairs$Laa <- NA
  pairs$p.aa <- NA
  pairs$Lbb <- NA
  pairs$p.bb <- NA
  pairs$Lixn <- NA
  pairs$p.ixn <- NA
  pairs$notes <- NA
  
  #get column name of ID column
  idcol <- mt_track_id_column(mtraj)
  
  #---- home range checks for SPATIAL method  ----
  if (method == 'spatial'){
    #IF no hr is specified
    if (is.null(hr)){
      hr <- mtraj |>
        dplyr::group_by_at(idcol) |>
        dplyr::summarise() |>
        st_convex_hull()
    } else {
      #Check if the polygon has a column with the correct name
      if (any(idcol %in% names(hr))==FALSE){
        print(paste0("hr object does not have correctly named id column: ",idcol))
        return(NULL)
      }
    }
  }
  
  #---- method specific computations for frequency method  ----
  if (method == 'frequency'){
    
    #check overlap zone
    chk3 <-'id1' %in% names(OZ)
    chk4 <-'id2' %in% names(OZ)
    if (chk3 == FALSE || chk4== FALSE) {
      print('There are no columns named id1 and id2 in the OZ. This is required for method = frequency.')
      return(NULL)
    }
  }
  
  #Loop through the pairs
  for (i in 1:n.pairs){
    traj1 <- mtraj[mt_track_id(traj)==pairs$ID1[i],]
    traj2 <- mtraj[mt_track_id(traj)==pairs$ID2[i],]
    
    trajs <- GetSimultaneous(traj1,traj2,tc)
    
    tr1 <- trajs[mt_track_id(trajs)==pairs$ID1[i],]
    tr2 <- trajs[mt_track_id(trajs)==pairs$ID2[i],]
    
    n1 <- nrow(traj1)
    n2 <- nrow(traj2)
    n <- nrow(tr1)
    
    #Method specific Calculations
    if (method == 'spatial'){
      
      hr1 <- hr[ hr[[idcol]] == pairs$ID1[i], ]
      hr2 <- hr[ hr[[idcol]] == pairs$ID2[i], ]
      
      #could check that hr1 and hr2 exist and are spatial polygons here.
      chk1 <- st_contains(hr1,hr2,sparse=FALSE)
      chk2 <- st_contains(hr2,hr1,sparse=FALSE)
      
      if (chk1 == TRUE){
        pairs$notes[i]<- 'ContainsA'
        next()
      } else if (chk2 == TRUE){
        pairs$notes[i]<- 'ContainsB'
        next()
      }
      areaA <- suppressWarnings(st_difference(hr1,hr2))
      areaB <- suppressWarnings(st_difference(hr2,hr1))
      areaAB <- suppressWarnings(st_intersection(hr1,hr2))
      
      if (is.null(areaAB)){ 
        pairs$notes[i] <- 'No overlap'
        next() 
      }
      #get intersected polygon for simultaneous fixes
      A1 <- st_intersects(tr1,areaA,sparse=FALSE)
      AB1 <- st_intersects(tr1,areaAB,sparse=FALSE)
      B2 <- st_intersects(tr2,areaB,sparse=FALSE)
      AB2 <- st_intersects(tr2,areaAB,sparse=FALSE)
      
      #check that the intersection vectors are same length before computing marginal values
      l.vec <- c(length(A1),length(B2),length(AB1),length(AB2))
      if (diff(range(l.vec)) == 0){
        #compute marginal values for simultaneous fixes
        n11 <- sum(AB1*AB2)
        n22 <- sum(A1*B2)
        n12 <- sum(A1*AB2)
        n21 <- sum(AB1*B2)
      } else {n11 <- 0; n12 <- 0; n21 <- 0; n22 <- 0}
      
      #compute expected values -- spatial
      a <- st_area(hr1)
      b <- st_area(hr2)
      ab <- st_area(areaAB)
      
      units(a) <- NULL
      units(b) <- NULL
      units(ab) <- NULL
      
      p11 <- (ab^2)/(a*b)
      p12 <- (1 - (ab/a))*(ab/b)
      p21 <- (ab/a)*(1 - (ab/b))
      p22 <- (1-(ab/a))*(1-(ab/b))
      
    } else if (method == 'frequency'){
      
      areaAB <- OZ[ OZ$id1 == pairs$ID1[i] & OZ$id2 == pairs$ID2[i] , ]
      areaA <- suppressWarnings(st_difference(st_convex_hull(st_union(tr1)),areaAB))
      areaB <- suppressWarnings(st_difference(st_convex_hull(st_union(tr2)),areaAB))
      
      #get intersected polygon for simultaneous fixes
      A1 <- st_intersects(tr1,areaA,sparse=FALSE)
      AB1 <- st_intersects(tr1,areaAB,sparse=FALSE)
      B2 <- st_intersects(tr2,areaB,sparse=FALSE)
      AB2 <- st_intersects(tr2,areaAB,sparse=FALSE)
      
      #check that the intersection vectors are same length before computing marginal values
      l.vec <- c(length(A1),length(B2),length(AB1),length(AB2))
      if (diff(range(l.vec)) == 0){
        #compute marginal values for simultaneous fixes
        n11 <- sum(AB1*AB2)
        n22 <- sum(A1*B2)
        n12 <- sum(A1*AB2)
        n21 <- sum(AB1*B2)
      } else {n11 <- 0; n12 <- 0; n21 <- 0; n22 <- 0}
      
      #get pts that intersect areaAB for ALL fixes
      rAB <- sum(AB1)
      sAB <- sum(AB2)
      p11 <- (rAB*sAB)/(n1*n2)
      p12 <- (1 - (rAB/n1))*(sAB/n2)
      p21 <- (rAB/n1)*(1-(sAB/n2))
      p22 <- (1 - (rAB/n1))*(1 - (sAB/n2))
      
    } # END OF METHOD SPECIFIC CALCULATIONS
    
    #Compute chi-square results 
    #compute summary statistics
    w <- (n11*n22)/(n12*n21)
    L <- log(w)
    se.L <- sqrt((1/n11) + (1/n12) + (1/n21) + (1/n22))
    #Chi-square of marginal frequency values
    n.1 <- n11 + n21
    n.2 <- n12 + n22
    n1. <- n11 + n12
    n2. <- n21 + n22 
    #compute summary and chi-values for Extrinsic Hypothesis.
    p.1 <- p11+p21
    p.2 <- p12+p22
    p1. <- p11+p12
    p2. <- p21+p22

    # Intrinsic Hypothesis
    #------------ NOT RETURNED -----------------------
    #chiINT <- (n*((n11*n22 - n12*n21)^2))/(as.numeric(n1.)*as.numeric(n2.)*as.numeric(n.1)*as.numeric(n.2))
    #phiINT <- sqrt(chiINT/n)
    #--------------------------------------------------------------------

    #see email from Eric Howe
    #Not returned, not sure value of Chi.tot
    chi.tot <- (((n11-p11*n)^2)/p11*n) + (((n12-p12*n)^2)/p12*n) + (((n21-p21*n)^2)/(p21*n)) + (((n22-p22*n)^2)/(p22*n))
    Laa <- log((p.2*n.1)/(p.1*n.2))
    chiAA <- ((n.1-p.1*n)^2)/(p.2*p.1*n)
    Lbb <- log((p2.*n1.)/(p1.*n2.))
    chiBB <- ((n1.-p1.*n)^2)/(p2.*p1.*n)
    Lixn <- log(((n11/p11)+(n22/p22))/((n12/p12)+(n21/p21)))
    chiIXN <- (((n11/p11) + (n22/p22) - (n12/p12) - (n21/p21))^2)/(n*((1/p11) + (1/p12) + (1/p21) + (1/p22)))
    #odds for each cell
    o11 <- n11/(p11*n)
    o12 <- n12/(p12*n)
    o21 <- n21/(p21*n)
    o22 <- n22/(p22*n)
    
    #compute chi-square p-values for easier interpretation.
    #Intrinsic Hypothesis - Not Returned
    #p.INT <- 1 - pchisq(chiINT,df=1)
    #p.tot <- 1 - pchisq(chi.tot,df=3)  #something off, see email from Eric Howe, how to interpret?
    
    #Extrinsic Hypothesis
    p.AA <- 1 - pchisq(chiAA,df=1)
    p.BB <- 1 - pchisq(chiBB,df=1)
    p.IXN <- 1 - pchisq(chiIXN,df=1)      #Fixed Typo Here...
    #contingency matrices not returned, but could be useful for something
    pTable <- matrix(c(p11,p12,p21,p22),ncol=2,byrow=T,dimnames=list(c("B","b"),c("A","b")))
    nTable <- matrix(c(n11,n12,n21,n22),ncol=2,byrow=T,dimnames=list(c("B","b"),c("A","b")))
    oTable <- matrix(c(o11,o12,o21,o22),ncol=2,byrow=T,dimnames=list(c("B","b"),c("A","b")))
    
    #Not Returning IntHyp as not easily interpreted and ExtHyp is better
    #IntHyp <- list(n=n,L=L,se.L=se.L,p.INT=p.INT,phi.INT=phiINT)
    #ExtHyp <- list(Laa=Laa,p.AA=p.AA,Lbb=Lbb,p.BB=p.BB,Lixn=Lixn,p.IXN=p.IXN)
    
    pairs$Laa[i] <- Laa
    pairs$p.aa[i] <- p.AA
    pairs$Lbb[i] <- Lbb
    pairs$p.bb[i] <- p.BB
    pairs$Lixn[i] <- Lixn
    pairs$p.ixn[i] <- p.IXN
     
  }

  return(pairs)
}
#====================End of Minta Function =====================================
