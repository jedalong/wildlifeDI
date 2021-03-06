% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/HAI.R
\name{HAI}
\alias{HAI}
\title{Half-weight Association Index}
\usage{
HAI(traj1, traj2, OZ, tc = 0, dc = 50)
}
\arguments{
\item{traj1}{an object of the class \code{ltraj} which contains the time-stamped
movement fixes of the first object. Note this object must be a \code{type II
ltraj} object. For more information on objects of this type see \code{
help(ltraj)}.}

\item{traj2}{same as \code{traj1}.}

\item{OZ}{spatial polygon sf object associated with the home range (or some other form of) spatial overlap between \code{traj1} and \code{traj2}. Required to be an \code{sp} polygon object.}

\item{tc}{time threshold for determining simultaneous fixes -- see function: \code{GetSimultaneous}.}

\item{dc}{distance tolerance limit (in appropriate units) for defining when 
two fixes are spatially together.}
}
\value{
This function returns the numeric value of the HAI statistic. Values near 1 indicate
attraction within the shared home range area, while values near 0 indicate avoidance
within this shared area.
}
\description{
This function computes the Half-weight Association Index for examining the presence
of dynamic interaction in wildlife telemetry studies. This implementation follows
that outlined in the paper Atwood and Weeks (2003).
}
\details{
This function can be used to test for the presence of dynamic interaction within 
the shared area (often termed the overlap zone) of the two animals home ranges. Specifically,
HAI is calculated in identical fashion to that for \code{Ca}, but considers only those fixes 
in the shared area. Typically, the overlap zone (OZ) is easily obtained by taking the spatial
intersection of two polygon home ranges.
}
\examples{
\dontrun{
data(deer)
deer37 <- deer[1]
deer38 <- deer[2]
library(adehabitatHR)
library(sp)
library(rgeos)
#use minimum convex polygon for demonstration...
hr37 <- mcp(SpatialPoints(ld(deer37)[,1:2]))
hr38 <- mcp(SpatialPoints(ld(deer38)[,1:2]))
OZ <- gIntersection(hr37,hr38)
#tc = 7.5 minutes, dc = 50 meters
HAI(deer37, deer38, OZ=OZ, tc=7.5*60, dc=50)
}
}
\references{
Atwood, T.C. and Weeks Jr., H.P. (2003) Spatial home-range overlap and temporal
interaction in eastern coyotes: The influence of pair types and fragmentation.
\emph{Canadian Journal of Zoology}, \bold{81}: 1589-1597.\cr\cr
}
\seealso{
GetSimultaneous, Ca
}
\keyword{indices}
