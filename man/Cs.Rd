% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Cs.r
\name{Cs}
\alias{Cs}
\title{Coefficient of Sociality}
\usage{
Cs(traj, traj2, tc = 0)
}
\arguments{
\item{traj}{an object of the class \code{move2} which contains the time-stamped movement fixes of at least two individuals. For more information on objects of this type see \code{help(mt_as_move2)}.}

\item{traj2}{(optional) same as traj, but for the second group of individuals. See \code{checkTO}}

\item{tc}{time threshold for determining simultaneous fixes -- see function: \code{GetSimultaneous}.}
}
\value{
This function returns a list of objects representing the calculated values from the
 Cs statistic and associated \emph{p}-values from the signed rank test.
 \itemize{
   \item Do -- The mean distance of simultaneous fixes.
   \item De -- The mean expected distance, from all fixes.
   \item Cs -- The coefficient of sociality, see \bold{Details}.
   \item p.Attract -- One sided \emph{p}-value from signed rank test, testing for attraction.
   \item p.Avoid -- One sided \emph{p}-value from signed rank test, testing for avoidance.
   }
}
\description{
The function \code{Cs} computes the coefficient of sociality between two moving
   objects following the methods outlined by Kenward et al. (1993). It also uses a
   signed Wilcoxon-rank test to test for significance.
}
\details{
This function can be used to calculate the Kenward et al. (1993) coefficient of sociality (Cs)
 between two animals. The Cs statistic tests the observed mean
 distance between simultaneous fixes against that expected by the overall
 distribution of distances between all fixes.
 \deqn{Cs=\frac{D_E-D_O}{D_O+D_E}}{Cs=(De-Do)/(De+Do)}
 Where \eqn{D_O}{Do} is the mean observed distance between simultaneous fixes, and \eqn{D_E}{De}
 is the mean expected distance between all fixes. Kenward et al. (1993) propose Cs
 as a useful metric for exploring attraction or avoidance behaviour.
 Values for Cs closer to 1 indicate
 attraction, while values for Cs closer to -1 indicate avoidance. Values of Cs
 near 0 indicate that the two animals' movements have no influence on one another.
 \cr \cr
 Further, the difference between the observed and expected distances are compared
 using a paired signed-rank test (both one-sided tests, indicative of attraction
 or avoidance). See the function \code{GetSimultaneous} for details on how
 simultaneous fixes are determined from two trajectories.
}
\examples{
data(deer)
#tc = 7.5 minutes
Cs(deer, tc = 7.5*60) 

}
\references{
Kenward, R.E., Marcstrom, V. and Karlbom, M. (1993) Post-nestling behaviour in
 goshawks, \emph{Accipiter gentilis: II}. Sex differences in sociality and nest-switching.
 \emph{Animal Behaviour}. \bold{46}, 371--378.
}
\seealso{
GetSimultaneous
}
\keyword{indices}
