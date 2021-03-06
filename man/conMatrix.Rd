% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/conMatrix.R
\name{conMatrix}
\alias{conMatrix}
\title{conMatrix}
\usage{
conMatrix(mtraj, idcol = "burst", output = "count")
}
\arguments{
\item{mtraj}{an object of the class \code{ltraj} which is output from the function \code{conProcess}.}

\item{idcol}{column id associated with IDs of the individuals, default is the 'burst'}

\item{output}{('count' or 'rate') whether to compute the counts of contacts in the contact matrix (default) or the contact-rate.}
}
\value{
A matrix, with the contact rates between individuals.
}
\description{
Create a matrix that can be used for social network analysis.
}
\details{
This function is used to calculate the contact rates between individuals and output them in the form of a matrix.
NOTE: This function is designed to be used only when a single ltraj object is input into the conProcess function. I.e., single-species contact networks.
}
\examples{
\dontrun{
data(does)
doecons <- conProcess(does,tc=15*60,dc=50)
doemat <- conMatrix(doecons)
doemat_rate <- conMatrix(doecons,output='rate')
}

}
\seealso{
conProcess, conPairs
}
\keyword{contacts}
