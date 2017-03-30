\name{mortGambAq}

\alias{mortGambAq}

\title{
Mortality model for the Anopheles gambiae s.s. aquatic stages
}

\description{
This is a temperature, and fractional dependent mortality model for the Anopheles gambiae s.s. in the aquatic stages
}

\usage{
mortGambAq(tslb, arab, gamb)
}

\arguments{
  \item{tslb}{Water temperature}
  \item{arab}{Number of An. arabiensis larvae}
  \item{gamb}{Number of An. gambiae s.s. larvae}
}

\details{

\code{mortGambAq} calculate mortaility rate according temperature and the number of An. arabiensis and An. gambiae s.s.
}

\value{
A vector with mortality rates in the aquatic stages
}

\references{

}

\author{
Torleif Markussen Lunde (Centre for International Health/Geophysical Institute, University of Bergen, Norway)
}

\note{}

\seealso{
}

\examples{
myMort <- mortGambAq(20:27, rep(1, length(20:27)), rep(1, length(20:27)))
}

\keyword{malaria}
