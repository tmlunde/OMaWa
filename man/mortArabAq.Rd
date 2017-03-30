\name{mortArabAq}

\alias{mortArabAq}

\title{
Mortality model for the Anopheles arabiensis aquatic stages
}

\description{
This is a temperature, and fractional dependent mortality model for the Anopheles arabiensis in the aquatic stages
}

\usage{
mortArabAq(tslb, arab, gamb)
}

\arguments{
  \item{tslb}{Water temperature}
  \item{arab}{Number of An. arabiensis larvae}
  \item{gamb}{Number of An. gambiae s.s. larvae}
}

\details{

\code{mortArabAq} calculate mortaility rate according temperature and the number of An. arabiensis and An. gambiae s.s.
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
myMort <- mortArabAq(20:27, rep(1, length(20:27)), rep(1, length(20:27)))
}

\keyword{malaria}
