\name{arabExtraDev}

\alias{arabExtraDev}

\title{
Modified development rate per day for Anopheles arabiensis in the aquatic stages dependent of the fraction of Anopheles gambiae s.s. and Anopheles arabiensis
}

\description{
Calculates the number of extra days required to develop from one stage to the next in the aquatic stages.
}

\usage{
arabExtraDev(arab, gamb)
}

\arguments{
  \item{arab}{Number of An. arabiensis larva}
  \item{gamb}{Number of An. gambiae s.s. larva}
}

\details{

\code{arabExtraDev} calculate a modification constant for the development rate in the aquatic stages.
}

\value{
A vector with modification constants
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
devRatePupaToAdultAt20 <- 1/aqDev(20, 6)
devRatePupaToAdultAt20Mod <- devRatePupaToAdultAt20 * arabExtraDev(1, 10)
}

\keyword{malaria}
