\name{aqDev}

\alias{aqDev}

\title{
Development time in days for Anopheles gambiae s.l. in the aquatic stages
}

\description{
Calculates the number of days required to develop from one stage to the next in the aquatic stages. Valid for Anopheles gambiae s.l.
}

\usage{
aqDev(Temperature, stage)
}

\arguments{
  \item{Temperature}{Water temperature in C}
  \item{stage}{Which stage}
}

\details{

\code{aqDev} calculate development time in days according temperature. \code{stage = 1} is from egg to 1st instar larva, \code{stage = 2}, from 1st instar larva to 2st instar larva, ...,
\code{stage = 6} from pupa to adult.
}

\value{
A vector with development time in days
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
}

\keyword{malaria}
