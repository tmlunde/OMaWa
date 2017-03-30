\name{calcArabMort}

\alias{calcArabMort}

\title{
Mortality model for the Anopheles arabiensis
}

\description{
This is a daily maximum temperature dependent mortality model for the Anopheles arabiensis
}

\usage{
calcArabMort(tMax, bovine, human, HBI)
}

\arguments{
  \item{tMax}{Air temperature mosquitoes are exposed to}
  \item{bovine}{The probability of finding bovine blood}
  \item{human}{The probability of finding human blood}
  \item{HBI}{Human blood index}
}

\details{

\code{calcArabMort} calculate mortaility rate according to age groups. 
}

\value{
A list (\code{obj}) containig mortality rate for the different age groups. \code{obj$m1} is the mortality rate for age group 1, \code{obj$m2} for age group 2, ..., and \code{obj$m9} for age group 9.
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
tMax <- 25
prms <- calcArabMort(tMax, 1, 1, 1)
unlist(sapply(1:9, function (x) prms[paste("m", x, sep = "")]))

}

\keyword{malaria}
