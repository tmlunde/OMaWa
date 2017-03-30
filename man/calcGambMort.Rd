\name{calcGambMort}

\alias{calcGambMort}

\title{
Mortality model for the Anopheles gambiae complex
}

\description{
This is a temperature, relative humidity and size dependent mortality model for the Anopheles gambiae complex 
}

\usage{
calcGambMort(t2, RH, bovine, human, HBI, d = 0:199, mSiz = NULL)
}

\arguments{
  \item{t2}{Air temperature mosquitoes are exposed to}
  \item{RH}{Relative humidity mosquitoes are exposed to}
  \item{bovine}{The probability of finding bovine blood}
  \item{human}{The probability of finding human blood}
  \item{HBI}{Human blood index}
  \item{d}{Days at which the mortality curve should be estimated}
  \item{mSiz}{Mosquito size measured as wing length in mm}
}

\details{

\code{calcGambMort} calculate mortaility rate according to age groups. 
}

\value{
A list (\code{obj}) containg mortality rate for the different age groups. \code{obj$m1} is the mortality rate for age group 1, \code{obj$m2} for age group 2, ..., and \code{obj$m9} for age group 9.
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
Tair <- 25
RH <- 60
siz <- 3
prms <- calcGambMort(Tair, RH, 1, 1, 1, d = 0:199, mSiz = siz)
unlist(sapply(1:9, function (x) prms[paste("m", x, sep = "")]))

}

\keyword{malaria}
