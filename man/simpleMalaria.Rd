\name{simpleMalaria}

\alias{simpleMalaria}

\title{
A simple age structured model describing transmission of malaria from humans to mosquitoes
}

\description{
This function shows how mosquitoes are become infected with malaria depending on the temperature, relative humidity, mosquito size and the fraction of infectious humans.
}

\usage{
simpleMalaria(Tair, RH, siz, fracHumInf, times)
}

\arguments{
  \item{Tair}{Air temperature mosquitoes are exposed to}
  \item{RH}{Relative humidity mosquitoes are exposed to}
  \item{siz}{Mosquito size measured as wing length (mm)}
  \item{fracHumInf}{Fraction of infectious humans}
  \item{times}{Time model should be evaluated in days}
}

\details{

\code{simpleMalaria} calculate the number of suceptible, infected and infectous mosquitoes at times. The initial mosquito population is 1000 individuals at time zero. 
The simple model includes temperature EIP and gonotrophic cycle, as well as temperature, relative humidity, and size dependent mortality.
}

\value{
\code{deSolve} object (\code{obj}) with \code{100*rowSums(obj[,19:26])/1000} is the fraction of infectious mosquitoes at time times, \code{100*rowSums(obj[,11:18])/1000} is the fraction of infected mosquitoes at time times, and \code{100*rowSums(obj[,2:9])/1000} is the fraction of suceptible mosquitoes at time times. 
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
plot(100*rowSums(simpleMalaria(Tair=25, fracHumInf=.01)[,19:26])/1000, type = "l")

}

\keyword{malaria}
