\name{morP}

\alias{morP}

\title{
Mordecai Daily adult survival probability (p)
}

\description{
Calculates the daily adult survival probability (p) for Anopheles gambiae s.s.
}

\usage{
morP(T)
}

\arguments{
  \item{T}{Temperature in C}
}

\details{
R implemenation of the daily adult survival probability calculation in Mordecai et al.
}

\value{
A vector daily adult survival probability rate
}

\references{
Bayoh, Nabie: Studies on the development and survival of Anopheles gambiae sensu stricto at various temperatures and relative humidities, University of Durham, 2001

Mordecai, Erin A, Paaijmans, Krijn P, Johnson, Leah R, Balzer, Christian, Ben-Horin, Tal, de Moor, Emily, McNally, Amy, Pawar, Samraat, Ryan, Sadie J, Smith, Thomas C, Lafferty, Kevin D, and Thrall, Peter: Optimal temperature for malaria transmission is dramatically lower than previously predicted., Ecol Lett, October 2012
}

\author{
Torleif Markussen Lunde (Centre for International Health/Geophysical Institute, University of Bergen, Norway)
}

\note{}

\seealso{
}

\examples{
morP(10:40)
}

\keyword{malaria}
