\name{morR0}

\alias{morR0}

\title{
Mordecai basic reproductive number
}

\description{
Calculates the temperature dependent basic reproductive number
}

\usage{
morR0(T, N, r)
}

\arguments{
  \item{T}{Temperature in C}
  \item{N}{Human density}
  \item{r}{Rate at which infected humans recover and acquire immunity}
}

\details{
R implemenation of the temperature dependent basic reproductive number calculation in Mordecai et al.
}

\value{
A vector the temperature dependent basic reproductive number
}

\references{
Mordecai, Erin A, Paaijmans, Krijn P, Johnson, Leah R, Balzer, Christian, Ben-Horin, Tal, de Moor, Emily, McNally, Amy, Pawar, Samraat, Ryan, Sadie J, Smith, Thomas C, Lafferty, Kevin D, and Thrall, Peter: Optimal temperature for malaria transmission is dramatically lower than previously predicted., Ecol Lett, October 2012
}

\author{
Torleif Markussen Lunde (Centre for International Health/Geophysical Institute, University of Bergen, Norway)
}

\note{}

\seealso{
}

\examples{
morR0(10:40)
}

\keyword{malaria}
