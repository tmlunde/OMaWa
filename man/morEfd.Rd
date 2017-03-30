\name{morEfd}

\alias{morEfd}

\title{
Mordecai eggs laid per adult female per day
}

\description{
Calculates the number of eggs laid per adult female per day (Aedes albopictus)
}

\usage{
morEfd(T)
}

\arguments{
  \item{T}{Temperature in C}
}

\details{
R implemenation of the number of eggs laid per adult female per day calculation in Mordecai et al.
}

\value{
A vector of number of eggs laid per adult female per day
}

\references{
Delatte, H, Gimonneau, G, Triboire, A, and Fontenille, D: Influence of temperature on immature development, survival, longevity, fecundity, and gonotrophic cycles of Aedes albopictus, vector of chikungunya and dengue in the Indian Ocean., J Med Entomol 46(1), volume 46, 33–41, January 2009

Mordecai, Erin A, Paaijmans, Krijn P, Johnson, Leah R, Balzer, Christian, Ben-Horin, Tal, de Moor, Emily, McNally, Amy, Pawar, Samraat, Ryan, Sadie J, Smith, Thomas C, Lafferty, Kevin D, and Thrall, Peter: Optimal temperature for malaria transmission is dramatically lower than previously predicted., Ecol Lett, October 2012
}

\author{
Torleif Markussen Lunde (Centre for International Health/Geophysical Institute, University of Bergen, Norway)
}

\note{}

\seealso{
}

\examples{
morEfd(10:40)
}

\keyword{malaria}
