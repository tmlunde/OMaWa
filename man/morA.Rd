\name{morA}

\alias{morA}

\title{
Mordecai biting rate (mean oviposition time)^-1
}

\description{
Calculates the biting rate (mean oviposition time)^-1 for Anopheles pseudopunctipennis
}

\usage{
morA(T)
}

\arguments{
  \item{T}{Temperature in C}
}

\details{
R implemenation of the biting rate calculation in Mordecai et al.
}

\value{
A vector with biting rate
}

\references{
Lardeux, Frederic J, Tejerina, Rosenka H, Quispe, Vicente, and Chavez, Tamara K: A physiological time analysis of the duration of the gonotrophic cycle of Anopheles pseudopunctipennis and its implications for malaria transmission in Bolivia., Malar J 7, volume 7, 141, 2008

Mordecai, Erin A, Paaijmans, Krijn P, Johnson, Leah R, Balzer, Christian, Ben-Horin, Tal, de Moor, Emily, McNally, Amy, Pawar, Samraat, Ryan, Sadie J, Smith, Thomas C, Lafferty, Kevin D, and Thrall, Peter: Optimal temperature for malaria transmission is dramatically lower than previously predicted., Ecol Lett, October 2012
}

\author{
Torleif Markussen Lunde (Centre for International Health/Geophysical Institute, University of Bergen, Norway)
}

\note{}

\seealso{
}

\examples{
morA(10:40)
}

\keyword{malaria}
