\name{morBc}

\alias{morBc}

\title{
Mordecai Vector competence
}

\description{
Calculates the vector competence for Anopheles quadrimaculatus
}

\usage{
morBc(T)
}

\arguments{
  \item{T}{Temperature in C}
}

\details{
R implemenation of the vector competence calculation in Mordecai et al.
}

\value{
A vector with biting rate
}

\references{
Stratman-Thomas, W.K. (1940). The influence of temperature on Plasmodium vivax. Am. J. Trop. Med. Hyg., 20, 703–715.

Mordecai, Erin A, Paaijmans, Krijn P, Johnson, Leah R, Balzer, Christian, Ben-Horin, Tal, de Moor, Emily, McNally, Amy, Pawar, Samraat, Ryan, Sadie J, Smith, Thomas C, Lafferty, Kevin D, and Thrall, Peter: Optimal temperature for malaria transmission is dramatically lower than previously predicted., Ecol Lett, October 2012
}

\author{
Torleif Markussen Lunde (Centre for International Health/Geophysical Institute, University of Bergen, Norway)
}

\note{}

\seealso{
}

\examples{
morBc(10:40)
}

\keyword{malaria}
