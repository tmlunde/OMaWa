\name{morPea}

\alias{morPea}

\title{
Mordecai Egg-to-adult survival probability
}

\description{
Calculates the egg-to-adult survival probability for Anopheles gambiae s.s.
}

\usage{
morPea(T)
}

\arguments{
  \item{T}{Temperature in C}
}

\details{
R implemenation of the egg-to-adult survival probability calculation in Mordecai et al.
}

\value{
A vector of egg-to-adult survival probability 
}

\references{
Bayoh, M.N., and Lindsay, S.W.: Effect of temperature on the development of the aquatic stages of Anopheles gambiae sensu stricto (Diptera: Culicidae), Bulletin of Entomological Research 93, volume 93, 375–381, 2003

Mordecai, Erin A, Paaijmans, Krijn P, Johnson, Leah R, Balzer, Christian, Ben-Horin, Tal, de Moor, Emily, McNally, Amy, Pawar, Samraat, Ryan, Sadie J, Smith, Thomas C, Lafferty, Kevin D, and Thrall, Peter: Optimal temperature for malaria transmission is dramatically lower than previously predicted., Ecol Lett, October 2012
}

\author{
Torleif Markussen Lunde (Centre for International Health/Geophysical Institute, University of Bergen, Norway)
}

\note{}

\seealso{
}

\examples{
morPea(10:40)
}

\keyword{malaria}
