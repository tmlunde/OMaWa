\name{morPdr}

\alias{morPdr}

\title{
Mordecai parasite development rate (PDR)
}

\description{
Calculates the parasite development rate (PDR) (An. gambiae, Anopheles culicifacies, Anopheles stephensi, An. quadrimaculatus, Anopheles atroparvus)
}

\usage{
morPdr(T)
}

\arguments{
  \item{T}{Temperature in C}
}

\details{
R implemenation of the parasite development rate calculation in Mordecai et al.
}

\value{
A vector of parasite development rate
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
morPdr(10:40)
}

\keyword{malaria}
