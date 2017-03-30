\name{gonHosh}

\alias{gonHosh}

\title{
Calculate biting rate for Anopheles gambiae s.l.
}

\description{
Calculate biting rate for Anopheles gambiae s.l. based on the number of days it takes to complete a gonothropic cycle
}

\usage{
gonHosh(T)
}

\arguments{
  \item{T}{Air temperature in C}
}

\details{

\code{gonHosh} calculate the biting rate for Anopheles gambiae s.l.
}

\value{
A vector with biting rate per day
}

\references{
Hoshen, M.B., and Morse, A.P.: A model structure for estimating malaria risk, Environmental Change and Malaria Risk: Global and Local Implications, volume 9, 10, Eds: Takken, Willem, Martens, Pim, and Bogers, Robert J., 2005
}

\author{
Torleif Markussen Lunde (Centre for International Health/Geophysical Institute, University of Bergen, Norway)
}

\note{}

\seealso{
}

\examples{
bite <- gonHosh(20)
}


\keyword{malaria}
