\name{pfF}

\alias{pfF}

\title{
Calculate EIP for Plasmodium falciparum
}

\description{
Time from gametocyte ingestion by a mosquito to release of sporozoites from
the salivary glands (extrinsic incubation period, EIP)
}

\usage{
pfF(T)
}

\arguments{
  \item{T}{Air temperature in C}
}

\details{

\code{pfF} returns EIP in days.
}

\value{
A vector with EIP
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
EIP <- pfF(20)
}


\keyword{malaria}
