\name{norm}
\alias{norm}
\alias{logIntensities}

\title{Transform a SNPMaP object}
\description{
  Quantile normalize or (natural) log the intensities in a 'raw' format \code{SNPMaP} object.
}
\usage{
norm(x, Force = FALSE, lowMemory = x@lowMemory, tidy = FALSE, \dots)
logIntensities(x, lowMemory = x@lowMemory, tidy = FALSE, \dots)
}
\arguments{
  \item{x}{An object of class \code{\link[SNPMaP:SNPMaP-class]{SNPMaP}} in 'raw' format.}
  \item{Force}{logical; \code{Force=TRUE} will allow \code{SNPMaP} objects in formats other than 'raw' to be quantile normalized.}
  \item{lowMemory}{logical; should the data be stored on disk rather than in memory?}
  \item{tidy}{logical; should the old \code{FileDoubleMatrix} be \code{\link{unlink}ed}?}
  \item{\dots}{Additional arguments.}
}
\value{
 An object of class \code{\link[SNPMaP:SNPMaP-class]{SNPMaP}}.
}

\seealso{\code{\link[SNPMaP:SNPMaP-package]{SNPMaP-package}}}

\keyword{manip}

