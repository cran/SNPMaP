\name{minusMismatch}
\alias{minusMismatch}

\title{Subtract mismatch intensities}
\description{
  Subtract average (within a quartet) mismatch probe (MM) intensities from perfect match probe (PM) intensities
  in a 'short' format \code{SNPMaP} object.
}
\usage{
minusMismatch(x, tidy = FALSE, \dots)
}
\arguments{
  \item{x}{An object of class \code{\link[SNPMaP:SNPMaP-class]{SNPMaP}} in 'short' format.}
  \item{tidy}{logical; should the old \code{FileDoubleMatrix} be \code{\link{unlink}ed}?}
  \item{\dots}{Additional arguments.}
}
\value{
 An object of class \code{\link[SNPMaP:SNPMaP-class]{SNPMaP}}.
}

\seealso{\code{\link[SNPMaP:SNPMaP-package]{SNPMaP-package}}}

\keyword{manip}

