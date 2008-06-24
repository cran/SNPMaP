\name{disk2memory}
\alias{disk2memory}
\alias{memory2disk}

\title{Transferring SNPMaP objects between disk and memory}
\description{
  Transfer the data in a \code{SNPMaP} object from \code{FileDoubleMatrix} to memory or from memory
  to \code{FileDoubleMatrix}.
}
\usage{
disk2memory(x, tidy = FALSE, ...)
memory2disk(x, ...)
}
\arguments{
  \item{x}{An object of class \code{\link[SNPMaP:SNPMaP-class]{SNPMaP}}.}
  \item{tidy}{logical; should the old \code{FileDoubleMatrix} be \code{\link{unlink}ed}?}
  \item{\dots}{Additional arguments.}
}
\value{
 An object of class \code{\link[SNPMaP:SNPMaP-class]{SNPMaP}}.
}

\seealso{\code{\link[SNPMaP:SNPMaP-package]{SNPMaP-package}}}

\keyword{manip}

