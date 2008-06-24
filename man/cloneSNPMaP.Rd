\name{cloneSNPMaP}
\alias{cloneSNPMaP}

\title{Copy SNPMaP objects on disk}
\description{
  If you try to copy a SNPMaP object on disk to another object, it will just copy the \code{FileDoubleMatrix}
  filehandles. To copy the data, you need to \code{cloneSnpmap()}.
}
\usage{
cloneSNPMaP(x, ...)
}
\arguments{
  \item{x}{An object of class \code{SNPMaP}.}
  \item{\dots}{Additional arguments.}
}

\value{
  An object of class \code{SNPMaP}: a copy of \code{x}.
}

\seealso{\code{\link[SNPMaP:SNPMaP-class]{SNPMaP-class}}}

\keyword{manip}
