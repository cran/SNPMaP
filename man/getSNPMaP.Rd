\name{getSNPMaP}
\alias{getSNPMaP}
\alias{setSNPMaP}

\title{Get or set slots of a SNPMaP object}
\description{
  Get or set the values of the slots in a \code{SNPMaP} object.
}
\usage{
getSNPMaP(x, \dots)
setSNPMaP(x, namedList)
}
\arguments{
  \item{x}{An object of class \code{SNPMaP}.}
  \item{namedList}{A list of the form \code{list(useMM=TRUE, normalize=FALSE)}; the exact name of a slot followed by the new value.}
  \item{\dots}{Slots to retrieve data from (may be partial matches).}
}
\details{
  Currently, \code{getSNPMaP()} will retrieve the values of the following slots: \code{useMM, normalize, logInt, summary,
  lowMemory, tempDir, table, chiptype, celDim, set, snps, chps, cols, width, transformation, experiment,
  created, version, majorHistory}. \code{setSNPMaP()} will set the following slots: \code{useMM, normalize, logInt, summary,
  lowMemory, tempDir, set, experiment}. Comments can be retrived and set using \code{comment()}.
}

\value{
  For \code{getSNPMaP()}, a named list of slot values (suitable as an argument for \code{setSNPMaP()}).
  For \code{setSNPMaP()}, a named list of the old slot values (suitable as an argument for \code{setSNPMaP()}).
}

\section{Warning}{
  \code{setSNPMaP()} does not clone the data in \code{SNPMaP} objects stored on disk (see \code{\link[SNPMaP:cloneSNPMaP]{cloneSNPMaP}()}.)
}

\seealso{\code{\link[SNPMaP:SNPMaP-class]{SNPMaP-class}}}

\examples{
\dontrun{
 ## Retrieve slot values
 getSNPMaP(x, 'experiment', 'snps')
 ## Set new slot values
 x <- setSNPMaP(x, list(useMM=TRUE, normalize=FALSE))
}
}

\keyword{manip}
