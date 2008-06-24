\name{writeSNPMaP}
\alias{writeSNPMaP}
\alias{readSNPMaP}

\title{Write a SNPMaP object to a text file}
\description{
  Transfer the data in a SNPMaP object to a text file (and read the file back into \R efficiently).
}
\usage{
writeSNPMaP(x, file="", mm=FALSE, sep="\t", dec=".", transpose=TRUE, \dots)
readSNPMaP(file="", sep="\t", dec=".", transpose=TRUE, \dots)
}
\arguments{
  \item{x}{object of class \code{SNPMaP}.}
  \item{file}{character; name of the file that should be created or appended to. The default, "", causes the output to be written to the console.}
  \item{mm}{logical indicating whether the data in the \code{mismatch} slot should be returned instead.}
  \item{sep}{character; the separator to be used between columns in the file.}
  \item{dec}{character; the decimal point to be used.}
  \item{transpose}{logical; should the data frame be transposed for writing to file (or transposed when reading from file)?}
  \item{\dots}{additional arguments passed to methods.}
}
\details{
 \code{transpose} defaults to \code{TRUE} because it is \emph{much} faster to write one
 (or a few) row(s) per chip than it is to write one row per SNP. It also produces smaller text files.
 Set \code{transpose=FALSE} to preserve the orientation of the matrix. Likewise, \code{readSNPMaP} assumes
 the file was written using the default \code{transpose=TRUE}; change this to \code{FALSE} to read in an
 untransposed file.
}
\value{
 \code{writeSNPMaP} returns \code{TRUE} invisibly on success. \code{readSNPMaP} returns a named matrix.
}

\seealso{\code{\link[SNPMaP:SNPMaP-package]{SNPMaP-package}}}
\examples{
\dontrun{
 ## Write to a tab-delimited file
 writeSNPMaP(x, file='mySNPMaP.dat')
 ## Write to a csv file
 writeSNPMaP(x, file='mySNPMaP.csv', sep=',')
 ## Write a semicolon delimited file with commas for decimal points
 ## (Standard in some parts of Western Europe)
 writeSNPMaP(x, file='mySNPMaP.csv', sep=';', dec=',')
 ## Read a tab-delimited file written using writeSNPMaP()
 y<-readSNPMaP(file='mySNPMaP.dat')
}
}
\keyword{manip}

