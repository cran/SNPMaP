\name{cel2ras}
\alias{cel2ras}
\alias{cel2rasS}
\alias{cel2raw}
\alias{cel2long}
\alias{cel2short}
\alias{raw2long}
\alias{long2short}
\alias{short2ras}
\alias{ras2rasS}
\alias{raw2rasS}
\alias{raw2short}
\alias{raw2ras}
\alias{long2ras}
\alias{long2rasS}
\alias{short2rasS}

\title{SNPMaP workflows}
\description{
  Workflow functions to process \code{SNPMaP} objects. Called by \code{\link[SNPMaP:snpmap]{snpmap}()}.
}
\usage{
cel2raw(x, cels = x@chps, lowMemory = x@lowMemory, ...)
cel2long(x, cels = x@chps, lowMemory = x@lowMemory, set = x@set, normalize = x@normalize,
        log.intensities = x@logInt, ...)
cel2short(x, cels = x@chps, lowMemory = x@lowMemory, set = x@set, normalize = x@normalize,
        log.intensities = x@logInt, ...)
cel2ras(x, cels = x@chps, lowMemory = x@lowMemory, set = x@set, normalize = x@normalize,
        log.intensities = x@logInt, subtractMismatch = x@useMM, ...)
cel2rasS(x, cels = x@chps, lowMemory = x@lowMemory, set = x@set, normalize = x@normalize,
        log.intensities = x@logInt, subtractMismatch = x@useMM, FUN = x@summary, ...)
raw2long(x, lowMemory = x@lowMemory, set = x@set, normalize = x@normalize,
        log.intensities = x@logInt, tidy = FALSE, ...)
raw2short(x, lowMemory = x@lowMemory, set = x@set, normalize = x@normalize,
        log.intensities = x@logInt, tidy = FALSE, ...)
raw2ras(x, lowMemory = x@lowMemory, set = x@set, normalize = x@normalize,
        log.intensities = x@logInt, subtractMismatch = x@useMM, tidy = FALSE, ...)
raw2rasS(x, lowMemory = x@lowMemory, set = x@set, normalize = x@normalize,
        log.intensities = x@logInt, subtractMismatch = x@useMM, FUN = x@summary,
        tidy = FALSE, ...)
long2short(x, lowMemory = x@lowMemory, tidy = FALSE, ...)
long2ras(x, lowMemory = x@lowMemory, subtractMismatch = x@useMM, tidy = FALSE, ...)
long2rasS(x, lowMemory = x@lowMemory, subtractMismatch = x@useMM, FUN = x@summary,
		tidy = FALSE, ...)
short2ras(x, lowMemory = x@lowMemory, subtractMismatch = x@useMM, tidy = FALSE, ...)
short2rasS(x, lowMemory = x@lowMemory, subtractMismatch = x@useMM, FUN = x@summary,
		tidy = FALSE, ...)
ras2rasS(x, lowMemory = x@lowMemory, FUN = x@summary, tidy = FALSE, ...)
}
\arguments{
  \item{x}{An object of class \code{\link[SNPMaP:SNPMaP-class]{SNPMaP}}.}
  \item{cels}{character; a vector of the CEL files to be included.}
  \item{lowMemory}{logical; should the SNP data be stored on disk rather than in memory?}
  \item{tidy}{logical; if the data is stored on disk, should the old \code{FileDoubleMatrix} be \code{\link[base:unlink]{unlink}ed}?}
  \item{set}{numeric; the set of probesets to include. See \code{\link[SNPMaP.cdm:SNPMaP.cdm-package]{SNPMaP.cdm-package}}.}
  \item{normalize}{logical; quantile normalize probe intensities across arrays.}
  \item{log.intensities}{logical; take natural log of probe intensities.}
  \item{subtractMismatch}{logical; subtract mismatch probe intensities if available from perfect match intensities.}
  \item{FUN}{function; function to compute summary statistic for RAS on each chip. See \code{\link[SNPMaP:SNPMaP-class]{SNPMaP-class}}.}
  \item{\dots}{additional arguments passed to other workflow and internal functions.}
}
\details{
  The workflow functions work by calling \code{cel2raw()}, \code{raw2long()}, \code{long2short()},
  \code{short2ras()}, \code{ras2rasS()}.
}
\value{
  Object of class \code{\link[SNPMaP:SNPMaP-class]{SNPMaP}}.
}
\seealso{
 \code{\link[SNPMaP:SNPMaP-package]{SNPMaP-package}},
 \code{\link[SNPMaP.cdm:SNPMaP.cdm-package]{SNPMaP.cdm-package}},
 \code{\link[SNPMaP:SNPMaP-class]{SNPMaP-class}},
 \code{\link[SNPMaP:snpmap]{snpmap}()},
 \code{\link[SNPMaP:norm]{norm}()},
 \code{\link[SNPMaP:logIntensities]{logIntensities}()}
}

\examples{
\dontrun{
 ## Getting started
 ## Creates the 'raw' SNPMaP object x on disk with mismatch probes included
 x<-snpmap(useMM=TRUE, RUN='cel2raw', lowMemory=TRUE)
 ## Print a summary of the SNPMaP object
 summary(x)
 ## Add a comment (prints in the summary)
 comment(x)<-'High and low extreme pools from January'
 ## View pseudo image to screen for artefacts
 image(x)
 ## Plot probe intensities
 plot(x, FUN=log)
 boxplot(x, FUN=log)
 ## tidy=TRUE removes the FileDoubleMatrix from the old x to keep the disk tidy
 x<-raw2ras(x, tidy=TRUE)
 ## Plot Relative Allele Scores
 plot(x)
 ## Default tidy=FALSE does not remove the original FileDoubleMatrix from disk
 ## Useful if you want to keep x (no side effects)
 y<-ras2rasS(x)
 ## View the first ten rows
 as.matrix(y[1:10,])
 ## View a set of SNPs
 as.matrix(y[c("SNP_A-4192909", "SNP_A-4192918"),])
 ## Transfer the SNPMaP object from disk to memory
 y<-disk2memory(y, tidy=TRUE)
 ## Run the analysis again from CEL files to RAS summaries without viewing intermediate stages
 ## This time in memory (may require a lot of RAM)
 z<-snpmap(useMM=TRUE, RUN='cel2rasS', lowMemory=FALSE)
 plot(z)
 ## Get the RAS summary scores as a standard matrix
 rasSummaries<-as.matrix(z)
 ## Read all the sets into a list
 allSets<-msnpmap(set=0)
}
}
\keyword{manip}

