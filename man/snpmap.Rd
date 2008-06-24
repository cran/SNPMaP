\name{snpmap}
\alias{snpmap}
\alias{msnpmap}
\title{SNP Microarrays and Pooling}
\description{
  Functions to process SNPMaP data from CEL files to Relative Allele Score summaries (rasS).
}
\usage{
snpmap(cels = dir(pattern = ".[cC][eE][lL]$"), lowMemory = TRUE, set = 1,
        useMM = FALSE, normalize = FALSE, log.intensities = FALSE, ras.summary = qcMean,
        tempDir = ".", RUN = "cel2rasS", interactive=FALSE, ...)
msnpmap(cels = dir(pattern = ".[cC][eE][lL]$"), lowMemory = TRUE, set = 0,
        useMM = FALSE, normalize = FALSE, log.intensities = FALSE, ras.summary = qcMean,
        tempDir = ".", RUN = "cel2rasS", interactive=FALSE, ...)
}
\arguments{
  \item{cels}{character; a vector of CEL files to use in the analysis}
  \item{lowMemory}{logical; should the SNP data be stored on disk rather than in memory?}
  \item{set}{numeric, a single number for \code{snpmap} or a vector for \code{msnpmap}; the probesets to include see \code{\link[SNPMaP.cdm:SNPMaP.cdm-package]{SNPMaP.cdm-package}}. The special \code{set=0} corresponds to all the \code{sets} on the array.}
  \item{useMM}{logical; should mismatch probes be used in the analysis (if they are present on the array)?}
  \item{normalize}{logical; should the probe intensities be quantile normalised across arrays?}
  \item{log.intensities}{logical; use the natural log of the probe intensities?}
  \item{ras.summary}{function to compute summary statistic for RAS scores on each chip. See \code{\link[SNPMaP:SNPMaP-class]{SNPMaP-class}}.}
  \item{tempDir}{character; writable directory to store the probe data if \code{lowMemory=TRUE}.}
  \item{RUN}{character, name of a workflow function \code{cel2*}; see \code{\link[SNPMaP:cel2ras]{cel2ras}()}.}
  \item{interactive}{logical; should \code{snpmap()} prompt for CEL files (Windows only)?}
  \item{\dots}{addtional arguments passed on to workflow function.}
}
\details{
  \code{snpmap()} sets up the SNPMaP analysis and calls a workflow function to extract the data from the
  CEL files to 'raw' intensity data for all probes on the array (at any \code{set}), to a 'long' format
  matrix (one column per array), to a 'short' format matrix (one row per probeset), to 'ras' (a matrix of
  Relative Allele Scores, one per quartet, one row per probeset), or to 'rasS' (a matrix of RAS, one per
  probeset, one column per array). Specifying \code{lowMemory=TRUE} stores the data as a
  \code{\link[R.huge:FileMatrix]{FileDoubleMatrix}} on disk rather than in memory.
}
\value{
  Object of class \code{\link[SNPMaP:SNPMaP-class]{SNPMaP}} for \code{snpmap()}; list of \code{SNPMaP}
  objects for \code{msnpmap()}.
}
\section{Warning}{
  \code{snpmap()} returns the SNPs from a single \code{\link[SNPMaP.cdm:SNPMaP.cdm-package]{set}} on the array. For
  some arrays you will need to call \code{snpmap()} several times to access all the SNPs, or \code{msnpmap(set=0)}.
}
\seealso{
 \code{\link[SNPMaP:SNPMaP-package]{SNPMaP-package}}.\cr
 \code{\link[SNPMaP.cdm:SNPMaP.cdm-package]{SNPMaP.cdm-package}}.\cr
 \code{\link[SNPMaP:SNPMaP-class]{SNPMaP-class}} for SNPMaP objects and methods.\cr
 \code{\link[SNPMaP:cel2ras]{cel2ras}()} for workflow functions.\cr
 \code{\link[SNPMaP:disk2memory]{disk2memory}()} to transfer \code{SNPMaP} objects between disk and memory.\cr
 \code{\link[SNPMaP:cloneSNPMaP]{cloneSNPMaP}()} to copy a \code{SNPMaP} object on disk.\cr
 \code{\link[SNPMaP:writeSNPMaP]{writeSNPMaP}()} to write \code{SNPMaP} objects to text files.\cr
 \code{\link[SNPMaP:norm]{norm}()}, the quantile normalization function.\cr
 \code{\link[SNPMaP:logIntensities]{logIntensities}()}, the function that logs the raw intensities.
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
