\name{SNPMaP-class}
\docType{class}
\alias{SNPMaP-class}
\alias{SNPMaPdata-class}
\alias{FileDoubleMatrix-class}
\alias{as.matrix,SNPMaP-method}
\alias{boxplot,SNPMaP-method}
\alias{close,SNPMaP-method}
\alias{image,SNPMaP-method}
\alias{open,SNPMaP-method}
\alias{plot,SNPMaP,missing-method}
\alias{summary,SNPMaP-method}
\alias{plot,SNPMaP-method}
\alias{[,SNPMaP,missing,missing,missing-method}
\alias{[,SNPMaP,missing,ANY,missing-method}
\alias{[,SNPMaP,ANY,missing,missing-method}
\alias{[,SNPMaP,ANY,ANY,missing-method}
\alias{print.summary.SNPMaP}
\alias{initialize,SNPMaP-method}

\title{Class SNPMaP}
\description{A class of objects for SNP Microarrays and Pooling using Affymetrix arrays.}

\section{Slots}{
  \describe{
    \item{\code{snpdata}:}{Object of class \code{SNPMaPdata}, either
    a \code{matrix} or a \code{FileDoubleMatrix} with probe intensities or Relative Allele Score
    (RAS) estimates from the SNPMaP experiment.}
    \item{\code{mismatch}:}{As \code{snpdata}. Contains mismatch probe intensities if they are
    present on the microarray and required for analysis.}
    \item{\code{useMM}:}{Object of class \code{logical}, if \code{TRUE} then the mismatch probe intensities
    will be subtracted from the perfect match probe intensities during the analysis. Defaults to \code{FALSE}.}
    \item{\code{normalize}:}{Object of class \code{logical}, if \code{TRUE} then the probe intensities will
    be quantile normalised across arrays during the analysis. Defaults to \code{FALSE}.}
    \item{\code{logInt}:}{Object of class \code{logical}, if \code{TRUE} then the raw probe intensities from the
    array will be replaced with the natural log of the intensities during the analysis. Defaults to \code{FALSE}.}
    \item{\code{summary}:}{Object of class \code{function}, the function that will be used to summarise
    the vector of RAS scores from all the quartets corresponding to a single SNP within an array into a single
    statistic per SNP per array. Should accept a \code{numeric} vector comprising doubles between 0 and 1 and
    \code{NA} and return a single value. Defaults to \code{\link[SNPMaP:qcMean]{qcMean}}.}
    \item{\code{lowMemory}:}{Object of class \code{logical}, if \code{TRUE} then the analysis will
    attempt to minimise the amount of RAM used by storing the \code{snpdata} and \code{mismatch} matricies
    as objects of class \code{FileDoubleMatrix} (see the \code{link[R.huge:R.huge-package]{R.huge}} package)
    on disk at the location named in the \code{tempDir} slot. Defaults to \code{TRUE}.}
    \item{\code{tempDir}:}{Object of class \code{character}. Names the location of a writable
    directory where \code{lowMemory} versions of the \code{snpdata} and \code{mismatch} matricies may be
    stored. Defaults to '.' (ie the working directory).}
    \item{\code{table}:}{Object of class \code{character}. Gives the current format of the \code{snpdata}.
    May be (in order) 'empty', 'raw', 'long', 'short', 'ras' or 'rasS'.}
    \item{\code{chiptype}:}{Object of class \code{character} naming the type of array involved in this
    analysis ('Mapping250K\_Sty', for example). The arrays in a \code{SNPMaP} object must be of a single
    \code{chiptype}.}
    \item{\code{celDim}:}{Object of class \code{numeric} vector of two giving the number of rows
    and columns of probes on the chip. Read from the CEL header by \code{\link[SNPMaP:cel2raw]{cel2raw}()},
    used for pseudo \code{images} of arrays.}
    \item{\code{set}:}{Object of class \code{numeric}. Several of the Affymetrix chips have groups of
    probesets with different numbers of quartets. For example, the Mapping250K\_Sty array has probesets
    comprising 6 or 10 quartets. \code{set} specifies which group of probes are to be included in the
    current \code{SNPMaP} object. It defaults to 1, which corresponds to the largest group of probes on
    the array (ie those with 6 quartets on the Mapping250K\_Sty). Further groups can be accessed by setting
    \code{set} to 2 or higher. For more information, see \code{\link[SNPMaP.cdm:SNPMaP.cdm-package]{SNPMaP.cdm-package}}.}
    \item{\code{snps}:}{Object of class \code{character}. A vector of the SNPs being analysed in the
    current \code{SNPMaP} object.}
    \item{\code{chps}:}{Object of class \code{character}. A vector of the CEL files in the current analysis.
    Defaults to all the files ending '.CEL' in the working directory.}
    \item{\code{cols}:}{Object of class \code{character}. A vector identifying the columns of the \code{snpdata}
    and \code{mismatch} slots of the current \code{SNPMaP} object.}
    \item{\code{width}:}{Object of class \code{numeric}. The number of quartets per SNP in the current group
    of probesets.}
    \item{\code{transformation}:}{Object of class \code{character}. A vector to record transformations carried out
    on the data.}
    \item{\code{experiment}:}{Object of class \code{factor}. A factor identifying the experimental group each array
    belongs to. Not used at present.}
    \item{\code{annotation}:}{Object of class \code{list}. A slot to hold annotation to be used in the
    analysis. Not used at present.}
    \item{\code{created}:}{Object of class \code{character}; the date and time when this \code{SNPMaP} object
    was created.}
    \item{\code{version}:}{Object of class \code{character}; the version number of the \code{SNPMaP} package that created this object.}
    \item{\code{majorHistory}:}{Object of class \code{list}. A record of the major functions carried out
    on this \code{SNPMaP} object, along with the date and time they were called.}
  }
}
\usage{
\S4method{as.matrix}{SNPMaP}(x, mm=FALSE, \dots)
\S4method{plot}{SNPMaP,missing}(x, FUN=function(x){x}, xlim=c(loX, hiX), ylim=c(0, hiY), xlab="guess",
                ylab="Density", main=x@table, col=rainbow(length(x@chps)),
                legend.position="left", legend.bty="n", lty=1:length(x@chps),
                zero.line=TRUE, \dots)
\S4method{boxplot}{SNPMaP}(x, FUN=function(x){x}, ylab="guess", main=x@table, \dots)
\S4method{image}{SNPMaP}(x, chips=x@chps, prompt=FALSE, FUN=log, col=grey(seq(0,1,0.01)), fastRender=4,
				rows=x@celDim[1], cols=x@celDim[2], \dots)
}
\arguments{
  \item{x}{object of class \code{SNPMaP}.}
  \item{mm}{logical indicating whether the data in the \code{mismatch} slot should be returned instead.}
  \item{FUN}{function specifying how each data point should be transformed before plotting. Try \code{log} for probe intensity data.}
  \item{xlim}{the x limits (x1, x2) of the plot. Note that x1 > x2 is allowed and leads to a 'reversed axis'.}
  \item{ylim}{the y limits of the plot.}
  \item{xlab}{character; a label for the x axis, defaults to a description of x.}
  \item{ylab}{character; a label for the y axis, defaults to a description of y.}
  \item{main}{character; a main title for the plot.}
  \item{col}{the colours for the lines.}
  \item{legend.position}{position of the legend in the pane to the right of the graph.}
  \item{legend.bty}{the type of box to be drawn around the legend. The allowed values are "o" and "n".}
  \item{lty}{the line type.}
  \item{zero.line}{logical; if TRUE, add a base line at y = 0}
  \item{chips}{character or numeric; the arrays to be imaged from the \code{chps} slot.}
  \item{prompt}{logical indicating whether the user should be asked before moving on to the next image.}
  \item{fastRender}{numeric \emph{n} indicating 1/\emph{n} rows and columns should be \code{imaged}.}
  \item{rows}{numeric; the number of rows of probes on the array.}
  \item{cols}{numeric; the number of columns of probes on the array.}
  \item{\dots}{additional arguments passed to methods.}
}
\section{Methods}{
  \describe{
    \item{as.matrix}{\code{signature(x="SNPMaP")}: returns the contents of the \code{snpdata} slot as a \code{matrix}, or, if \code{mm=TRUE}, returns the contents of the \code{mismatch} slot.}
    \item{plot}{\code{signature(x="SNPMaP", y="missing")}: plots the \code{density} of the contents of the \code{snpdata} slot, with each array represented by a single line identified in a legend. Uses \code{layout}.}
    \item{boxplot}{\code{signature(x="SNPMaP")}: box plot for each array in the \code{snpdata} slot.}
    \item{image}{\code{signature(x="SNPMaP")}: if the object is in 'raw' format, draws a pseudo image of the array as it was scanned to check for artefacts.}
    \item{open}{\code{signature(con="SNPMaP")}: if the data in the \code{SNPMaP} object is stored on disk, open a connection.}
    \item{close}{\code{signature(con="SNPMaP")}: if the data is on disk, close the connection.}
    \item{summary}{\code{signature(object="SNPMaP")}: summary of the \code{SNPMaP} object.}
    \item{[}{\code{signature(x="SNPMaP", i="missing", j="missing")}: Extract subsets of the data contained in the \code{snpdata} slot of the \code{SNPMaP} object. Using a single index (\code{[]} rather than \code{[,]}) returns \emph{rows}.}
    \item{[}{\code{signature(x="SNPMaP", i="ANY", j="missing")}: see above.}
    \item{[}{\code{signature(x="SNPMaP", i="missing", j="ANY")}: see above.}
    \item{[}{\code{signature(x="SNPMaP", i="ANY", j="ANY")}: see above.}
    \item{initialize}{\code{signature(.Object="SNPMaP")}: default initialize method.}
  }
}

\seealso{
  \code{\link[SNPMaP:SNPMaP-package]{SNPMaP-package}}.\cr
  \code{\link[SNPMaP.cdm:SNPMaP.cdm-package]{SNPMaP.cdm-package}}.\cr
  \code{\link[SNPMaP:snpmap]{snpmap}()} to set up a SNPMaP analysis.\cr
  \code{\link[SNPMaP:disk2memory]{disk2memory}()} to transfer \code{SNPMaP} objects between disk and memory.\cr
  \code{\link[SNPMaP:cel2ras]{cel2ras}()} for workflow functions.\cr
  \code{\link[SNPMaP:cloneSNPMaP]{cloneSNPMaP}()} to copy a \code{SNPMaP} object on disk.\cr
  \code{\link[SNPMaP:writeSNPMaP]{writeSNPMaP}()} to write \code{SNPMaP} objects to text files.
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
\keyword{classes}


