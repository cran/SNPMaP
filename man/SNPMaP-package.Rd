\name{SNPMaP-package}
\alias{SNPMaP-package}
\docType{package}
\title{
SNP Microarrays and Pooling in \R
}
\description{
 Pooling DNA on SNP microarrays is a cost-effective way to carry out genome-wide association studies
 for heritable disorders or traits. The SNPMaP package provides formal SNPMaP objects and methods in \R as
 a base for these analyses using Affymetrix genotyping arrays.
}
\details{
\tabular{ll}{
Package: \tab SNPMaP\cr
Type: \tab Package\cr
Version: \tab 1.0.2\cr
Date: \tab 2008-06-24\cr
License: \tab GPL (>= 3)\cr
}
 A SNPMaP study involves selecting individuals based on their disease status (case or control) or their
 score on a quantitative trait (high or low extremes). DNA from these individuals is pooled in biological
 replicate pools (20 pools of 50 different individuals each, for example). Each pool is then genotyped
 according to the standard protocols on genotyping microarrays. The probe intensities from these arrays
 can be extracted from the CEL files into \code{SNPMaP} objects using the \code{snpmap()} function and
 Relative Allele Scores (RAS) can be calculated as estimates of allele frequency in each pool. By comparing
 the allele frequencies in high and low or case and control pools, hundreds of thousands of SNPs across
 the genome can be screened for association with a trait or disorder. The \code{SNPMaP} package also
 provides methods for visualising the data at each stage of the analysis. Using the \code{lowMemory}
 option allows this to be done on a standard desktop computer (albeit slower than if all data is kept
 in memory). This package is an evolution of the scripts referred to in Meaburn et al (2006) and is described
 in Davis, Plomin and Schalkwyk (submitted for publication); please cite this paper if you find the package useful.
 Additional supporting material is available at \url{http://sgdp.iop.kcl.ac.uk/snpmap/}.
}
\section{Array types}{
Details of the arrays supported by the current version of \code{SNPMaP} can be found in the \code{\link[SNPMaP.cdm:SNPMaP.cdm-package]{SNPMaP.cdm}} package.
}
\section{Future work}{
 Future releases of the package will build on the range of chips covered as well as introducing pre-packaged
 methods for carrying out the genome-wide screen for association.
}
\author{
Oliver SP Davis and Leo C Schalkwyk

Maintainer: Oliver SP Davis \email{snpmap@iop.kcl.ac.uk}
}
\references{
Davis, OSP, Plomin, R, and Schalkwyk, LC. (submitted for publication) The SNPMaP package for R:
A framework for genome-wide association using DNA pooling on microarrays.\cr
\cr
Meaburn E, Butcher LM, Schalkwyk LC, and Plomin R. (2006) Genotyping pooled DNA using 100K SNP microarrays:
a step towards genomewide association scans. \emph{Nucleic Acids Research, 34(4):e28}. \url{http://dx.doi.org/10.1093/nar/gnj027}
}
\keyword{package}
\seealso{
 \code{\link[SNPMaP:snpmap]{snpmap}()} to set up a SNPMaP analysis.\cr
 \code{\link[SNPMaP:SNPMaP-class]{SNPMaP-class}} to represent a SNPMaP study.\cr
 \code{\link[SNPMaP.cdm:SNPMaP.cdm-package]{SNPMaP.cdm-package}} for the \code{cdm} matrices that interpret the 'raw' format \code{SNPMaP} objects.\cr
 \code{\link[affxparser:affxparser-package]{affxparser-package}} that reads the CEL files.\cr
 \code{\link[R.huge:R.huge-package]{R.huge-package}} that provides \code{FileDoubleMatrices} for the \code{lowMemory} option.\cr
 \code{\link[methods:methods-package]{methods-package}} for S4 formal classes.
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

