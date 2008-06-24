\name{qcMean}
\alias{qcMean}

\title{Summarise Relative Allele Scores}
\description{
  The default function for summarising the RAS from the quartets in a probeset into one value per
  probeset per array.
}
\usage{
qcMean(x)
}

\arguments{
  \item{x}{numeric; a vector of RAS from the quartets comprising a single probeset on a single array.}
}

\value{
  Returns the mean of the RAS scores (na.rm=T) unless more than a third of the quartets are NA, in which
  case returns NA.
}

\seealso{
 \code{\link[SNPMaP:cel2ras]{ras2rasS}} for the workflow function that calls the RAS summary function.
 \code{\link[SNPMaP:SNPMaP-package]{SNPMaP-package}}.
}

\keyword{manip}
