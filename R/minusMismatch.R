## minusMismatch.R
## Copyright (c) 2008 Oliver SP Davis and Leo Schalkwyk

## This file is part of the SNPMaP package.

## SNPMaP is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.

## SNPMaP is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.

## You should have received a copy of the GNU General Public License
## along with SNPMaP.  If not, see <http://www.gnu.org/licenses/>.

## Subtract mismatch intensities from perfect match intensities
## Arguments: SNPMaP object, logical. Return: SNPMaP object
`minusMismatch` <-
		function (x, tidy = FALSE, ...) 
{
	ca <- deparse(match.call(), width.cutoff = 500)
	da <- date()
	ca <- c(ca, da)
	if (!is(x, "SNPMaP")) 
		stop("this is not a 'SNPMaP' object")
	if (any(x@transformation == 'subtractMismatch'))
		stop("these data have already been mismatch-subtracted")
	if (all(dim(x@mismatch)==c(0,0))) 
		stop("cannot find mismatch data")
	if (!all(dim(x@mismatch) == dim(x@snpdata))) 
		stop("dimensions of mismatch data do not match dimensions of perfect match data")
	cat("Subtracting mean mismatch intensities ...\n")
	flush.console()
	x <- .meanMismatch(x, tidy = tidy)
	x <- .checkOpen(x)
	d <- dim(x@snpdata)
	mat <- matrix()
	if (x@lowMemory) {
		mat <- .makefile(x@table, d[1], d[2], x@tempDir)
	}
	else {
		mat <- .gws(d[1], d[2])
		if (!is.matrix(mat)) {
			warning("out of memory; switching to disk ...")
			mat <- .makefile(x@table, d[1], d[2], x@tempDir)
			x@lowMemory <- TRUE
		}
	}
	pb<-ProgressBar(max = 40, stepLength = 40/(dim(x@snpdata)[2]))
	reset(pb)
	flush.console()
	for (i in 1:dim(x@snpdata)[2]) {
		y <- x@snpdata[, i] - x@mismatch[, i]
		y[y < 0] <- 0
		mat[, i] <- y
		increase(pb)
		flush.console()
	}
	if (x@lowMemory) {
		if (is(x@snpdata, "FileDoubleMatrix")) {
			close(x@snpdata)
			if (tidy) 
				unlink(x@snpdata$pathname)
		}
		if (is(x@mismatch, "FileDoubleMatrix")) {
			close(x@mismatch)
			unlink(x@mismatch$pathname)
		}
	}
	slot(x, 'snpdata', check=FALSE) <- mat
	slot(x, 'mismatch', check=FALSE) <- new('matrix')
	x@useMM <- FALSE
	x@majorHistory[[length(x@majorHistory) + 1]] <- ca
	x@transformation <- append(x@transformation, 'subtractMismatch')
	close(x)
	return(x)
}

