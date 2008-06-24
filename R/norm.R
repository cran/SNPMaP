## norm.R
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

## Natural log of intensities
## Arguments: SNPMaP object  Return:  SNPMaP object

`logIntensities` <-
		function (x, lowMemory = x@lowMemory, tidy = FALSE, ...) 
{
	ca <- deparse(match.call(), width.cutoff = 500)
	da <- date()
	ca <- c(ca, da)
	if (!is(x, "SNPMaP")) 
		stop("this is not a 'SNPMaP' object")
	if (any(x@transformation == 'logInt'))
		stop("these data have already been logged")
	cat("Logging intensities ...\n")
	flush.console()
	x <- .checkOpen(x)
	d <- dim(x@snpdata)
	mat <- matrix()
	if (lowMemory) {
		mat <- .makefile(x@table, d[1], d[2], x@tempDir)
	}else {
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
		mat[, i] <- log(x@snpdata[, i])
		increase(pb)
		flush.console()
	}
	if (is(x@snpdata, "FileDoubleMatrix")) {
		close(x@snpdata)
		if (tidy) 
			unlink(x@snpdata$pathname)
	}
	slot(x, 'snpdata', check=FALSE) <- mat
	if (!all(dim(x@mismatch)==c(0,0))) {
		mm <- matrix()
		if (lowMemory) {
			mm <- .makefile(paste(x@table, "MM", sep = ""), d[1], 
					d[2], x@tempDir)
		}else{
			mm <- .gws(d[1], d[2])
			if (!is.matrix(mm)) {
				warning("out of memory; switching to disk ...")
				mm <- .makefile(paste(x@table, "MM", sep = ""), 
						d[1], d[2], x@tempDir)
				x@lowMemory <- TRUE
			}
		}
		pb2<-ProgressBar(max = 40, stepLength = 40/(dim(x@mismatch)[2]))
		reset(pb2)
		flush.console()
		for (i in 1:dim(x@mismatch)[2]) {
			mm[, i] <- log(x@mismatch[, i])
			increase(pb2)
			flush.console()
		}
		if (is(x@mismatch, "FileDoubleMatrix")) {
			close(x@mismatch)
			if (tidy) 
				unlink(x@mismatch$pathname)
		}
		slot(x, 'mismatch', check=FALSE) <- mm
	}
	x@logInt <- FALSE
	x@majorHistory[[length(x@majorHistory) + 1]] <- ca
	x@transformation <- append(x@transformation, 'logInt')
	close(x)
	return(x)
}

## Quantile normalisation
## Arguments: SNPMaP object of type 'raw', logical, logical
## Return: SNPMaP object
`norm` <-
		function (x, Force = FALSE, lowMemory = x@lowMemory, tidy = FALSE, ...) 
{
	ca <- deparse(match.call(), width.cutoff = 500)
	da <- date()
	ca <- c(ca, da)
	if (!is(x, "SNPMaP")) 
		stop("this is not a 'SNPMaP' object")
	if (x@table != "raw") {
		if (!Force) {
			stop(paste("the data is in", x@table, "form. norm() probably only makes sense for raw data.\n"))
		}else {
			warning(paste("the data is in", x@table, "form. norm() probably only makes sense for raw data.\n"))
		}
	}
	if (any(x@transformation == 'norm'))
		stop("these data have already been normalized")
	n <- dim(x@snpdata)
	cat("Quantile normalizing ...\n")
	flush.console()
	x <- .checkOpen(x)
	cat("Checking memory to make quantile normalization order table ...\n")
	flush.console()
	ords <- .gws(n[1], n[2], "integer")
	gc()
	if (is.matrix(ords)) {
		cat("Memory OK ...\n")
		flush.console()
	}else {
		cat("Making order table on disk ...\n")
		flush.console()
		ords <- .makefile("orders", n[1], n[2], x@tempDir)
	}
	mat <- matrix()
	if (lowMemory) {
		mat <- .makefile(x@table, n[1], n[2], x@tempDir)
	}else {
		mat <- .gws(n[1], n[2])
		gc()
		if (!is.matrix(mat)) {
			warning("out of memory; switching to disk ...")
			mat <- .makefile(x@table, n[1], n[2], x@tempDir)
			x@lowMemory <- TRUE
		}
	}
	toto <- 0
	pb<-ProgressBar(max = 40, stepLength = 40/(n[2]))
	reset(pb)
	flush.console()
	for (i in 1:n[2]) {
		ords[, i] <- rank(x@snpdata[, i])
		toto <- toto + sort(x@snpdata[, i])
		increase(pb)
		flush.console()
	}
	toto <- toto/n[2]
	pb2<-ProgressBar(max = 40, stepLength = 40/(n[2]))
	reset(pb2)
	flush.console()
	for (i in 1:n[2]) {
		mat[, i] <- toto[ords[, i]]
		increase(pb2)
		flush.console()
	}
	if (is(ords, "FileDoubleMatrix")) {
		close(ords)
		unlink(ords$pathname)
	}
	if (is(x@snpdata, "FileDoubleMatrix")) {
		close(x@snpdata)
		if (tidy) 
			unlink(x@snpdata$pathname)
	}
	slot(x, 'snpdata', check=FALSE) <- mat
	x@normalize <- FALSE
	x@majorHistory[[length(x@majorHistory) + 1]] <- ca
	x@transformation <- append(x@transformation, 'norm')
	close(x)
	return(x)
}


