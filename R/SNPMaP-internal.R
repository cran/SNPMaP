## SNPMaP-internal.R
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

## Internal function.  Arguments: SNPMaP object  Return: SNPMaP object with open 
## connection to file-mapped snpdata if present

`.checkOpen` <-
		function (x) 
{
	if (is(x@snpdata, "FileDoubleMatrix")) {
		if (!isOpen(x@snpdata)) {
			open(x@snpdata)
			isOpen(x@snpdata) || stop("cannot open connection to ", 
					x@snpdata$pathname)
		}
		if (is(x@mismatch, "FileDoubleMatrix")) {
			if (!isOpen(x@mismatch)) {
				open(x@mismatch)
				isOpen(x@mismatch) || stop("cannot open connection to ", 
						x@mismatch$pathname)
			}
		}
	}
	return(x)
}

## Internal function.  Arguments: length and mode of desired object
## Return: desired vector or FALSE if there is not enough memory available
## used for 'great white shark' method of memory checking
## usage:
## d   <- dim(x@snpdata)
## mat <- .gws(d[1], d[2])
## if (length(mat) == 1) {
##    warning("out of memory; switching to disk ...")
##    mat <- .makefile(x@table, d[1], d[2], x@tempDir)
## }

`.gws` <-
		function (rows, cols, mode = "numeric") 
{
	gc()
	tryCatch({
				fill <- vector(mode, 1)
				res <- matrix(fill, nrow=rows, ncol=cols)
				invisible(res)
			}, error = function(e) {
				warning("out of memory")
				return(FALSE)
			})
}


## Internal function. Test if cel file exists and if name has .cel ext 
## isCelFile() tweaked to close connections on exit
## Arguments: character  Return:  logical

`.isCel` <-
		function (filename, ...) 
{
	if (!file.exists(filename)) {
		stop("Cannot check file format. File not found: ", filename)
	}
	con <- NULL
	on.exit({
				if (inherits(con, "connection") && isOpen(con)) close(con)
			})
	header <- NULL
	for (ver in c("4", "3", "1")) {
		tryCatch({
					if (inherits(con, "connection") && isOpen(con)) {
						close(con)
						con <- NULL
					}
					con <- file(filename, open = "rb")
					if (ver == "4") {
						header <- affxparser:::.readCelHeaderV4(con)
					}
					else if (ver == "3") {
						header <- affxparser:::.readCelHeaderV3(con)
					}
					else {
						header <- readCcgHeader(con)
						dataTypeId <- header$dataHeader$dataTypeId
						if (!identical(dataTypeId, "affymetrix-calvin-intensity")) 
							header <- NULL
					}
				}, error = function(ex) {
				})
	}
	isCelFile <- (!is.null(header))
	isCelFile
}

#{
#	if (!file.exists(filename)) {
#		stop("cannot find ", filename)
#	}
#	isCel <- logical(0)
#	if (length(grep(".[cC][eE][lL]$", filename)) < 1) {
#		isCel <- FALSE
#	}
#	else {
#		isCel <- TRUE
#	}
#	return(isCel)
#}


## Internal function to create memory-mapped file, wrapper for FileDoubleMatrix
## Arguments: filename prefix, number of rows, number of columns, path
## Return: a connection

`.makefile` <-
		function (prefix, rowN, colN, location = ".") 
{
	tryCatch({
				path <- paste(location, "temp_SNPMaP", sep = "/")
				if (!any(regexpr("^temp_SNPMaP$", dir(location)) > 0)) {
					mkdirs(path) || stop("cannot create a temporary directory '", 
							path, "\n")
				}
				else if (file.info(path)[2] == FALSE) {
					mkdirs(path) || stop("cannot create a temporary directory '", 
							path, "'\n")
				}
				fileName <- character(0)
				fileName <- tempfile(prefix, path)
				mat <- FileDoubleMatrix(fileName, nrow = rowN, ncol = colN)
				return(mat)
			}, error = function(e) {
				print(e)
				warning("there was a problem creating the temporary file '", 
						fileName, "'; attempting to finalize.")
				if (exists("mat")) {
					if (isOpen(mat)) {
						close(mat)
					}
					rm(mat)
					unlink(fileName)
				}
			})
}


## Internal function. Average mismatch intensities across A and B probes
## Arguments: SNPMaP object with mismatch data, logical
## Return: SNPMaP object

`.meanMismatch` <-
		function (x, tidy = FALSE, ...) 
{
	ca <- deparse(match.call())
	da <- date()
	ca <- c(ca, da)
	if (!is(x, "SNPMaP")) 
		stop("this is not a 'SNPMaP' object")
	x <- .checkOpen(x)
	if (all(dim(x@mismatch)==c(0,0))) 
		stop("cannot find mismatch data")
	if (x@table != "short") 
		stop("table must be in 'short' format")
	cat("Calculating mean mismatch intensities for A and B alleles ...\n")
	flush.console()
	w <- x@width * 2
	d <- dim(x@mismatch)
	mm <- matrix()
	if (x@lowMemory) {
		mm <- .makefile(paste(x@table, "MM", sep = ""), d[1], 
				d[2], x@tempDir)
	}
	else {
		mm <- .gws(d[1], d[2])
		if (!is.matrix(mm)) {
			warning("out of memory; switching to disk ...")
			mm <- .makefile(paste(x@table, "MM", sep = ""), d[1], 
					d[2], x@tempDir)
			x@lowMemory <- TRUE
		}
	}
	pb<-ProgressBar(max = 40, stepLength = 40/length(x@chps))
	reset(pb)
	flush.console()
	for (i in 1:length(x@chps)) {
		mm[, (1 + ((i - 1) * w)):(i * w)] <- .mm(x@mismatch[, (1 + ((i - 1) * w)):(i * w)])
		increase(pb)
		flush.console()
	}
	if (x@lowMemory & is(x@mismatch, "FileDoubleMatrix")) {
		close(x@mismatch)
		if (tidy) 
			unlink(x@mismatch$pathname)
	}
	slot(x, 'mismatch', check=FALSE) <- mm
	x@majorHistory[[length(x@majorHistory) + 1]] <- ca
	close(x)
	return(x)
}

## Internal function. Mean mismatch intensities across A and B alleles.
## Arguments: mismatch matrix corresponding to one array
## Return: mismatch matrix
`.mm` <-
		function (x) 
{
	v <- dim(x)[2]/2
	rep((x[, 1:v] + x[, (v + 1):(2 * v)])/2, 2)
}


## Internal function. Calculate RAS from A and B allele intensities
## Arguments: matrix corresponding to one array
## Return: matrix
`.ra` <-
		function (x) 
{
	v <- dim(x)[2]/2
	x[, 1:v]/(x[, 1:v] + x[, (v + 1):(2 * v)])
}


