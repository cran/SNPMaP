## cel2ras.R
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

## Workflow function wrapper. Converts SNPMaP objects from CEL files to 'long' format
## Arguments: SNPMaP object, [character, numeric, logical, logical, ...]
## Returns: SNPMaP object
`cel2long` <- 
		function (x, cels = x@chps, lowMemory = x@lowMemory, set = x@set, normalize = x@normalize,
				log.intensities = x@logInt, ...)
{
	## Catch call for majorHistory slot, add datetime
	ca<-deparse(match.call(), width.cutoff=500)
	d<-date()
	ca<-c(ca, d)
	if(!is(x, "SNPMaP")) stop(deparse(substitute(x)), " is not a 'SNPMaP' object")
	x@majorHistory[[length(x@majorHistory)+1]]<-ca
	## Use workflow functions in order
	tryCatch({
				x<-cel2raw(x, cels=cels, lowMemory = lowMemory, ...)
				x<-raw2long(x, set = set, normalize = normalize,
						log.intensities = log.intensities, tidy=TRUE, ...)
				return(x)
			}, error=function(e){
				print(e)
			}, interrupt=function(i){
				print(i)
			})
}

## Workflow function wrapper. Converts SNPMaP objects from CEL files to 'ras' format
## Arguments: SNPMaP object, [character, numeric, logical, logical, logical, ...]
## Returns: SNPMaP object
`cel2ras` <-
		function (x, cels = x@chps, lowMemory = x@lowMemory, set = x@set, normalize = x@normalize,
				log.intensities = x@logInt, subtractMismatch = x@useMM, ...) 
{
	## Catch call for majorHistory slot, add datetime
	ca <- deparse(match.call(), width.cutoff = 500)
	d <- date()
	ca <- c(ca, d)
	if (!is(x, "SNPMaP")) 
		stop(deparse(substitute(x)), " is not a 'SNPMaP' object")
	x@majorHistory[[length(x@majorHistory) + 1]] <- ca
	## Use workflow functions in order
	tryCatch({
				x <- cel2raw(x, cels = cels, lowMemory = lowMemory, ...)
				x <- raw2long(x, set = set, normalize = normalize,
						log.intensities = log.intensities, tidy = TRUE, ...)
				x <- long2short(x, tidy = TRUE, ...)
				x <- short2ras(x, subtractMismatch = subtractMismatch, tidy = TRUE, ...)
				return(x)
			}, error = function(e) {
				print(e)
			}, interrupt = function(i) {
				print(i)
			})
}

## Workflow function wrapper. Converts SNPMaP objects from CEL files to 'rasS' format
## Arguments: SNPMaP object, [character, numeric, logical, logical, function, ...]
## Returns: SNPMaP object
`cel2rasS` <-
		function (x, cels = x@chps, lowMemory = x@lowMemory, set = x@set, normalize = x@normalize,
				log.intensities = x@logInt, subtractMismatch = x@useMM, FUN = x@summary, ...) 
{
	## Catch call for majorHistory slot, add datetime
	ca <- deparse(match.call(), width.cutoff = 500)
	d <- date()
	ca <- c(ca, d)
	if (!is(x, "SNPMaP")) 
		stop(deparse(substitute(x)), " is not a 'SNPMaP' object")
	x@majorHistory[[length(x@majorHistory) + 1]] <- ca
	## Use workflow functions in order
	tryCatch({
				x <- cel2raw(x, cels = cels, lowMemory = lowMemory, ...)
				x <- raw2long(x, set = set, normalize = normalize,
						log.intensities = log.intensities, tidy = TRUE, ...)
				x <- long2short(x, tidy = TRUE, ...)
				x <- short2ras(x, subtractMismatch = subtractMismatch, tidy = TRUE, ...)
				x <- ras2rasS(x, FUN = FUN, tidy = TRUE, ...)
				return(x)
			}, error = function(e) {
				print(e)
			}, interrupt = function(i) {
				print(i)
			})
}

## Workflow function. Converts SNPMaP objects from CEL files to 'raw' format
## Arguments: SNPMaP object, [character, ...]
## Returns: SNPMaP object
`cel2raw` <- function (x, cels = x@chps, lowMemory = x@lowMemory, ...) 
{
	## Catch call for majorHistory slot, add datetime
	ca <- deparse(match.call(), width.cutoff = 500)
	da <- date()
	ca <- c(ca, da)
	if (!is(x, "SNPMaP")) 
		stop(deparse(substitute(x)), " is not a 'SNPMaP' object")
	cat("Checking CEL files ...\n")
	flush.console()
	## Check CEL files
	ex <- file.exists(cels)
	if (any(!ex)) {
		if (length(cels[!ex]) == 1) {
			cat(cels[!ex], "does not exist, skipping ...\n")
			flush.console()
		}
		else {
			cat(cels[!ex], "do not exist, skipping ...\n")
			flush.console()
		}
		cels <- cels[ex]
	}
	cheque <- sapply(cels, function(y) {
				if (.isCel(y)) {
					z <- readCelHeader(y)
					cat(y, " OK; chiptype = ", z$chiptype, ", outliers = ", 
							z$noutliers, "\n", sep = "")
					flush.console()
					good <- TRUE
				}
				else {
					z <- list(chiptype = "unknown")
					cat(y, "is not a cel file; skipping ...\n")
					flush.console()
					good <- FALSE
				}
				c(z$chiptype, good)
			})
	if (any(cheque[1, as.logical(cheque[2, ])] != cheque[1, 1])) 
		stop("all CEL files must be of the same type")
	## Read CEL files
	cat("Reading ...\n")
	flush.console()
	cels <- cels[as.logical(cheque[2, ])]
	heado<-readCelHeader(cels[1])
	rowN <- heado$total
	colN <- length(cels)
	mat <- matrix()
	## Decide whether to make a matrix in memory or on disk
	if (lowMemory) {
		mat <- .makefile("raw", rowN, colN, x@tempDir)
		x@lowMemory <- TRUE
	}
	else {
		mat <- .gws(rowN, colN)
		if (!is.matrix(mat)) {
			warning("out of memory; switching to disk ...")
			mat <- .makefile("raw", rowN, colN, x@tempDir)
			x@lowMemory <- TRUE
		}
		else {
			x@lowMemory <- FALSE
		}
	}
	## Set progress bar
	pb<-ProgressBar(max = 40, stepLength = 40/length(cels))
	reset(pb)
	flush.console()
	## Loop over files and read intensities to matrix
	for (i in 1:length(cels)) {
		mat[, i] <- readCelIntensities(cels[i])
		increase(pb)
		flush.console()
	}
	## Correct slot values
	slot(x, 'snpdata', check=FALSE) <- mat
	rm(mat)
	gc()
	x@chiptype <- as.character(cheque[1, 1])
	x@chps <- cels
	x@cols <- cels
	x@celDim<-c(heado$rows, heado$cols)
	x@majorHistory[[length(x@majorHistory) + 1]] <- ca
	x@table <- "raw"
	close(x)
	return(x)
}

## Workflow function wrapper. Converts SNPMaP objects from CEL files to 'short' format
## Arguments: SNPMaP object, [character, numeric, logical, logical, ...]
## Returns: SNPMaP object
`cel2short` <-
		function (x, cels = x@chps, lowMemory = x@lowMemory, set = x@set, normalize = x@normalize,
				log.intensities = x@logInt, ...) 
{
	## Catch call for majorHistory slot, add datetime
	ca <- deparse(match.call(), width.cutoff = 500)
	d <- date()
	ca <- c(ca, d)
	if (!is(x, "SNPMaP")) 
		stop(deparse(substitute(x)), " is not a 'SNPMaP' object")
	x@majorHistory[[length(x@majorHistory) + 1]] <- ca
	## Use workflow functions in order
	tryCatch({
				x <- cel2raw(x, cels = cels, lowMemory = lowMemory, ...)
				x <- raw2long(x, set = set, normalize = normalize,
						log.intensities = log.intensities, tidy = TRUE, ...)
				x <- long2short(x, tidy = TRUE, ...)
				return(x)
			}, error = function(e) {
				print(e)
			}, interrupt = function(i) {
				print(i)
			})
}

## Workflow function wrapper. Converts SNPMaP objects from 'long' to 'ras' format
## Arguments: SNPMaP object, [logical, logical, ...]
## Returns: SNPMaP object
`long2ras` <-
		function (x, lowMemory = x@lowMemory, subtractMismatch = x@useMM, tidy = FALSE, ...)
{
	## Catch call for majorHistory slot, add datetime
	ca<-deparse(match.call(), width.cutoff=500)
	d<-date()
	ca<-c(ca, d)
	if(!is(x, "SNPMaP")) stop(deparse(substitute(x)), " is not a 'SNPMaP' object")
	x@majorHistory[[length(x@majorHistory)+1]]<-ca
	## Use workflow functions in order
	tryCatch({
				x<-long2short(x, lowMemory = lowMemory, tidy=tidy, ...)
				x<-short2ras(x, subtractMismatch = subtractMismatch, tidy=TRUE, ...)
				return(x)
			}, error=function(e){
				print(e)
			}, interrupt=function(i){
				print(i)
			})
}

## Workflow function wrapper. Converts SNPMaP objects from 'long' to 'rasS' format
## Arguments: SNPMaP object, [logical, function, logical, ...]
## Returns: SNPMaP object
`long2rasS` <-
		function (x, lowMemory = x@lowMemory, subtractMismatch = x@useMM, FUN = x@summary, tidy = FALSE, ...) 
{
	## Catch call for majorHistory slot, add datetime
	ca <- deparse(match.call(), width.cutoff = 500)
	d <- date()
	ca <- c(ca, d)
	if (!is(x, "SNPMaP")) 
		stop(deparse(substitute(x)), " is not a 'SNPMaP' object")
	if (x@table != "long") {
		warning("the data in ", deparse(substitute(x)), " is not in the 'long' format.\nAttempting to reformat ...")
		x <- raw2long(x, tidy = tidy, ...)
	}
	x@majorHistory[[length(x@majorHistory) + 1]] <- ca
	## Use workflow functions in order
	tryCatch({
				x <- long2short(x, lowMemory = lowMemory, tidy = tidy, ...)
				x <- short2ras(x, subtractMismatch = subtractMismatch, tidy = TRUE, ...)
				x <- ras2rasS(x, FUN = FUN, tidy = TRUE, ...)
				return(x)
			}, error = function(e) {
				print(e)
			}, interrupt = function(i) {
				print(i)
			})
}

## Workflow function. Converts SNPMaP objects from 'long' to 'short' format
## Arguments: SNPMaP object, [logical, ...]
## Returns: SNPMaP object
`long2short` <-
		function (x, lowMemory = x@lowMemory, tidy = FALSE, ...) 
{
	## Catch call for majorHistory slot, add datetime
	ca <- deparse(match.call(), width.cutoff = 500)
	da <- date()
	ca <- c(ca, da)
	if (!is(x, "SNPMaP")) 
		stop(deparse(substitute(x)), " is not a 'SNPMaP' object")
	if (x@table != "long") {
		warning("the data in ", deparse(substitute(x)), " is not in the 'long' format.\nAttempting to reformat ...")
		x <- raw2long(x, tidy = tidy, ...)
	}
	x <- .checkOpen(x)
	cat("Reformatting data to 'short' ...\n")
	flush.console()
	d <- dim(x@snpdata)
	co <- x@chps
	w <- 2 * x@width
	mat <- matrix()
	mm <- matrix()
	d2 <- c(d[1]/w, d[2] * w)
	rowN <- d2[1]
	colN <- d2[2]
	if (lowMemory) {
		mat <- .makefile("short", rowN, colN, x@tempDir)
		x@lowMemory <- TRUE
	}
	else {
		mat <- .gws(rowN, colN)
		if (!is.matrix(mat)) {
			warning("out of memory; switching to disk ...")
			mat <- .makefile("short", rowN, colN, x@tempDir)
			x@lowMemory <- TRUE
		}
		else {
			x@lowMemory <- FALSE
		}
	}
	pb<-ProgressBar(max = 40, stepLength = 40/(d[2]))
	reset(pb)
	flush.console()
	## Loop over matrix, changing dimensions to one row per SNP
	for (i in 1:d[2]) {
		mat[, (1 + ((i - 1) * w)):(i * w)] <- x@snpdata[, i]
		increase(pb)
		flush.console()
	}
	## Tidy up matrices on disk
	if (lowMemory & is(x@snpdata, "FileDoubleMatrix")) {
		close(x@snpdata)
		if (tidy) 
			unlink(x@snpdata$pathname)
	}
	slot(x, 'snpdata', check=FALSE) <- mat
	## Repeat with mismatch matrix if necessary
	if (!all(dim(x@mismatch)==c(0,0))) {
		if (lowMemory) {
			mm <- .makefile("shortMM", rowN, colN, x@tempDir)
		}
		else {
			mm <- .gws(rowN, colN)
			if (!is.matrix(mm)) {
				warning("out of memory; switching to disk ...")
				mm <- .makefile("shortMM", rowN, colN, x@tempDir)
				x@lowMemory <- TRUE
			}
		}
		if (d[2] != dim(x@mismatch)[2]) 
			stop("dimensions of mismatch data do not match dimensions of perfect match data in ", 
					deparse(substitute(x)))
		pb2<-ProgressBar(max = 40, stepLength = 40/(d[2]))
		reset(pb2)
		flush.console()
		for (i in 1:d[2]) {
			mm[, (1 + ((i - 1) * w)):(i * w)] <- x@mismatch[, i]
			increase(pb2)
			flush.console()
		}
		if (lowMemory & is(x@mismatch, "FileDoubleMatrix")) {
			close(x@mismatch)
			if (tidy) 
				unlink(x@mismatch$pathname)
		}
		slot(x, 'mismatch', check=FALSE) <- mm
	}
	## Generate column names
	w <- 2 * x@width
	x@cols <- paste(rep(x@chps, each = w), gl(2, w/2, w, c("a", 
							"b")), 1:(w/2), sep = "_")
	x@majorHistory[[length(x@majorHistory) + 1]] <- ca
	x@table <- "short"
	close(x)
	return(x)
}

## Workflow function. Converts SNPMaP objects from 'ras' to 'rasS' format
## Arguments: SNPMaP object, [function, logical, ...]
## Returns: SNPMaP object
`ras2rasS` <-
		function (x, lowMemory = x@lowMemory, FUN = x@summary, tidy = FALSE, ...) 
{
	## Catch call for majorHistory slot, add datetime
	ca <- deparse(match.call(), width.cutoff = 500)
	da <- date()
	ca <- c(ca, da)
	if (!is(x, "SNPMaP")) 
		stop(deparse(substitute(x)), " is not a 'SNPMaP' object")
	if (x@table != "ras") {
		warning("the data in ", deparse(substitute(x)), " is not in the 'ras' format.\nAttempting to reformat ...")
		x <- short2ras(x, tidy = tidy, ...)
	}
	cat("Summarizing Relative Allele Scores ...\n")
	flush.console()
	x <- .checkOpen(x)
	v <- x@width
	d <- dim(x@snpdata)
	rowN <- d[1]
	colN <- d[2]/v
	## Check summary function to ensure rasS consistency
	summ <- FUN(runif(v))
	if (length(summ) != 1) 
		stop("summary function does not return a single value from a vector of RAS scores")
	mat <- matrix()
	if (lowMemory) {
		mat <- .makefile("rasS", rowN, colN, x@tempDir)
		x@lowMemory <- TRUE
	}
	else {
		mat <- .gws(rowN, colN)
		if (!is.matrix(mat)) {
			warning("out of memory; switching to disk ...")
			mat <- .makefile("rasS", rowN, colN, x@tempDir)
			x@lowMemory <- TRUE
		}
		else {
			x@lowMemory <- FALSE
		}
	}
	pb<-ProgressBar(max = 40, stepLength = 40/colN)
	reset(pb)
	flush.console()
	## Loop over matrix, applying summary function to RAS scores
	## Return one column per array
	for (i in 1:colN) {
		mat[, i] <- apply((x@snpdata[, (1 + ((i - 1) * v)):(i * v)]), 1, FUN, ...)
		increase(pb)
		flush.console()
	}
	x@cols <- x@chps
	if (lowMemory & is(x@snpdata, "FileDoubleMatrix")) {
		close(x@snpdata)
		if (tidy) 
			unlink(x@snpdata$pathname)
	}
	x@summary <- FUN
	slot(x, 'snpdata', check=FALSE) <- mat
	x@majorHistory[[length(x@majorHistory) + 1]] <- ca
	x@table <- "rasS"
	close(x)
	return(x)
}

## Workflow function. Converts SNPMaP objects from 'raw' to 'long' format
## Arguments: SNPMaP object, [numeric, logical, logical, logical, ...]
## Returns: SNPMaP object
`raw2long` <-
		function (x, lowMemory = x@lowMemory, set = x@set, normalize = x@normalize, 
				log.intensities = x@logInt, tidy = FALSE, ...) 
{
	## Catch call for majorHistory slot, add datetime
	ca <- deparse(match.call(), width.cutoff = 500)
	da <- date()
	ca <- c(ca, da)
	if (!is(x, "SNPMaP")) 
		stop(deparse(substitute(x)), " is not a 'SNPMaP' object")
	if (x@table != "raw") 
		stop("the data in ", deparse(substitute(x)), " is not in the 'raw' format.\n")
	## Quantile normalise and log intensities if necessary
	if (normalize) {
		x <- norm(x, tidy = tidy)
	}
	if (log.intensities) {
		x <- logIntensities(x, tidy = ifelse(normalize, TRUE, 
						tidy))
	}
	x <- .checkOpen(x)
	## LazyLoad cdm matrix named by chiptype
	cdm <- eval(parse(text = x@chiptype))
	res <- matrix()
	mm <- matrix()
	rowN <- length(cdm[[set]]$perfect)
	colN <- length(x@chps)
	if (lowMemory) {
		res <- .makefile("long", rowN, colN, x@tempDir)
		x@lowMemory <- TRUE
	}
	else {
		res <- .gws(rowN, colN)
		if (!is.matrix(res)) {
			warning("out of memory; switching to disk ...")
			res <- .makefile("long", rowN, colN, x@tempDir)
			x@lowMemory <- TRUE
		}
		else {
			x@lowMemory <- FALSE
		}
	}
	cat("Reformatting data to 'long' ...\n")
	flush.console()
	if (set != x@set) x@set <- set
	pb<-ProgressBar(max = 40, stepLength = 40/colN)
	reset(pb)
	flush.console()
	## Loop over matrix extracting and reordering probe intensities
	## to group by SNP
	for (i in 1:colN) {
		y <- x@snpdata[, i]
		res[, i] <- y[cdm[[set]]$perfect]
		increase(pb)
		flush.console()
	}
	## Set 'snps' slot
	x@snps <- rownames(cdm[[set]]$perfect)
	x@width <- dim(cdm[[set]]$perfect)[2]/2
	## Repeat for mismatch matrix if necessary
	if (x@useMM && !is.null(cdm[[set]]$mismatch)) {
		if (lowMemory) {
			mm <- .makefile("longMM", rowN, colN, x@tempDir)
		}
		else {
			mm <- .gws(rowN, colN)
			if (!is.matrix(mm)) {
				warning("out of memory; switching to disk ...")
				mm <- .makefile("longMM", rowN, colN, x@tempDir)
				x@lowMemory <- TRUE
			}
		}
		pb2<-ProgressBar(max = 40, stepLength = 40/colN)
		reset(pb2)
		flush.console()
		for (i in 1:colN) {
			y <- x@snpdata[, i]
			mm[, i] <- y[cdm[[set]]$mismatch]
			increase(pb2)
			flush.console()
		}
		slot(x, 'mismatch', check=FALSE) <- mm
	}else{
		x@useMM <- FALSE
	}
	if (is(x@snpdata, "FileDoubleMatrix")) {
		close(x@snpdata)
		if (normalize | log.intensities | tidy) 
			unlink(x@snpdata$pathname)
	}
	slot(x, 'snpdata', check=FALSE) <- res
	x@majorHistory[[length(x@majorHistory) + 1]] <- ca
	x@table <- "long"
	rm(res, mm)
	gc()
	close(x)
	return(x)
}

## Workflow function wrapper. Converts SNPMaP objects from 'raw' to 'ras' format
## Arguments: SNPMaP object, [numeric, logical, logical, logical, logical, ...]
## Returns: SNPMaP object
`raw2ras` <-
		function (x, lowMemory = x@lowMemory, set = x@set, normalize = x@normalize,
				log.intensities = x@logInt, subtractMismatch = x@useMM, tidy = FALSE, ...) 
{
	## Catch call for majorHistory slot, add datetime
	ca <- deparse(match.call(), width.cutoff = 500)
	d <- date()
	ca <- c(ca, d)
	if (!is(x, "SNPMaP")) 
		stop(deparse(substitute(x)), " is not a 'SNPMaP' object")
	if (x@table != "raw") 
		stop("the data in ", deparse(substitute(x)), " is not in the 'raw' format.\n")
	x@majorHistory[[length(x@majorHistory) + 1]] <- ca
	## Use workflow functions in order
	tryCatch({
				x <- raw2long(x, lowMemory = lowMemory, set = set, normalize = normalize,
						log.intensities = log.intensities, tidy = tidy, ...)
				x <- long2short(x, tidy = TRUE, ...)
				x <- short2ras(x, subtractMismatch = subtractMismatch, tidy = TRUE, ...)
				return(x)
			}, error = function(e) {
				print(e)
			}, interrupt = function(i) {
				print(i)
			})
}

## Workflow function wrapper. Converts SNPMaP objects from 'raw' to 'rasS' format
## Arguments: SNPMaP object, [numeric, logical, logical, logical, function, logical, ...]
## Returns: SNPMaP object
`raw2rasS` <-
		function (x, lowMemory = x@lowMemory, set = x@set, normalize = x@normalize,
				log.intensities = x@logInt, subtractMismatch = x@useMM, FUN = x@summary, tidy = FALSE, ...) 
{
	## Catch call for majorHistory slot, add datetime
	ca <- deparse(match.call(), width.cutoff = 500)
	d <- date()
	ca <- c(ca, d)
	if (!is(x, "SNPMaP")) 
		stop(deparse(substitute(x)), " is not a 'SNPMaP' object")
	if (x@table != "raw") 
		stop("the data in ", deparse(substitute(x)), " is not in the 'raw' format.\n")
	x@majorHistory[[length(x@majorHistory) + 1]] <- ca
	## Use workflow functions in order
	tryCatch({
				x <- raw2long(x, lowMemory = lowMemory, set = set, normalize = normalize,
						log.intensities = log.intensities, tidy = tidy, ...)
				x <- long2short(x, tidy = TRUE, ...)
				x <- short2ras(x, subtractMismatch = subtractMismatch, tidy = TRUE, ...)
				x <- ras2rasS(x, FUN = FUN, tidy = TRUE, ...)
				return(x)
			}, error = function(e) {
				print(e)
			}, interrupt = function(i) {
				print(i)
			})
}

## Workflow function wrapper. Converts SNPMaP objects from 'raw' to 'short' format
## Arguments: SNPMaP object, [numeric, logical, logical, logical, logical, ...]
## Returns: SNPMaP object
`raw2short` <-
		function (x, lowMemory = x@lowMemory, set = x@set, normalize = x@normalize,
				log.intensities = x@logInt, tidy = FALSE, ...) 
{
	## Catch call for majorHistory slot, add datetime
	ca <- deparse(match.call(), width.cutoff = 500)
	d <- date()
	ca <- c(ca, d)
	if (!is(x, "SNPMaP")) 
		stop(deparse(substitute(x)), " is not a 'SNPMaP' object")
	if (x@table != "raw") 
		stop("the data in ", deparse(substitute(x)), " is not in the 'raw' format.\n")
	x@majorHistory[[length(x@majorHistory) + 1]] <- ca
	## Use workflow functions in order
	tryCatch({
				x <- raw2long(x, lowMemory = lowMemory, set = set, normalize = normalize,
						log.intensities = log.intensities, tidy = tidy, ...)
				x <- long2short(x, tidy = TRUE, ...)
				return(x)
			}, error = function(e) {
				print(e)
			}, interrupt = function(i) {
				print(i)
			})
}

## Workflow function. Converts SNPMaP objects from 'short' to 'ras' format
## Arguments: SNPMaP object, [logical, logical, ...]
## Returns: SNPMaP object
`short2ras` <-
		function (x, lowMemory = x@lowMemory, subtractMismatch = x@useMM, tidy = FALSE, ...) 
{
	## Catch call for majorHistory slot, add datetime
	ca <- deparse(match.call(), width.cutoff = 500)
	da <- date()
	ca <- c(ca, da)
	if (!is(x, "SNPMaP")) 
		stop(deparse(substitute(x)), " is not a 'SNPMaP' object")
	if (x@table != "short") {
		warning("the data in ", deparse(substitute(x)), " is not in the 'short' format.\nAttempting to reformat ...")
		x <- long2short(x, tidy = tidy, ...)
	}
	## Subtract mismatch probe intensities from perfect match if necessary
	if (subtractMismatch) {
		x <- minusMismatch(x, tidy = tidy)
	}
	cat("Calculating Relative Allele Scores (RAS) ...\n")
	flush.console()
	x <- .checkOpen(x)
	d <- dim(x@snpdata)
	rowN <- d[1]
	colN <- d[2]/2
	mat <- matrix()
	v <- x@width
	w <- 2 * v
	arrays <- d[2]/w
	if (lowMemory) {
		mat <- .makefile("ras", rowN, colN, x@tempDir)
		x@lowMemory <- TRUE
	}
	else {
		mat <- .gws(rowN, colN)
		if (!is.matrix(mat)) {
			warning("out of memory; switching to disk ...")
			mat <- .makefile("ras", rowN, colN, x@tempDir)
			x@lowMemory <- TRUE
		}
		else {
			x@lowMemory <- FALSE
		}
	}
	pb<-ProgressBar(max = 40, stepLength = 40/arrays)
	reset(pb)
	flush.console()
	## Loop over matrix chip by chip calculating RAS: A/A+B
	for (i in 1:arrays) {
		mat[, (1 + ((i - 1) * v)):(i * v)] <- .ra(x@snpdata[, (1 + ((i - 1) * w)):(i * w)])
		increase(pb)
		flush.console()
	}
	x@cols <- paste(rep(x@chps, each = w/2), 1:(w/2), sep = "_")
	if (lowMemory & is(x@snpdata, "FileDoubleMatrix")) {
		close(x@snpdata)
		if (subtractMismatch | tidy) 
			unlink(x@snpdata$pathname)
	}
	slot(x, 'snpdata', check=FALSE) <- mat
	x@majorHistory[[length(x@majorHistory) + 1]] <- ca
	x@table <- "ras"
	close(x)
	return(x)
}

## Workflow function wrapper. Converts SNPMaP objects from 'short' to 'rasS' format
## Arguments: SNPMaP object, [logical, function, logical, ...]
## Returns: SNPMaP object
`short2rasS` <-
		function (x, lowMemory = x@lowMemory, subtractMismatch = x@useMM, FUN = x@summary, tidy = FALSE, ...) 
{
	## Catch call for majorHistory slot, add datetime
	ca <- deparse(match.call(), width.cutoff = 500)
	d <- date()
	ca <- c(ca, d)
	if (!is(x, "SNPMaP")) 
		stop(deparse(substitute(x)), " is not a 'SNPMaP' object")
	if (x@table != "short") {
		warning("the data in ", deparse(substitute(x)), " is not in the 'short' format.\nAttempting to reformat ...")
		x <- long2short(x, tidy = tidy, ...)
	}
	x@majorHistory[[length(x@majorHistory) + 1]] <- ca
	## Use workflow functions in order
	tryCatch({
				x <- short2ras(x, lowMemory = lowMemory, subtractMismatch = subtractMismatch, tidy = tidy, ...)
				x <- ras2rasS(x, FUN = FUN, tidy = TRUE, ...)
				return(x)
			}, error = function(e) {
				print(e)
			}, interrupt = function(i) {
				print(i)
			})
}



