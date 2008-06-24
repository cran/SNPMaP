## disk2memory.R
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

## Function to move SNPMaP object from disk to memory
## Arguments: SNPMaP object [logical, ...]
## Returns: SNPMaP object in memory
`disk2memory` <-
		function (x, tidy = FALSE, ...) 
{
	ca <- deparse(match.call(), width.cutoff = 500)
	da <- date()
	ca <- c(ca, da)
	if (!is(x, "SNPMaP")) 
		stop(deparse(substitute(x)), " is not a 'SNPMaP' object")
	if (!x@lowMemory) 
		stop(deparse(substitute(x)), " is already in memory")
	x <- .checkOpen(x)
	x@lowMemory <- FALSE
	mat <- x@snpdata[, ]
	close(x@snpdata)
	if (tidy) 
		unlink(x@snpdata$pathname)
	slot(x, 'snpdata', check=FALSE) <- mat
	if (!all(dim(x@mismatch)==c(0,0))) {
		mat2 <- x@mismatch[, ]
		close(x@mismatch)
		if (tidy) 
			unlink(x@mismatch$pathname)
		slot(x, 'mismatch', check=FALSE) <- mat2
	}
	x@majorHistory[[length(x@majorHistory) + 1]] <- ca
	return(x)
}

## Function to move SNPMaP object from memory to disk
## Arguments: SNPMaP object [logical, ...]
## Returns: SNPMaP object on disk
`memory2disk` <-
		function (x, ...) 
{
	ca <- deparse(match.call(), width.cutoff = 500)
	da <- date()
	ca <- c(ca, da)
	if (!is(x, "SNPMaP")) 
		stop(deparse(substitute(x)), " is not a 'SNPMaP' object")
	if (x@lowMemory) 
		stop(deparse(substitute(x)), " is already on disk")
	x@lowMemory <- TRUE
	d <- dim(x@snpdata)
	mat <- .makefile(x@table, d[1], d[2], x@tempDir)
	mat[, ] <- x@snpdata[, ]
	slot(x, 'snpdata', check=FALSE) <- mat
	if (!all(dim(x@mismatch)==c(0,0))) {
		mat2 <- .makefile(paste(x@table, "MM", sep = ""), d[1], 
				d[2], x@tempDir)
		mat2[, ] <- x@mismatch[, ]
		slot(x, 'mismatch', check=FALSE) <- mat2
	}
	x@majorHistory[[length(x@majorHistory) + 1]] <- ca
	close(x)
	return(x)
}

