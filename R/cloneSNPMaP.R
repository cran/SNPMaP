## cloneSNPMaP.R
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

## Function to copy data from SNPMaP object on disk, rather than
## only copying filehandles. Arguments: SNPMaP object
## Returns: copy of SNPMaP object
`cloneSNPMaP` <-
		function (x, ...) 
{
	ca <- deparse(match.call(), width.cutoff = 500)
	da <- date()
	ca <- c(ca, da)
	if (!is(x, "SNPMaP")) 
		stop(deparse(substitute(x)), " is not a 'SNPMaP' object")
	if (!x@lowMemory) 
		return(x)
	x <- .checkOpen(x)
	y <- x
	mat <- matrix()
	mm <- matrix()
	if (is(x@snpdata, "FileDoubleMatrix")) {
		d <- dim(x@snpdata)
		mat <- .makefile(y@table, d[1], d[2], y@tempDir)
		for (i in 1:d[2]) {
			mat[, i] <- x@snpdata[, i]
		}
		slot(y, 'snpdata', check=FALSE) <- mat
	}
	if (is(x@mismatch, "FileDoubleMatrix")) {
		d <- dim(x@mismatch)
		mm <- .makefile(paste(y@table, "MM", sep = ""), d[1], 
				d[2], y@tempDir)
		for (i in 1:d[2]) {
			mm[, i] <- x@mismatch[, i]
		}
		slot(y, 'mismatch', check=FALSE) <- mm
	}
	y@majorHistory[[length(y@majorHistory) + 1]] <- ca
	close(x)
	close(y)
	return(y)
}
