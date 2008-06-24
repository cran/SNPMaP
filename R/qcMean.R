## qcMean.R
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

## Default function to summarise RAS scores to a single
## number per SNP per array. Arguments: numeric; RAS scores from a
## single SNP on one array. Returns: single numeric (a mean if
## two thirds of the probes are not NA)
`qcMean` <-
		function (x) 
{
	d <- length(x)
	if ((sum(!is.na(x))/d) < (2/3)) {
		y <- NA
	}
	else {
		y <- mean(x, na.rm = TRUE)
	}
	return(y)
}

