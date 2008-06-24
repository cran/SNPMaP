## getSNPMaP.R
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

## Function to access SNPMaP slots. Arguments: SNPMaP object
## Returns: list of slot values
`getSNPMaP` <- function (x, ...) 
{
	slots <- as.character(c(...))
	slots <- match.arg(slots, c('useMM', 'normalize', 'logInt', 'summary',
					'lowMemory', 'tempDir', 'table', 'chiptype', 'celDim', 'set',
					'snps', 'chps', 'cols', 'width', 'transformation', 'experiment',
					'created', 'version', 'majorHistory'), several.ok=TRUE)
	if(length(slots)==0) stop('no valid arguments')
	y <- list()
	for(i in 1:length(slots)){
		y[[i]] <- slot(x, slots[i])
		names(y)[i] <- slots[i]
	}
	return(y)
}

## Function to set SNPMaP slots. Arguments: SNPMaP object
## Returns: a SNPMaP object
`setSNPMaP` <- function (x, namedList) 
{
	if(!is.list(namedList)) stop("'namedList' argument is not a list")
	slots <- names(namedList)
	slots <- match.arg(slots, c('useMM', 'normalize', 'logInt', 'summary',
					'lowMemory', 'tempDir', 'set', 'experiment'), several.ok=TRUE)
	namedList <- namedList[slots]
	if(length(namedList)==0) stop('no valid arguments')
	y <- list()
	for(i in 1:length(namedList)){
		y[[i]] <- slot(x, names(namedList)[i])
		names(y)[i] <- names(namedList)[i]
		slot(x, names(namedList)[i], check = TRUE) <- namedList[[i]]
	}
	return(x)
}
