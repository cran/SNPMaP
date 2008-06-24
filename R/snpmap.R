## snpmap.R
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

## Wrapper function to prepare SNPMaP data for analysis. Arguments: [character,
## logical, numeric, logical, logical, logical, function, character, character,
## logical, ...]. Returns: SNPMaP object (warns and retruns list if more than
## one set is requested)
`snpmap` <-
		function (cels = dir(pattern = ".[cC][eE][lL]$"), lowMemory = TRUE, 
				set = 1, useMM = FALSE, normalize = FALSE, log.intensities = FALSE, 
				ras.summary = qcMean, tempDir = ".", RUN = 'cel2rasS', interactive=FALSE, ...) 
{
	## Catch call for majorHistory slot
	ca <- deparse(match.call(), width.cutoff = 500)
	da <- date()
	ca <- c(ca, da)
	## If more than one set is requested, throw warning, use msnpmap instead and
	## return a list of SNPMaP objects
	if(set==0 | length(set)>1){
		warning("more than one 'set' requested; calling msnpmap() and returning a list")
		object<-msnpmap(cels=cels, lowMemory=lowMemory, set=set, useMM=useMM,
				normalize=normalize, log.intensities=log.intensities, ras.summary=ras.summary,
				tempDir=tempDir, RUN=RUN, interactive=interactive, ...)
		return(object)
	}
	object <- new("SNPMaP")
	## Set slots
	object@majorHistory[[length(object@majorHistory) + 1]] <- ca
	object@lowMemory <- lowMemory
	object@tempDir <- tempDir
	object@set <- set
	object@useMM <- useMM
	object@normalize <- normalize
	object@logInt <- log.intensities
	## Interactive file choice only available on Windows at the moment
	if(interactive){
		if(exists('choose.files')){
			cels<-choose.files(
					caption='choose CEL files for processing',
					filters=rbind(c("CEL files (*.CEL, *.cel)", "*.cel; *.CEL"), c("All files (*.*)", "*.*")),
					index=1
					)
		}else{
			stop("interactive=TRUE is Windows-specific")
		}
	}
	object@chps <- cels
	object@summary <- ras.summary
	object@table <- "empty"
	object@chiptype <- "unknown"
	## Catch workflow function and evaluate
	RUN<-match.arg(RUN, c('cel2raw', 'cel2long', 'cel2short', 'cel2ras', 'cel2rasS'))
	RUN<-eval(parse(text=RUN))
	object <- RUN(object, ...)
	return(object)
}

## Wrapper function to prepare SNPMaP data for analysis; can handle more than one set
## Arguments: [character, logical, numeric, logical, logical, logical, function, character, character,
## logical, ...]. Returns: list of SNPMaP objects
`msnpmap` <-
		function (cels = dir(pattern = ".[cC][eE][lL]$"), lowMemory = TRUE, 
				set = 0, useMM = FALSE, normalize = FALSE, log.intensities = FALSE, 
				ras.summary = qcMean, tempDir = ".", RUN = 'cel2rasS', interactive=FALSE, ...) 
{
	## Catch call for majorHistory slot
	ca <- deparse(match.call(), width.cutoff = 500)
	da <- date()
	ca <- c(ca, da)
	object <- new("SNPMaP")
	## Set slots
	object@majorHistory[[length(object@majorHistory) + 1]] <- ca
	object@lowMemory <- lowMemory
	object@tempDir <- tempDir
	object@set <- set
	object@useMM <- useMM
	object@normalize <- normalize
	object@logInt <- log.intensities
	## Interactive file choice only available on Windows at the moment
	if(interactive){
		if(exists('choose.files')){
			cels<-choose.files(
					caption='choose CEL files for processing',
					filters=rbind(c("CEL files (*.CEL, *.cel)", "*.cel; *.CEL"), c("All files (*.*)", "*.*")),
					index=1
					)
		}else{
			warning("interactive=TRUE is Windows-specific; defaulting to 'cels' argument")
		}
	}
	object@chps <- cels
	object@summary <- ras.summary
	object@table <- "empty"
	object@chiptype <- "unknown"
	object<-cel2raw(object, ...)
	## Catch workflow function (consistent with snpmap()))
	RUN<-match.arg(RUN, c('cel2raw', 'cel2long', 'cel2short', 'cel2ras', 'cel2rasS'))
	if(RUN=='cel2raw') return(list(all=object))
	## Decide on workflow function long2*
	if(normalize) object<-norm(object, tidy=TRUE)
	if(log.intensities) object<-logIntensities(object, tidy=TRUE)
	RUN<-switch(RUN, 'cel2long'=function(x){x}, 'cel2short'=long2short,
			'cel2ras'=long2ras, 'cel2rasS'=long2rasS)
	## Deal with the special 'set=0' (all sets)
	if(set==0){
		posits<-1:length(eval(parse(text=object@chiptype)))
	}else{
		posits<-set
	}
	snpmplist<-list()
	## Loop over sets, converting to 'long' format
	for(i in 1:length(posits)){
		cat('Set ', posits[i], ':\n', sep='')
		object@set<-posits[i]
		snpmplist[[i]] <- raw2long(object, tidy=FALSE, log.intensities=FALSE, normalize=FALSE, ...)
		names(snpmplist)[i]<-names(eval(parse(text=object@chiptype)))[posits[i]]
	}
	## Tidy up 'raw' object on disk
	if(is(object@snpdata, 'FileDoubleMatrix')){
		if(isOpen(object@snpdata)) close(object@snpdata)
		unlink(object@snpdata$pathname)
	}
	## Loop over sets running workflow long2*
	for(i in 1:length(snpmplist)){
		cat('Set ', snpmplist[[i]]@set, ':\n', sep='')
		snpmplist[[i]]<-RUN(snpmplist[[i]], tidy=TRUE, ...)
	}
	return(snpmplist)
}

