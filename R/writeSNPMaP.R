## writeSNPMaP.R
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

## write.table function
## Write SNPMaP data to a file or connection
## Arguments: SNPMaP object, [character, logical, logical, character, character, logical]
## Returns: TRUE invisibly
`writeSNPMaP` <- function(x, file="", mm=FALSE, sep="\t", dec=".", transpose=TRUE, ...)
{
		x<-.checkOpen(x)
		if(file!=""){
			file <- file(
					file,
					open='w'
					)
		}
		if(mm){
			data <- x@mismatch
		}else{
			data <- x@snpdata
		}
		if(all(dim(data)==c(0,0))) stop('no data to write')
		op <- options(OutDec=dec)
		if(!transpose){
			cat("", as.character(x@cols), file = file, sep = sep, append = TRUE)
			flush.console()
			for(i in 1:dim(data)[1]){
				cat('\n', file = file, append = TRUE)
				cat(shQuote(as.character(x@snps[i])), as.character(data[i,]), file = file, sep = sep, append = TRUE)
				flush.console()
			}
		}else{
			cat("", as.character(x@snps), file = file, sep = sep, append = TRUE)
			flush.console()
			for(i in 1:dim(data)[2]){
				cat('\n', file = file, append = TRUE)
				cat(shQuote(as.character(x@cols[i])), as.character(data[,i]), file = file, sep = sep, append = TRUE)
				flush.console()
			}
		}
		cat('\n', file = file, append = TRUE)
		flush.console()
		if(inherits(file, 'connection')){
			flush(file)
			close(file)
		}
		options(op)
		return(invisible(TRUE))
}

`readSNPMaP` <- function(file="", sep="\t", dec=".", transpose=TRUE, ...)
{
	if(transpose){
		con<-file(file, open="r")
		head<-scan(file=con, what='character', nlines=1, sep=sep, quiet=TRUE)
		body<-matrix(scan(file=con, what='character', sep=sep, quiet=TRUE),nrow=length(head))
		close(con)
		top<-body[1,]
		body<-body[-1,]
		if(dec!='.') body<-sub(dec, '.', body)
		storage.mode(body)<-'double'
		dims<-list(head[-1], top)
		dimnames(body)<-dims
		return(body)
	}else{
		con<-file(file, open="r")
		head<-scan(file=con, what='character', nlines=1, sep=sep, quiet=TRUE)
		body<-matrix(scan(file=con, what='character', sep=sep, quiet=TRUE), byrow=TRUE, ncol=length(head))
		close(con)
		top<-body[,1]
		body<-body[,-1]
		if(dec!='.') body<-sub(dec, '.', body)
		storage.mode(body)<-'double'
		dims<-list(top, head[-1])
		dimnames(body)<-dims
		return(body)
	}
}
