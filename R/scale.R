## scale.R
## Copyright (c) 2009 Leo Schalkwyk and Oliver SP Davis

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

scaleq <- function (x, targets=list(), quant=c(.90,.95), lowMemory = x@lowMemory, tidy = FALSE, 
    ...) 
{
   ca <- deparse(match.call(), width.cutoff = 500)
   da <- date()
   ca <- c(ca, da)
   if (!is(x, "SNPMaP")) 
      stop("this is not a 'SNPMaP' object")
   n <- dim(x@snpdata)
   x <- .checkOpen(x)
   mat <- mate(x, lowMemory)
   cat("extracting quantiles ...\n")
   flush.console()
   toto <- 0
   pb <- ProgressBar(max = 40, stepLength = 40/(n[2]))
   reset(pb)
   flush.console()
   for (i in 1:n[2]) {  # don't know if another method is more 
                      # efficient for a filedoublematrix 
      toto <- toto + x@snpdata[, i]
      increase(pb)
      flush.console()
   }
   dim(toto) <- NULL
   qq <- quantile(toto, quant)# / n[2]
   (toto > qq[1]) & 
   (toto < qq[2])-> stnu 
   sta<- sum(toto[stnu])/n[2]
   stnu <- (1:length(toto))[stnu] # can't index a filedoublematrix with 
   qqm <- qq/n[2]                 # a logical, apparently
   mat <- sca(x,mat,sta,stnu)
   if(length(targets) ==0 ){
      if (is(x@snpdata, "FileDoubleMatrix")) {
         close(x@snpdata)
         if (tidy) 
            unlink(x@snpdata$pathname)
      }
      slot(x, "snpdata", check = FALSE) <- mat
      x@normalize <- FALSE
      x@majorHistory[[length(x@majorHistory) + 1]] <- ca
      x@transformation <- append(x@transformation, "scale")
      close(x)
      return(x)
      }
      else{
         cat ('multiple snpmap objects to be scaled - will return a list\n')
         retthing <- lapply(targets,function(tar){
         m <- mate(tar, lowMemory)
         m <- sca(tar,m,sta,stnu)
         if (is(tar@snpdata, "FileDoubleMatrix")) {
            close(tar@snpdata)
            if (tidy) 
                unlink(tar@snpdata$pathname)
         }
         slot(tar, "snpdata", check = FALSE) <- m
         tar@normalize <- FALSE
         tar@majorHistory[[length(tar@majorHistory) + 1]] <- ca
         tar@transformation <- append(tar@transformation, "scale")
         close(tar)
         tar
         }
      )
   }
   append(x,retthing)  
}



mate <- function(x,lowMem){
    mat <- matrix()
    n <- dim(x@snpdata)
    if (lowMem) {
        mat <- .makefile(x@table, n[1], n[2], x@tempDir)
    }
    else {
        mat <- .gws(n[1], n[2])
        gc()
        if (!is.matrix(mat)) {
            warning("out of memory; switching to disk ...")
            mat <- .makefile(x@table, n[1], n[2], x@tempDir)
            x@lowMemory <- TRUE
        }
    }
mat
}

sca <- function (x,mat,sta,stnu){
    x <- .checkOpen(x) # todo: sanity checks on x:
                       # chiptype, table same as initial obj
    n<- dim(mat)
    cat("scaling", "...\n")
    pb2 <- ProgressBar(max = 40, stepLength = 40/(n[2]))
    reset(pb2)
    flush.console()
    for (i in 1:n[2]) {
        mult <- sta/ sum(x@snpdata[stnu,i])
        mat[, i] <- x@snpdata[,i] * mult
        increase(pb2)
        flush.console()
    }
    mat
}

