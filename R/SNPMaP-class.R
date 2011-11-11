## SNPMaP-class.R
## Copyright (c) 2008 Oliver SP Davis and Leo Schalkwyk
##
## Distributed under the GNU General Public licence version 3 or later
##
## This program is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program.  If not, see <http://www.gnu.org/licenses/>.
##
##
##
## S4 classes for the SNP Microarrays and Pooling (SNPMaP) technique
##
## TODO: 1)    Quality control?
## TODO: 2)    Adapt ras2rasS to allow calculation of arbitary stats in additional columns
## TODO: 3)    Adapt cdx files for other arrays
## TODO: 4)    Add genomic position vector to cdx to enable sliding window analyses?

## Incorporate the S3 class FileDoubleMatrix from R.huge to make low RAM matrices
setOldClass("FileDoubleMatrix")

## Make a datatype for the snpdata slot that will take either
## a standard matrix or a memory-mapped matrix
setClassUnion("SNPMaPdata", c("matrix", "FileDoubleMatrix"))

## Define the SNPMaP class that will form the basis of the objects
setClass("SNPMaP",
    representation(
        snpdata=    "SNPMaPdata",
        mismatch=    "SNPMaPdata",
        useMM=        "logical",
        normalize=    "logical",
        logInt=        "logical",
        summary=    "function",
        lowMemory=    "logical",
        tempDir=    "character",
        table=        "character",
        chiptype=    "character",
        celDim=      "numeric",
        set=    "numeric",
        snps=        "character",
        chps=        "character",
        cols=        "character",
        width=        "numeric",
        transformation="character",
        experiment=     "factor",
        annotation=    "list",
        created=    "character",
        version=  "character",
        majorHistory= "list"
        ),
    prototype(
        snpdata=    new("matrix"),
        mismatch=    new("matrix"),
        useMM=        logical(0),
        summary=    new("function"),
        normalize=    logical(0),
        logInt=        logical(0),
        lowMemory=    logical(0),
        tempDir=    character(0),
        table=        character(0),
        chiptype=    character(0),
        celDim=      numeric(0),
        set=    numeric(0),
        snps=        character(0),
        chps=        character(0),
        cols=        character(0),
        width=        numeric(0),
        transformation=character(0),
        experiment=new("factor"),
        annotation=  list(),
        created=    date(),
        version=  sessionInfo(package='SNPMaP')$other[[1]]$Version,
        majorHistory= list()
        )
    )


##setClass(Class='mSNPMaP', representation='list')

## Set methods for indexing a SNPMaP object
## Test for the type of data structure and returns a named matrix
## x[1:10,]
## Arguments: ANY, ANY, [logical, logical]
## Returns: SNPMaP object
setGeneric("[")

## x[,]
setMethod("[",
    signature(x="SNPMaP", i="missing", j="missing", drop="missing"),
    function(x){
      if (!x@lowMemory) 
        return(x)
      x <- .checkOpen(x)
      y <- x
      mat <- matrix()
      mm <- matrix()
      if (is(x@snpdata, "FileDoubleMatrix")) {
        d <- dim(x@snpdata)
        mat <- .makefile(y@table, d[1], d[2], y@tempDir)
        for (k in 1:d[2]) {
          mat[, k] <- x@snpdata[, k]
        }
        y@snpdata <- mat
      }
      if (is(x@mismatch, "FileDoubleMatrix")) {
        d <- dim(x@mismatch)
        mm <- .makefile(paste(y@table, "MM", sep = ""), d[1], 
            d[2], y@tempDir)
        for (k in 1:d[2]) {
          mm[, k] <- x@mismatch[, k]
        }
        slot(y, 'mismatch', check=FALSE) <- mm
      }
      close(x)
      close(y)
      return(y)
    }
    )

## x[1:10,]
setMethod("[",
    signature(x="SNPMaP", i="ANY", j="missing", drop="missing"),
    function(x, i){
      ca <- paste('subset', length(i), 'rows')
      da <- date()
      ca <- c(ca, da)
      if(!any(x@table==c('short', 'ras', 'rasS')))
        stop("indexing rows currently works only with 'short', 'ras' or 'rasS' formats")
      x<-.checkOpen(x)
      if(is.character(i)){
        if(
            any(
                is.na(
                    i<-match(i, x@snps)
                    )
                )
            )stop("undefined rows selected")
      }
      if (!x@lowMemory) {
        slot(x, 'snpdata', check=FALSE)<-x@snpdata[i,,drop=FALSE]
        x@snps<-x@snps[i]
        if(!all(dim(x@mismatch)==c(0,0))){
          slot(x, 'mismatch', check=FALSE)<-x@mismatch[i,,drop=FALSE]
        }
        x@majorHistory[[length(x@majorHistory) + 1]] <- ca
        return(x)
      } else {
        x <- .checkOpen(x)
        y <- x
        mat <- matrix()
        mm <- matrix()
        if (is(x@snpdata, "FileDoubleMatrix")) {
          d <- dim(x@snpdata)
          mat <- .makefile(y@table, length(i), d[2], y@tempDir)
          for (k in 1:d[2]) {
            mat[i, k] <- x@snpdata[i, k]
          }
          y@snpdata <- mat
        }
        if (is(x@mismatch, "FileDoubleMatrix")) {
          d <- dim(x@mismatch)
          mm <- .makefile(paste(y@table, "MM", sep = ""), length(i), 
              d[2], y@tempDir)
          for (k in 1:d[2]) {
            mm[i, k] <- x@mismatch[i, k]
          }
          slot(y, 'mismatch', check=FALSE) <- mm
        }
        y@snps<-y@snps[i]
        y@majorHistory[[length(y@majorHistory) + 1]] <- ca
        close(x)
        close(y)
        return(y)
      }
    }
    )

## x[,1:10]
setMethod("[",
    signature(x="SNPMaP", i="missing", j="ANY", drop="missing"),
    function(x, j){
      ca <- paste('subset', length(j), 'columns')
      da <- date()
      ca <- c(ca, da)
      if(!any(x@table==c('raw', 'long', 'rasS')))
        stop("indexing columns currently works only with 'raw', 'long' or 'rasS' formats")
      x<-.checkOpen(x)
      if(is.character(j)){
        if(
            any(
                is.na(
                    j<-match(j, x@chps)
                    )
                )
            )stop("undefined columns selected")
      }
      if (!x@lowMemory) {
        slot(x, 'snpdata', check=FALSE)<-x@snpdata[,j,drop=FALSE]
        x@cols<-x@cols[j]
        x@chps<-x@chps[j]
        if(!all(dim(x@mismatch)==c(0,0))){
          slot(x, 'mismatch', check=FALSE)<-x@mismatch[,j,drop=FALSE]
        }
        x@majorHistory[[length(x@majorHistory) + 1]] <- ca
        return(x)
      } else {
        x <- .checkOpen(x)
        y <- x
        mat <- matrix()
        mm <- matrix()
        if (is(x@snpdata, "FileDoubleMatrix")) {
          d <- dim(x@snpdata)
          mat <- .makefile(y@table, d[1], length(j), y@tempDir)
          for (k in 1:length(j)) {
            mat[, k] <- x@snpdata[, j[k]]
          }
          y@snpdata <- mat
        }
        if (is(x@mismatch, "FileDoubleMatrix")) {
          d <- dim(x@mismatch)
          mm <- .makefile(paste(y@table, "MM", sep = ""), d[1], 
              length(j), y@tempDir)
          for (k in 1:length(j)) {
            mm[, k] <- x@mismatch[, j[k]]
          }
          slot(y, 'mismatch', check=FALSE) <- mm
        }
        y@cols<-y@cols[j]
        y@chps<-y@chps[j]
        y@experiment<-as.factor(as.character(y@experiment)[j])
        y@majorHistory[[length(y@majorHistory) + 1]] <- ca
        close(x)
        close(y)
        return(y)
      }
    }
    )

## x[1:10,1:10]
setMethod("[",
    signature(x="SNPMaP", i="ANY", j="ANY", drop="missing"),
    function(x, i, j){
      ca <- paste('subset', length(i), 'rows and', length(j), 'columns')
      da <- date()
      ca <- c(ca, da)
      if(!any(x@table==c('rasS')))
        stop("indexing both rows and columns currently works only with the 'rasS' format")
      x<-.checkOpen(x)
      if(is.character(i)){
        if(
            any(
                is.na(
                    i<-match(i, x@snps)
                    )
                )
            )stop("undefined rows selected")
      }
      if(is.character(j)){
        if(
            any(
                is.na(
                    j<-match(j, x@chps)
                    )
                )
            )stop("undefined columns selected")
      }
      if (!x@lowMemory) {
        slot(x, 'snpdata', check=FALSE)<-x@snpdata[i,j,drop=FALSE]
        x@snps<-x@snps[i]
        x@cols<-x@cols[j]
        x@chps<-x@chps[j]
        if(!all(dim(x@mismatch)==c(0,0))){
          slot(x, 'mismatch', check=FALSE)<-x@mismatch[i,j,drop=FALSE]
        }
        x@majorHistory[[length(x@majorHistory) + 1]] <- ca
        return(x)
      } else {
        x <- .checkOpen(x)
        y <- x
        mat <- matrix()
        mm <- matrix()
        if (is(x@snpdata, "FileDoubleMatrix")) {
          d <- dim(x@snpdata)
          mat <- .makefile(y@table, length(i), length(j), y@tempDir)
          for (k in 1:length(j)) {
            mat[i, k] <- x@snpdata[i, j[k]]
          }
          slot(y, 'snpdata', check=FALSE) <- mat
        }
        if (is(x@mismatch, "FileDoubleMatrix")) {
          d <- dim(x@mismatch)
          mm <- .makefile(paste(y@table, "MM", sep = ""), length(i), 
              length(j), y@tempDir)
          for (k in 1:length(j)) {
            mm[i, k] <- x@mismatch[i, j[k]]
          }
          slot(y, 'mismatch', check=FALSE) <- mm
        }
        y@snps<-y@snps[i]
        y@cols<-y@cols[j]
        y@chps<-y@chps[j]
        y@experiment<-as.factor(as.character(y@experiment)[j])
        y@majorHistory[[length(y@majorHistory) + 1]] <- ca
        close(x)
        close(y)
        return(y)
      }
    }
    )


## Specify a generator method to allow for easy subclassing in the future
## Arguments: SNPMaP object. Returns: SNPMaP object
setMethod("initialize",
    signature(.Object="SNPMaP"),
    function(.Object, ...){
      .Object<-callNextMethod()
      return(.Object)
    }
    )

## Define a summary method
## (returns an S3 object so the print method works correctly)
## Arguments: SNPMaP object
## Returns: summary.SNPMaP S3 object
setGeneric("summary")

setMethod("summary",
    signature(object="SNPMaP"),
    function(object, ...){
      x<-object
      y<-list()
      y$objName<-deparse(substitute(object))
      y$dimensions<-dim(x@snpdata)
      y$lowMemory<-x@lowMemory
      y$file<-character(0)
      y$mismatchFile<-character(0)
      
      if(x@lowMemory){
        if(is(x@snpdata, "FileDoubleMatrix")){
          y$file<-x@snpdata$pathname
        }
        
        if(is(x@mismatch, "FileDoubleMatrix")){
          y$mismatchFile<-x@mismatch$pathname
        }
        
      }
      
      y$table            <-x@table
      y$chiptype        <-x@chiptype
      y$snps            <-length(x@snps)
      y$numberOfArrays<-length(x@chps)
      y$arrayNames    <-x@chps
      y$created        <-x@created
      y$version    <-x@version
      y$history        <-x@majorHistory
      y$comment        <-comment(x)
      class(y)        <-"summary.SNPMaP"
      return(y)
    }
    )

## What prints when the summary is evaluated (a summary paragraph)
## Arguments: summary.SNPMaP S3 object
## Returns: summary.SNPMaP S3 object invisibly
`print.summary.SNPMaP` <-
    function (x, ...) 
{
  mem <- "in memory"
  if (length(x$file > 0)) {
    mem <- paste("on disk\nat", x$file)
  }
  if (length(x$mismatchFile) > 0) {
    mem <- paste(mem, "and on disk at", x$mismatchFile)
  }
  cat("\n")
  cat(x$objName, "is a version", x$version, "SNPMaP S4 object created on", x$created, 
      "\ncomprising data on", switch(x$table, empty = "no", 
          raw = "all the", x$snps), "SNPs from the following", 
      x$numberOfArrays, x$chiptype, ifelse(x$numberOfArrays == 
              1, "array", "arrays"), "\nin a", x$dimensions[1], 
      "by", x$dimensions[2], "table of", x$table, "format", 
      mem)
  cat(":\n\n")
  cat(x$arrayNames, sep = "\n")
  cat("\nMajor history:\n")
  if (length(x$history) > 0) {
    print(noquote(matrix(unlist(x$history), ncol = 2, byrow = TRUE, 
                dimnames = list(c(1:length(x$history)), c("Function", 
                        "Date")))))
  }
  else {
    cat("None.\n")
  }
  cat("\nComment:\n")
  cat(x$comment, "\n\n")
  invisible(x)
}

## Plot method draws density of intensities/RAS for each chip, PM only
## Each line is in a different style and is coloured differently
## See 'boxplot' for a boxplot method
## Arguments: SNPMaP object [function, numeric, numeric, character, character,
## character, character, character, numeric, logical, ...]
## Returns: SNPMaP object invisibly
setGeneric("plot")

setMethod("plot",
    signature(x="SNPMaP", y="missing"),
    function(
        x,
        FUN=function(x){x},
        xlim=c(loX, hiX),
        ylim=c(0, hiY),
        xlab="guess",
        ylab="Density",
        main="",
        col=rainbow(length(x@chps)),
        legend.position="left",
        legend.bty="n",
        lty=1:length(x@chps),
        zero.line=TRUE,
        ...
        ){
      if(x@table=="empty") stop(
            "this SNPMaP object is of 'empty' format"
            )
      x<-.checkOpen(x)
      y<-vector("list")
      
      if(xlab=="guess"){
        xlab<-switch(
            x@table,
            "raw"=,
            "long"=,
            "short"="Intensity",
            "ras"=,
            "rasS"="Relative Allele Score (RAS)"
            )
      }
      
      w<-switch(
          x@table,
          "raw"=,
          "long"=,
          "rasS"=1,
          "ras"=x@width,
          "short"=x@width*2
          )
      
      for(i in 1:length(x@chps)){
        ## All the columns corresponding to one chip
        y[[i]]<-density(FUN(x@snpdata[,(1+((i-1)*w)):(i*w)]), na.rm=TRUE)
      }
      
      X<-sapply(y, function(x) x$x)
      Y<-sapply(y, function(x) x$y)
      hiY<-max(Y)
      hiX<-max(X)
      loX<-min(X)
      
      cols<-character(length(x@chps))
      cols[]<-col
      
      ltps<-integer(length(x@chps))
      ltps[]<-lty
      
      layout(matrix(c(1,2), 1, 2), c(2,1))
      old <- par(mar=c(5,4,4,1)+0.1)
      plot(
          X[,1],
          Y[,1],
          xlim=xlim,
          ylim=ylim,
          type="l",
          xlab=xlab,
          ylab=ylab,
          main=main,
          lty=ltps[1],
          col=cols[1],
          ...
          )
      if(length(x@chps)>1){        
        for(i in 2:dim(X)[2]){
          lines(X[,i], Y[,i], lty=ltps[i], col=cols[i], ...)
        }
      }
      if(zero.line) abline(h=0, col="grey")
      par(old)
      tryCatch(
          {
            old <- par(mar=c(5,0,4,1)+0.1)
            plot.new()
            legend(
                legend.position,
                legend=x@chps,
                col=cols,
                lty=ltps,
                bty=legend.bty,
                ...
                )
          },
          finally={
            par(old)
          }
          )
      layout(1)
      close(x)
      
      invisible(x)
      
    }
    )

## Define a boxplot method for SNPMaP objects
## Plots intensity/RAS for each chip, PM only
## Arguments: SNPMaP object, [function, character, character, ...]
## Returns: SNPMaP object invisibly
setGeneric("boxplot")

setMethod("boxplot",
    signature(x="SNPMaP"),
    function(
        x,
        FUN=function(x){x},
        ylab="guess",
        main="",
        ...
        )
    {
      if(x@table=="empty") stop(
            "this SNPMaP object is of 'empty' format"
            )
      x<-.checkOpen(x)
      
      if(ylab=="guess"){
        ylab<-switch(
            x@table,
            "raw"=,
            "long"=,
            "short"="Intensity",
            "ras"=,
            "rasS"="Relative Allele Score (RAS)"
            )
      }
      
      y<-list(
          stats=matrix(nrow=5, ncol=length(x@chps)),
          n=numeric(length(x@chps)),
          conf=matrix(nrow=2, ncol=length(x@chps)),
          out=numeric(0),
          group=numeric(0),
          names=x@chps
          )
      
      w<-switch(
          x@table,
          "raw"=,
          "long"=,
          "rasS"=1,
          "ras"=x@width,
          "short"=x@width*2
          )
      
      for(i in 1:length(x@chps)){
        ## All the columns corresponding to one chip
        z<-boxplot.stats(FUN(x@snpdata[,(1+((i-1)*w)):(i*w)]))
        y$stats[,i]<-z$stats
        y$n[i]<-z$n
        y$conf[,i]<-z$conf
        y$out<-append(y$out, z$out)
        y$group<-append(y$group, rep(i, length(z$out)))
      }
      
      bxp(y, ylab=ylab, main=main, ...)
      
      close(x)
      
      invisible(x)
    }
    )

## Define an image method to produce images of the raw chip data
## Arguments: SNPMaP object, [character or numeric, logical, function,
## character, numeric, numeric, numeric, ...]
## Returns: SNPMaP object invisibly
setGeneric("image")

setMethod("image",
    signature(x="SNPMaP"),
    function(
        x,
        chips=x@chps,
        prompt=FALSE,
        FUN=log,
        col=grey(seq(0,1,0.01)),
        fastRender=4,
        rows=x@celDim[1],
        cols=x@celDim[2],
        ...
        )
    {    
      if(is.character(chips))chips<-match(chips, x@chps)
      if(x@table!="raw") stop(
            "the image() method can only produce chip images from data in the 'raw' format"
            )
      
      x<-.checkOpen(x)
      
      old<-par(mar=c(2,2,4,2)+0.1)
      
      j<-0
      
      for(i in chips){
        ## Get the raw chip data formatted as a matrix    
        z<-matrix(
            FUN(rev(x@snpdata[,i])),
            rows,
            cols,
            byrow=TRUE
            )
        ## Select one in 'fastRender' rows and columns
        rowN<-(row(z)%%fastRender==0)[,1]
        colN<-(col(z)%%fastRender==0)[1,]
        
        image(
            1:sum(rowN),
            1:sum(colN),
            z[rowN, colN],
            col=col,
            axes=FALSE,
            main=x@chps[i],
            xlab="",
            ylab="",
            ...
            )
        
        j<-j+1
        
        if(prompt==TRUE & j<length(chips)){
          answer<-readline(prompt = paste(
                  "press <enter> to advance to image",
                  j+1,
                  "..."
                  )
              )
          if(answer=="\n"){next}
        }
      }
      
      par(old)
      
      close(x)
      
      invisible(x)
    }
    )

## 'open' method
## Arguments: SNPMaP object. Returns: SNPMaP object
setGeneric("open")

setMethod("open",
    signature(con="SNPMaP"),
    function(con, ...){
      con<-.checkOpen(con)
      invisible(con)
    }
    )

## 'close' method
## Arguments: SNPMaP object. Returns: SNPMaP object
setGeneric("close")

setMethod("close",
    signature(con="SNPMaP"),
    function(con, ...){
      if(is(con@snpdata, "FileDoubleMatrix")){
        if(isOpen(con@snpdata)){
          close(con@snpdata)
          !isOpen(con@snpdata)||stop(
              "cannot close connection to ",
              con@snpdata$pathname
              )
        }
        if(is(con@mismatch, "FileDoubleMatrix")){
          if(isOpen(con@mismatch)){
            close(con@mismatch)
            !isOpen(con@mismatch)||stop(
                "cannot close connection to ",
                con@mismatch$pathname
                )
          }
        }
      }
      invisible(con)
    }
    )

## as.matrix method
## Export snpdata scores to a standard named R matrix for further analysis
## Arguments: SNPMaP object, [logical]. Returns: matrix
setGeneric("as.matrix")

setMethod("as.matrix",
    signature(x="SNPMaP"),
    function(x, mm=FALSE, ...){
      
      d<-numeric(0)        
      
      if(mm){
        d<-dim(x@mismatch)
      }else{
        d<-dim(x@snpdata)
      }
      
      x<-.checkOpen(x)
      
      mat<-.gws(d[1], d[2])
      
      if(!is.matrix(mat)) stop(
            "not enough memory to convert to matrix"
            )
      
      if(mm){
        for(i in 1:d[2]){
          mat[,i]<-x@mismatch[,i]
        }
      }else{
        for(i in 1:d[2]){
          mat[,i]<-x@snpdata[,i]
        }
      }
      
      
      if(is(x@mismatch, "FileDoubleMatrix")){
        close(x@mismatch)
      }
      if(is(x@snpdata, "FileDoubleMatrix")){
        close(x@snpdata)
      }
      
      dimnames(mat)<-list(x@snps, x@cols)
      
      return(mat)
    }
    )
