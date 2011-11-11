## extrA.R
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


`extrA` <- function(o, allele='A'){   # would be more elegant, poss faster to
                           # redimm array but may not play nice with
                           # filedoublematrix
  # extract allele A or B columns from a short matrix
  open(o)
  a <- 1:length(o@cols)
  j<-rep(c(rep(TRUE,o@width),rep(FALSE,o@width)), length(o@chps))
  if (allele=='B'){ j<- !j }
  A <-o@snpdata[,a[j]]
  close(o)
  A
}



`suco` <- function (m, chps) {
# sum columns by chip. assumes columns relevant to each chip are adjacent
# and equal in number. initial version takes & gives matrix such as
# extrA gives.  note that it only uses chps for the length --
  nchips      <- length(chps)
  width       <- dim(m)[2]
  nsnps       <- dim(m)[1]
  stopifnot ( width%%nchips==0 )
  valsperchip <- width/nchips
  m <- t(m)
  dim(m) <- c(valsperchip,nsnps * nchips) 
  m      <- colSums(m)
  matrix(m,ncol=nchips,byrow=TRUE)
}

