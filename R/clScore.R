#########################################################################/**
## @RdocFunction clScore
##
## @title "Computes a distance between two partitions of the same data"
##
## \description{ The function takes as input two partitions of a
## dataset into clusters, and returns a number which is small if the
## two partitions are close, large otherwise.}
##
## @synopsis
##
## \arguments{
##   \item{c1}{A @vector giving the assignment of the samples to cluster for the first partition}
##   \item{c2}{A @vector giving the assignment of the samples to cluster for the second partition}
## }
##
## \value{A number corresponding to the distance between c1 and c2}
##
## @examples "../inst/extdata/naiveRandRUV.Rex"
##
##*/########################################################################

clScore <- function(c1,c2){
  uc1 <- unique(c1)
  uc2 <- unique(c2)
  tmp <- 0
  for(c in uc1){
    for(cc in uc2){
      tmp <- tmp + (sum((c1 == c) & (c2 == cc))^2)/(sum(c1 == c)*sum(c2 == cc))
    }
  }
  return((length(uc1)+length(uc2))/2 - tmp)
}
