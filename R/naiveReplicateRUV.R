
#########################################################################/**
## @RdocFunction naiveReplicateRUV
##
## @title "Remove unwanted variation from a gene expression matrix using replicate samples"
##
## \description{ The function takes as input a gene expression matrix
## as well as the index of negative control genes and replicate
## samples. It estimates and remove unwanted variation from the gene
## expression.}
##
## @synopsis
##
## \arguments{
##   \item{Y}{Expression matrix where the rows are the samples and the columns are the genes.}
##   \item{cIdx}{Column index of the negative control genes in Y, for estimation of unwanted variation.}
##   \item{scIdx}{Matrix giving the set of replicates. Each row is a
## set of arrays corresponding to replicates of the same sample. The
## number of columns is the size of the largest set of replicates, and
## the smaller sets are padded with -1 values. For example if the sets
## of replicates are (1,11,21), (2,3), (4,5), (6,7,8), the scIdx should
## be
## 1 11 21
## 2 3  -1
## 4 5  -1
## 6 7   8}
##   \item{k}{Desired rank for the estimated unwanted variation
##   term. The returned rank may be lower if the replicate arrays and
##   control genes did not contain a signal of rank k.}
##   \item{rrem}{Optional, indicates which arrays should be removed
## from the returned result. Useful if the replicate arrays were not
## actual samples but mixtures of RNA which are only useful to
## estimate UV but which should not be included in the analysis.}
##   \item{p}{Optional. If given, the function returns an estimate of the term of interest
##   using rank-p restriction of the SVD of the corrected matrix.}
##   \item{tol}{Directions of variance lower than this value in the
##   replicate samples are dropped (which may result in an estimated
##   unwanted variation term of rank smaller than k).}}
##
## \value{A @list containing the following terms:
##  \item{X, b}{if p is not NULL, contains an estimate of the factor
##  of interest (X) and its effect (beta) obtained using rank-p
##  restriction of the SVD of Y - W alpha.}
##  \item{W, a}{Estimates of the unwanted variation factors (W) and their effect (alpha).}
##  \item{cY}{The corrected expression matrix Y - W alpha}
##  \item{Yctls}{The differences of replicate arrays which were used to
##  estimate W and alpha.}
##  }
##
## @examples "../inst/extdata/naiveReplicateRUV.Rex"
##
##*/########################################################################


naiveReplicateRUV <- function(Y, cIdx, scIdx, k, rrem=NULL, p=NULL, tol=1e-6)
{

  
  ##----------------
  ## Initialization
  ##----------------
  
  m <- nrow(Y)
  n <- ncol(Y)    

  ##-----------------------
  ## Build control samples
  ##-----------------------

  scIdx <- scIdx[rowSums(scIdx > 0) >= 2,,drop=FALSE]
  
  Yctls <- matrix(0,prod(dim(scIdx)),ncol(Y))
  c <- 0
  
  for(ii in 1:nrow(scIdx)){
    for(jj in 1:(ncol(scIdx) -1)){      
      if(scIdx[ii,jj] == -1)
        next
      c <- c+1
      Yctls[c,] <- Y[scIdx[ii,jj],,drop=FALSE] - colMeans(Y[scIdx[ii,((jj != 1:ncol(scIdx)) & (scIdx[ii,] > 0))],,drop=FALSE])
    }
  }
  Yctls <- Yctls[rowSums(Yctls) != 0,]  
  
  ## Remove replicates (optional)
  if(!is.null(rrem)){
    nremoved <- length(rrem)
    Y <- Y[-rrem,]
    m <- m - nremoved
  }
  
  Y <- rbind(Y, Yctls)
  sctl <- (m+1):(m+nrow(Yctls))
  ##return(Yctls)
  
  ##------------
  ## Estimation
  ##------------
  
    ## Get alphas from control samples
    svdRes <- svd(Y[sctl,], nu=0, nv=k)
    
    k.max <- sum(svdRes$d > tol)
    if(k > k.max){
        warning(sprintf('k larger than the rank of the control sample matrix Y[sctl, ]. Using k=%d instead', k.max))
        k <- k.max
    }
  ## k <- min(k, max(which(svdRes$d > tol))) # Don't return directions with 0 variance
  a <- t(as.matrix(svdRes$v[, 1:k]))

  ## Get W from alphas and control genes (by regression)
  W <- Y[, cIdx] %*% t(solve(a[, cIdx, drop=FALSE] %*% t(a[, cIdx, drop=FALSE]), a[, cIdx, drop=FALSE]))
  
  ## Get Xb from Wa (optional)

  if(!is.null(p)){
      Wa <- W %*% a     
      svdRes <- svd(Y[1:m,-cIdx] - Wa[1:m,-cIdx], nu=p, nv=p)
      X <- as.matrix(svdRes$u)
      b <- matrix(0, p, n)    
      b[,-cIdx] <- diag(svdRes$d[1:p],nrow=p) %*% t(as.matrix(svdRes$v))
  }
  else{    
    X <- b <- NULL
  }

  cY <- Y[1:m,] - W[1:m,] %*% a
  
  return(list(X=X, b=b, W=W, a=a, cY=cY, Yctls=Yctls))
}
