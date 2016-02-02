
#########################################################################/**
## @RdocFunction naiveRandRUV
##
## @title "Remove unwanted variation from a gene expression matrix using negative control genes"
##
## \description{ The function takes as input a gene expression matrix
## as well as the index of negative control genes. It estimates
## unwanted variation from these control genes, and removes them by
## regression, using ridge and/or rank regularization.}
##
## @synopsis
##
## \arguments{
##   \item{Y}{Expression matrix where the rows are the samples and the columns are the genes.}
##   \item{cIdx}{Column index of the negative control genes in Y, for estimation of unwanted variation.}
##   \item{nu.coeff}{Regularization parameter for the unwanted variation.}
##   \item{k}{Desired rank for the estimated unwanted variation term.}
##   \item{tol}{Smallest ratio allowed between a squared singular
## value of Y[, cIdx] and the largest of these squared singular
## values. All smaller singular values are discarded.}
## }
##
## \value{ A @matrix corresponding to the gene expression after
##  substraction of the estimated unwanted variation term.  }
##
## @examples "../inst/extdata/naiveRandRUV.Rex"
##
##*/########################################################################

naiveRandRUV <- function(Y, cIdx, nu.coeff=1e-3, k=min(nrow(Y), length(cIdx)), tol=1e-6){
    
    ## W is the square root of the empirical covariance on the control
    ## genes.
  
    svdYc <- svd(Y[, cIdx, drop=FALSE])
    k.max <- sum(svdYc$d^2/svdYc$d[1]^2 > tol)
    if(k > k.max){
        warning(sprintf('k larger than the rank of Y[, cIdx]. Using k=%d instead', k.max))
        k <- k.max
    }
    W <- svdYc$u[, 1:k, drop=FALSE] %*% diag(svdYc$d[1:k], nrow=k)
    
    ## Regularization heuristic: nu is a fraction of the largest eigenvalue of WW'
    
    nu <- nu.coeff*svdYc$d[1]^2 #/ (length(cIdx)+1)
    
    ## Naive correction: ridge regression of Y against W
    
    nY <- Y - W %*% solve(t(W)%*%W + nu*diag(k), t(W) %*% Y)
    
    return(nY)
}
