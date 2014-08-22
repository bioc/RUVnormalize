#########################################################################/**
## @RdocFunction iterativeRUV
##
## @title "Remove unwanted variation from a gene expression matrix
## using control genes, optionally replicate samples, and iterative
## estimates of the factor of interest"
##
## \description{ The function takes as input a gene expression matrix
## as well as the index of negative control genes and replicate
## samples. It estimates and remove unwanted variation from the gene
## expression. The major difference with naiveRandRUV and
## naiveReplicateRUV is that iterativeRUV jointly estimates the
## factor of interest and the unwanted variation term. It does so
## iteratively, by estimating each term using the current estimate of
## the other one.}
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
##   \item{paramXb}{A @list containing parameters for the estimation of the term of interest:
##                  K corresponds to the rank of X.
##                  lambda is the regularization parameter. Large values of lambda lead to sparser, more shrunk estimates of beta.
##                  D, batch, iter and mode should not be modified unless you are familiar with sparse dictionary learning algorithms.}
## 
##   \item{k}{Desired rank for the estimated unwanted variation
##   term. The returned rank may be lower if the replicate arrays and
##   control genes did not contain a signal of rank k.}
##   \item{nu.coeff}{Regularization parameter for the unwanted variation.}
##   \item{cEps}{tolerance for relative changes of Wa and Xb estimators at
##               each step. When both get smaller than cEps, the iterations stop.}
##   \item{maxIter}{Maximum number of iterations.}
##   \item{Wmethod}{'svd' or 'rep', depending whether W is estimated from
## control genes or replicate samples.}
##   \item{Winit}{Optionally provides an initial value for W.}
##   \item{wUpdate}{Number of iterations between two updates of W. By default,
## W is never updated. Make sure that enough iterations are done after
## the last update of W. E.g, setting W to maxIter will only allow for
## one iteration of estimating alpha given (Xb, W) and no
## re-estimation of Xb.}
##  }
##
## \value{
##  A @list containing the following terms:
##  \item{X, b}{if p is not NULL, contains an estimate of the factor of interest (X) and its effect (beta) obtained using rank-p restriction of the SVD of Y - W  alpha.}
##  \item{W, a}{Estimates of the unwanted variation factors (W) and their effect (alpha).}
##  \item{cY}{The corrected expression matrix Y - W  alpha.}
##  }
##
##
## @examples "../inst/extdata/iterativeRUV.Rex"
##
##*/########################################################################

iterativeRUV <- function(Y, cIdx, scIdx=NULL, paramXb, k, nu.coeff=0,
                         cEps=1e-8, maxIter=30,
                         Wmethod='svd', Winit=NULL, wUpdate=maxIter+1){

    if (!require(spams)) {
        stop('This function requires the spams package to be loaded. Spams can be obtained for free from http://spams-devel.gforge.inria.fr/downloads.html')
    }
    
  m <- nrow(Y)
  n <- ncol(Y)
  p <- paramXb$K

  if(paramXb$mode == 'kmeans'){
    paramXb$lambda <- 0
  }
  if(is.null(paramXb$lambda2)){
    paramXb$lambda2 <- 0
  }
  
  ##--------------------
  ## Optimize
  ##--------------------

  converged <- 0
  iter <- 0

  aRep <- NULL
  
  if(!is.null(Winit)){ # W is provided: estimate alpha only
    W <- Winit
    nu <- nu.coeff * svd(W)$d[1]^2
    a <- solve(t(W)%*%W + 2*nu*diag(ncol(W)), t(W) %*% Y)
  }
  else{
    if(Wmethod == 'rep'){ # Init by 1-shot, using replicates
      ## sRes <- shotRUVcs(Y, cIdx, scIdx, k, rrem=NULL, p=NULL)
      sRes <- naiveReplicateRUV(Y, cIdx, scIdx, k, rrem=NULL, p=NULL)
      W <- sRes$W[1:m, ]
      aRep <- sRes$a
      Yctls <- sRes$Yctls
      aRepProj <- t(solve(aRep %*% t(aRep), aRep))
      a <- aRep
      nu <- nu.coeff*svd(W)$d[1]^2
    }
    else if(Wmethod == 'svd'){ # Use random alpha model with control genes only
      svdYc <- svd(Y[, cIdx])
      W <- svdYc$u[,1:k] %*% diag(svdYc$d[1:k])
      nu <- nu.coeff * svdYc$d[1]^2
      a <- solve(t(W)%*%W + 2*nu*diag(k), t(W) %*% Y)
    }else{
      stop('Unknown method to evaluate W. Choose between svd and rep')
    }
  }

  Wa <- W %*% a

  X <- matrix(0,m,p)
  b <- matrix(0,p,n)
  Xb <- X %*% b

  while(!converged){

    iter <- iter + 1
    
    ## (X,b) from (W,a)    
    print('Estimating (X,b) from (W,a)')
    XbOld <- Xb

    if(paramXb$mode == 'kmeans'){
      kmres <- kmeans((Y[, -cIdx] - Wa[1:m, -cIdx]),centers=p,nstart=20)
      idx <- kmres$cluster
      for(kk in 1:p){
        X[, kk] <- cbind(as.numeric(idx==kk))
      }
      b[, -cIdx] <- kmres$centers
    }
    else if(paramXb$mode == 'PENALTY'){
      ## Sparse dictionary (L1) version
      X <- spams::spams.trainDL(X=Y[, -cIdx] - Wa[1:m, -cIdx], D=paramXb$D,
                          lambda1=paramXb$lambda, lambda2=paramXb$lambda2,
                          mode=paramXb$mode, batch=paramXb$batch, iter=paramXb$iter, K=paramXb$K)
      paramXb$D =  X
      b[, -cIdx] <- as.matrix(spams::spams.lasso(X=Y[, -cIdx] - Wa[1:m, -cIdx], D=X,
                                           lambda1=paramXb$lambda, lambda2=paramXb$lambda2))
    }else if(paramXb$mode == 'PENALTY2'){
      ## Sparse dictionary (L0) version    
      X <- spams::spams.trainDL(X=Y[, -cIdx] - Wa[1:m, -cIdx], D=paramXb$D,
                         lambda1=paramXb$lambda, mode=paramXb$mode, batchsize=paramXb$batch, iter=paramXb$iter, K=paramXb$K)
      paramXb$D = X
      b[, -cIdx] <- as.matrix(spams::spams.omp(Y[, -cIdx] - Wa[1:m, -cIdx], D=X, L=paramXb$L, eps=paramXb$eps))      
    }else if(paramXb$mode == 'SVD'){
      svdYmWa <- svd(Y[, -cIdx] - Wa[1:m, -cIdx], nu=p, nv=p)
      X <- svdYmWa$u
      b[, -cIdx] <- diag(svdYmWa$d[1:p]/(2*paramXb$lambda2 + 1)) %*% t(svdYmWa$v)
    }
    else{
      stop(sprintf('Unknown mode %s',paramXb$mode))
    }
    
    Xb <- X %*% b
    
    WaOld <- Wa
    WOld <- W    

    ## Update W every wUpdate iterations
    if(iter / wUpdate == iter %/% wUpdate){
      
      print('')
      print('*************************************')
      print('Re-estimating W from Y - Xb residuals')
      print('*************************************')
      print('')
      
      if(Wmethod == 'svd'){
        svdYmXb <- svd((Y - Xb), nu=k, nv=0)
         W <- svdYmXb$u %*% diag(svdYmXb$d[1:k])
      }
      else if(Wmethod == 'rep'){        
        aRep <- t(svd(rbind((Y - Xb), Yctls), nu=0, nv=k)$v)
        ## Temporary fix: occasionally LAPACK will return singular
        ## vectors with non-unit norm.
        
        if(any(abs(rowSums(aRep^2) - 1) > 1e-8)){ # LAPACK bug
          print('LAPACK bug: svd returned wrong singular vectors, using LINPACK')
          warning('LAPACK bug: svd returned wrong singular vectors, using LINPACK')
          aRep <- t(svd(rbind((Y - Xb), Yctls), nu=0, nv=nrow(rbind((Y - Xb), Yctls)), LINPACK=TRUE)$v[, 1:k, drop=FALSE])
          if(any(abs(rowSums(aRep^2) - 1) > 1e-8)){
            stop('LAPACK bug: svd returned wrong singular vectors, LINPACK did not work either')
          }
        }
        aRepProj <- t(solve(aRep %*% t(aRep), aRep))
        W <- (Y - Xb) %*% aRepProj
      }
      nu <- nu.coeff * svd(W)$d[1]^2
    }

    ## a from (X,b)
    print('Estimating a from (X,b)')

    a <- solve(t(W)%*%W + 2*nu*diag(ncol(W)), t(W) %*% (Y - Xb))
    
    Wa <- W %*% a
    
    print('Update done')
    
    l2Err <- (norm((Y - Xb - Wa[1:m, ]), 'F')^2)/(m*n)
        
    dXb <- norm(Xb-XbOld,'F')/norm(Xb,'F')
    dWa <- norm(Wa-WaOld,'F')/norm(Wa,'F')
    
    if(ncol(W) != ncol(WOld)){
      dW <- Inf}
    else{
      dW <- norm(W-WOld,'F')/norm(W,'F')      
    }  
    
    print(sprintf('iter %d, dXb/Xb=%g, dWa/Wa=%g, dW/W=%g, l2Err=%g', iter, dXb, dWa, dW, l2Err))

    if(iter >= maxIter || (!is.nan(max(dXb,dWa)) && max(dXb,dWa) < cEps)){
      converged = 1
    }
  }

  cY <- Y - Wa[1:m, ]
  
  return(list(W=W, a=a, X=X, b=b, cY=cY))
}
