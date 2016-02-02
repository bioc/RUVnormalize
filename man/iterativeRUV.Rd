\name{iterativeRUV}
\alias{iterativeRUV}


 \title{Remove unwanted variation from a gene expression matrix
 using control genes, optionally replicate samples, and iterative
 estimates of the factor of interest}

 \description{ The function takes as input a gene expression matrix
 as well as the index of negative control genes and replicate
 samples. It estimates and remove unwanted variation from the gene
 expression. The major difference with naiveRandRUV and
 naiveReplicateRUV is that iterativeRUV jointly estimates the
 factor of interest and the unwanted variation term. It does so
 iteratively, by estimating each term using the current estimate of
 the other one.}

 \usage{
iterativeRUV(Y, cIdx, scIdx=NULL, paramXb, k, nu.coeff=0, cEps=1e-08, maxIter=30,
  Wmethod="svd", Winit=NULL, wUpdate=maxIter + 1, tol=1e-6)
}

 \arguments{
   \item{Y}{Expression matrix where the rows are the samples and the columns are the genes.}
   \item{cIdx}{Column index of the negative control genes in Y, for estimation of unwanted variation.}
   \item{scIdx}{Matrix giving the set of replicates. Each row is a
 set of arrays corresponding to replicates of the same sample. The
 number of columns is the size of the largest set of replicates, and
 the smaller sets are padded with -1 values. For example if the sets
 of replicates are (1,11,21), (2,3), (4,5), (6,7,8), the scIdx should
 be
 1 11 21
 2 3  -1
 4 5  -1
 6 7   8}
   \item{paramXb}{A \code{\link[base]{list}} containing parameters for the estimation of the term of interest:
                  K corresponds to the rank of X.
                  lambda is the regularization parameter. Large values of lambda lead to sparser, more shrunk estimates of beta.
                  D, batch, iter and mode should not be modified unless you are familiar with sparse dictionary learning algorithms.}

   \item{k}{Desired rank for the estimated unwanted variation
   term. The returned rank may be lower if the replicate arrays and
   control genes did not contain a signal of rank k.}
   \item{nu.coeff}{Regularization parameter for the unwanted variation.}
   \item{cEps}{tolerance for relative changes of Wa and Xb estimators at
               each step. When both get smaller than cEps, the iterations stop.}
   \item{maxIter}{Maximum number of iterations.}
   \item{Wmethod}{'svd' or 'rep', depending whether W is estimated from
 control genes or replicate samples.}
   \item{Winit}{Optionally provides an initial value for W.}
   \item{wUpdate}{Number of iterations between two updates of W. By default,
 W is never updated. Make sure that enough iterations are done after
 the last update of W. E.g, setting W to maxIter will only allow for
 one iteration of estimating alpha given (Xb, W) and no
 re-estimation of Xb.}
   \item{tol}{Smallest ratio allowed between a squared singular
 value of Y[, cIdx] and the largest of these squared singular
 values. All smaller singular values are discarded.}
  }

\details{In terms of model, the rank k can be thought of as the
number of independent sources of unwanted variation in the data (i.e.,
if one source is a linear combination of other sources, it does not
increase the rank). The ridge nu.coeff should be inversely proportional
to the (expected) magnitude of the unwanted variation.

In practice, even if the real number of independent sources of unwanted
variation (resp. their magnitude) is known, using a smaller k (resp.,
larger ridge) could yield better corrections because one may not have
enough samples to effectively estimate all the effects.

More intuition and guidance on the practical choice of these parameters are available in the paper
(\url{http://biostatistics.oxfordjournals.org/content/17/1/16.full}) and its
supplement
(\url{http://biostatistics.oxfordjournals.org/content/suppl/2015/08/17/kxv026.DC1/kxv026supp.pdf}). In
particular:
- Equation 2.3 in the manuscript gives an interpretation of the
ridge parameter in terms of a probabilistic model.
- Section 5.1 of the manuscript provides guidelines to select both
parameters on real data.
- Section 3 of the supplement compares the effect of reducing the
  rank and increasing the ridge.
- Section 4 of the supplement gives a detailed discussion of how to
  select the ridge parameter on a real example.
}
  
 \value{
  A \code{\link[base]{list}} containing the following terms:
  \item{X, b}{if p is not NULL, contains an estimate of the factor of interest (X) and its effect (beta) obtained using rank-p restriction of the SVD of Y - W  alpha.}
  \item{W, a}{Estimates of the unwanted variation factors (W) and their effect (alpha).}
  \item{cY}{The corrected expression matrix Y - W  alpha.}
  }


  \examples{
if(require('RUVnormalizeData') && require('spams')){
    ## Load the spams library
    library(spams)
    
    ## Load the data
    data('gender', package='RUVnormalizeData')
    
    Y <- t(exprs(gender))
    X <- as.numeric(phenoData(gender)$gender == 'M')
    X <- X - mean(X)
    X <- cbind(X/(sqrt(sum(X^2))))
    chip <- annotation(gender)
    
    ## Extract regions and labs for plotting purposes
    lregions <- sapply(rownames(Y),FUN=function(s) strsplit(s,'_')[[1]][2])
    llabs <- sapply(rownames(Y),FUN=function(s) strsplit(s,'_')[[1]][3])
    
    ## Dimension of the factors
    m <- nrow(Y)
    n <- ncol(Y)
    p <- ncol(X)
    
    Y <- scale(Y, scale=FALSE) # Center gene expressions
    
    cIdx <- which(featureData(gender)$isNegativeControl) # Negative control genes
    
    ## Prepare plots
    annot <- cbind(as.character(sign(X)))
    colnames(annot) <- 'gender'
    plAnnots <- list('gender'='categorical')
    lab.and.region <- apply(rbind(lregions, llabs),2,FUN=function(v) paste(v,collapse='_'))
    gender.col <- c('-1' = "deeppink3", '1' = "blue")
    
    ## Remove platform effect by centering.
    
    Y[chip=='hgu95a.db',] <- scale(Y[chip=='hgu95a.db',], scale=FALSE)
    Y[chip=='hgu95av2.db',] <- scale(Y[chip=='hgu95av2.db',], scale=FALSE)
    
    ## Number of genes kept for clustering, based on their variance
    nKeep <- 1260
    
    ## Prepare control samples
    
    scIdx <- matrix(-1,84,3)
    rny <- rownames(Y)
    added <- c()
    c <- 0
    
    ## Replicates by lab
    for(r in 1:(length(rny) - 1)){
        if(r \%in\% added)
            next
        c <- c+1
        scIdx[c,1] <- r
        cc <- 2
        for(rr in seq(along=rny[(r+1):length(rny)])){
            if(all(strsplit(rny[r],'_')[[1]][-3] ==  strsplit(rny[r+rr],'_')[[1]][-3])){
                scIdx[c,cc] <- r+rr
                cc <- cc+1
                added <- c(added,r+rr)
            }
        }   
    }
    scIdxLab <- scIdx
    
    scIdx <- matrix(-1,84,3)
    rny <- rownames(Y)
    added <- c()
    c <- 0
    
    ## Replicates by region
    for(r in 1:(length(rny) - 1)){
        if(r \%in\% added)
            next
        c <- c+1
        scIdx[c,1] <- r
        cc <- 2
        for(rr in seq(along=rny[(r+1):length(rny)])){
            if(all(strsplit(rny[r],'_')[[1]][-2] ==  strsplit(rny[r+rr],'_')[[1]][-2])){
                scIdx[c,cc] <- r+rr
                cc <- cc+1
                added <- c(added,r+rr)
            }
        }
    }
    scIdx <- rbind(scIdxLab,scIdx)
    
    ## Number of genes kept for clustering, based on their variance
    nKeep <- 1260
    
    ## Prepare plots
    annot <- cbind(as.character(sign(X)))
    colnames(annot) <- 'gender'
    plAnnots <- list('gender'='categorical')
    lab.and.region <- apply(rbind(lregions, llabs),2,FUN=function(v) paste(v,collapse='_'))
    gender.col <- c('-1' = "deeppink3", '1' = "blue")
    
    ##---------------------------
    ## Iterative replicate-based
    ##---------------------------
    
    cEps <- 1e-6
    maxIter <- 30
    p <- 20
    
    paramXb <- list()
    paramXb$K <- p
    paramXb$D <- matrix(c(0.),nrow = 0,ncol=0)
    paramXb$batch <- TRUE
    paramXb$iter <- 1
    paramXb$mode <- 'PENALTY'
    paramXb$lambda <- 0.25
    
    ## Correction
    iRes <- iterativeRUV(Y, cIdx, scIdx, paramXb, k=20, nu.coeff=0,
                         cEps, maxIter,
                         Wmethod='rep', wUpdate=11)
    
    ucY <- iRes$cY
    
    ## Cluster the corrected data
    sdY <- apply(ucY, 2, sd)
    ssd <- sort(sdY,decreasing=TRUE,index.return=TRUE)$ix
    kmresIter <- kmeans(ucY[,ssd[1:nKeep]],centers=2,nstart=200)
    vclustIter <- kmresIter$cluster
    IterScore <- clScore(vclustIter,X)
    
    ## Plot the corrected data
    svdResIter <- NULL
    svdResIter <- svdPlot(ucY[, ssd[1:nKeep], drop=FALSE],
                          annot=annot,
                          labels=lab.and.region,
                          svdRes=svdResIter,
                          plAnnots=plAnnots,                    
                          kColors=gender.col, file=NULL)   
    
    ##--------------------------
    ## Iterated ridge
    ##--------------------------
    
    paramXb <- list()
    paramXb$K <- p
    paramXb$D <- matrix(c(0.),nrow = 0,ncol=0)
    paramXb$batch <- TRUE
    paramXb$iter <- 1
    paramXb$mode <- 'PENALTY' #2
    paramXb$lambda <- 1
    paramXb$lambda2 <- 0
    
    ## Correction
    iRes <- iterativeRUV(Y, cIdx, scIdx=NULL, paramXb, k=nrow(Y), nu.coeff=1e-2/2,
                         cEps, maxIter,
                         Wmethod='svd', wUpdate=11)
    
    nrcY <- iRes$cY
    
    ## Cluster the corrected data
    sdY <- apply(nrcY, 2, sd)
    ssd <- sort(sdY,decreasing=TRUE,index.return=TRUE)$ix
    kmresIter <- kmeans(nrcY[,ssd[1:nKeep]],centers=2,nstart=200)
    vclustIter <- kmresIter$cluster
    IterRandScore <- clScore(vclustIter,X)
    
    ## Plot the corrected data
    svdResIterRand <- NULL
    svdResIterRand <- svdPlot(nrcY[, ssd[1:nKeep], drop=FALSE],
                              annot=annot,
                              labels=lab.and.region,
                              svdRes=svdResIterRand,
                              plAnnots=plAnnots,                    
                              kColors=gender.col, file=NULL)   
}
}

