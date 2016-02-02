%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  ./clScore.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{clScore}
\alias{clScore}


\title{Computes a distance between two partitions of the same data}

\description{ The function takes as input two partitions of a
dataset into clusters, and returns a number which is small if the
two partitions are close, large otherwise.}

\usage{
clScore(c1, c2)
}

\arguments{
  \item{c1}{A \code{\link[base]{vector}} giving the assignment of the samples to cluster for the first partition}
  \item{c2}{A \code{\link[base]{vector}} giving the assignment of the samples to cluster for the second partition}
}

\value{A number corresponding to the distance between c1 and c2}

\examples{
if(require('RUVnormalizeData')){

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

    ##--------------------------
    ## Naive RUV-2 no shrinkage
    ##--------------------------

    k <- 20
    nu <- 0

    ## Correction
    nsY <- naiveRandRUV(Y, cIdx, nu.coeff=0, k=k)

    ## Clustering of the corrected data
    sdY <- apply(nsY, 2, sd)
    ssd <- sort(sdY,decreasing=TRUE,index.return=TRUE)$ix
    kmres2ns <- kmeans(nsY[,ssd[1:nKeep],drop=FALSE],centers=2,nstart=200)
    vclust2ns <- kmres2ns$cluster
    nsScore <- clScore(vclust2ns, X)

    ## Plot of the corrected data
    svdRes2ns <- NULL
    svdRes2ns <- svdPlot(nsY[, ssd[1:nKeep], drop=FALSE],
                         annot=annot,
                         labels=lab.and.region,
                         svdRes=svdRes2ns,
                         plAnnots=plAnnots,                    
                         kColors=gender.col, file=NULL)   

    ##--------------------------
    ## Naive RUV-2 + shrinkage
    ##--------------------------

    k <- m
    nu.coeff <- 1e-3

    ## Correction
    nY <- naiveRandRUV(Y, cIdx, nu.coeff=nu.coeff, k=k)

    ## Clustering of the corrected data
    sdY <- apply(nY, 2, sd)
    ssd <- sort(sdY,decreasing=TRUE,index.return=TRUE)$ix
    kmres2 <- kmeans(nY[,ssd[1:nKeep],drop=FALSE],centers=2,nstart=200)
    vclust2 <- kmres2$cluster
    nScore <- clScore(vclust2,X)

    ## Plot of the corrected data
    svdRes2 <- NULL
    svdRes2 <- svdPlot(nY[, ssd[1:nKeep], drop=FALSE],
                       annot=annot,
                       labels=lab.and.region,
                       svdRes=svdRes2,
                       plAnnots=plAnnots,                    
                       kColors=gender.col, file=NULL)   
}
}

