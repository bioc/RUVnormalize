%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  ./naiveRandRUV.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{naiveRandRUV}
\alias{naiveRandRUV}


\title{Remove unwanted variation from a gene expression matrix using negative control genes}

\description{ The function takes as input a gene expression matrix
as well as the index of negative control genes. It estimates
unwanted variation from these control genes, and removes them by
regression, using ridge and/or rank regularization.}

\usage{
naiveRandRUV(Y, cIdx, nuCoeff=0.001, k=nrow(Y))
}

\arguments{
  \item{Y}{Expression matrix where the rows are the samples and the columns are the genes.}
  \item{cIdx}{Column index of the negative control genes in Y, for estimation of unwanted variation.}
  \item{nuCoeff}{Regularization parameter for the unwanted variation.}
  \item{k}{Desired rank for the estimated unwanted variation term.}
}

\value{ A \code{\link[base]{matrix}} corresponding to the gene expression after
 substraction of the estimated unwanted variation term.  }

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
    nsY <- naiveRandRUV(Y, cIdx, nuCoeff=0, k=k)

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
    nu.coeff <- 1e-2

    ## Correction
    nY <- naiveRandRUV(Y, cIdx, nuCoeff=nu.coeff, k=k)

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

