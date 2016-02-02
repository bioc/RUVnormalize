%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  ./svdPlot.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{svdPlot}
\alias{svdPlot}


\title{Plot the data projected into the space spanned by their first two principal components}

\description{ The function takes as input a gene expression matrix
and plots the data projected into the space spanned by their first
two principal components.}

\arguments{
  \item{Y}{Expression matrix where the rows are the
  samples and the columns are the genes.}
  \item{annot}{A matrix
  containing the annotation to be plotted. Each row must correspond
  to a sample (row) of argument Y, each column must be a
  categorical or continuous descriptor for the sample. Optional.}
  \item{labels}{A vector with one element per sample (row) of
  argument Y. If this argument is specified, each sample is
  represented by its label. Otherwise, it is represented by a dot
  (if no annotation is provided) or by the value of the
  annotation. Optional.}
  \item{svdRes}{A list containing the result of svd(Y), possibly
  restricted to the first few singular values. Optional: if not
  provided, the function computes the SVD.}
  \item{plAnnots}{A list specifiying whether each column of the
  annot argument corresponds to a categorical or continuous
  factor. Each element of the list is named after a column of
  annot, and contains a string 'categorical' or 'continuous'. For
  each element of this list, a plot is produced where the samples
  are represented by colors corresponding to their annotation. Optional.}
  \item{kColors}{A vector of colors to be used to represent
  categorical factors. Optional: a default value is provided. If a
  categorical factors has more levels than the number of colors
  provided, colors are not used and the factor is represented in
  black.}
  \item{file}{A string giving the path to a pdf file for the plot. Optional.}
}

\value{
 A \code{\link[base]{list}} containing the result of svd(Y, nu=2, nv=0).
 }

\usage{
svdPlot(Y, annot=NULL, labels=NULL, svdRes=NULL, plAnnots=NULL, kColors=NULL, file=NULL)
}

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
    nu.coeff <- 1e-2
    
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

