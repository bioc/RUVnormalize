if(require('RUVnormalizeData') && require('spams')){
    ## Load the library
    library(RUVnormalize)
    library(spams)

    ## Load the data
    data('gender', package='RUVnormalizeData')

    ## Extract regions and labs for plotting purposes
    lregions <- sapply(rownames(Y),FUN=function(s) strsplit(s,'_')[[1]][2])
    llabs <- sapply(rownames(Y),FUN=function(s) strsplit(s,'_')[[1]][3])

    ## Dimension of the factors
    m <- nrow(Y)
    n <- ncol(Y)
    p <- ncol(X)

    Y <- scale(Y, scale=FALSE) # Center gene expressions

    cIdx <- which(isrlctl) # Negative control genes

    ## Prepare control samples

    scIdx <- matrix(-1,84,3)
    rny <- rownames(Y)
    added <- c()
    c <- 0

                                        # Replicates by lab
    for(r in 1:(length(rny) - 1)){
        if(r %in% added)
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
        if(r %in% added)
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

    ## Remove platform effect by centering.

    Y[Z[,4]==0,] <- scale(Y[Z[,4]==0,], scale=FALSE)
    Y[Z[,4]==1,] <- scale(Y[Z[,4]==1,], scale=FALSE)

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
    paramXb$lambda <- 8e-2
    paramXb$lambda2 <- 0

    ## Correction
    iRes <- iterativeRUV(Y, cIdx, scIdx=NULL, paramXb, k=nrow(Y), nu.coeff=1e-3/2,
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
