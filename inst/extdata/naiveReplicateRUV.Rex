if(require('RUVnormalizeData')){
    ## Load the library
    library(RUVnormalize)

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

    ## Correction
    sRes <- naiveReplicateRUV(Y, cIdx, scIdx, k=20)

    ## Clustering on the corrected data
    sdY <- apply(sRes$cY, 2, sd)
    ssd <- sort(sdY,decreasing=TRUE,index.return=TRUE)$ix
    kmresRep <- kmeans(sRes$cY[,ssd[1:nKeep],drop=FALSE],centers=2,nstart=200)
    vclustRep <- kmresRep$cluster
    RepScore <- clScore(vclustRep,X)

    ## Plot of the corrected data
    svdResRep <- NULL
    svdResRep <- svdPlot(sRes$cY[, ssd[1:nKeep], drop=FALSE], 
                         annot=annot, 
                         labels=lab.and.region, 
                         svdRes=svdResRep, 
                         plAnnots=plAnnots,                     
                         kColors=gender.col, file=NULL)    
}
