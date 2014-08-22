## Gender study data. Requires the RUVnormalizeData data package.

if(require('RUVnormalizeData')){
    ## library(spams)

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

    ## Prepare control samples

    scIdx <- matrix(-1,84,3)
    rny <- rownames(Y)
    added <- c()
    c <- 0

    ## Replicates by lab
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

    ##-----------------------
    ## Run different methods
    ##-----------------------

    ## No correction

    sdY <- apply(Y, 2, sd)
    ssd <- sort(sdY, decreasing=TRUE, index.return=TRUE)$ix

    kmres <- kmeans(Y[, ssd[1:nKeep], drop=FALSE], centers=2, nstart=200)
    vclust <- kmres$cluster
    uScore <- clScore(vclust,X)

    svdResUncorr <- NULL
    svdResUncorr <- svdPlot(Y[, ssd[1:nKeep], drop=FALSE], 
                            annot=annot, 
                            labels=lab.and.region, 
                            svdRes=svdResUncorr, 
                            plAnnots=plAnnots, 
                            kColors=gender.col, file=NULL)

    ## Centering by region-lab
    YmeanCorr <- Y
    for(rr in unique(lregions)){
        for(ll in unique(llabs)){
            YmeanCorr[(lregions==rr)&(llabs==ll),] <- scale(YmeanCorr[(lregions==rr)&(llabs==ll),],scale=FALSE)
        }
    }

    sdY <- apply(YmeanCorr, 2, sd)
    ssd <- sort(sdY,decreasing=TRUE,index.return=TRUE)$ix

    kmresMC <- kmeans(YmeanCorr[,ssd[1:nKeep],drop=FALSE],centers=2,nstart=200)
    vclustMC <- kmresMC$cluster
    MCScore <- clScore(vclustMC, X)

    svdResMC <- NULL
    svdResMC <- svdPlot(YmeanCorr[, ssd[1:nKeep], drop=FALSE], 
                        annot=annot, 
                        labels=lab.and.region, 
                        svdRes=svdResMC, 
                        plAnnots=plAnnots,                     
                        kColors=gender.col, file=NULL)    

    ## Naive RUV-2 no shrinkage

    k <- 20
    nu <- 0

    nsY <- naiveRandRUV(Y, cIdx, nuCoeff=0, k=k)

    sdY <- apply(nsY, 2, sd)
    ssd <- sort(sdY,decreasing=TRUE,index.return=TRUE)$ix

    kmres2ns <- kmeans(nsY[,ssd[1:nKeep],drop=FALSE],centers=2,nstart=200)
    vclust2ns <- kmres2ns$cluster
    nsScore <- clScore(vclust2ns, X)

    svdRes2ns <- NULL
    svdRes2ns <- svdPlot(nsY[, ssd[1:nKeep], drop=FALSE], 
                         annot=annot, 
                         labels=lab.and.region, 
                         svdRes=svdRes2ns, 
                         plAnnots=plAnnots,                     
                         kColors=gender.col, file=NULL)    

    ## Naive RUV-2 + shrinkage

    k <- m
    nu.coeff <- 1e-2

    nY <- naiveRandRUV(Y, cIdx, nuCoeff=nu.coeff, k=k)

    sdY <- apply(nY, 2, sd)
    ssd <- sort(sdY,decreasing=TRUE,index.return=TRUE)$ix

    kmres2 <- kmeans(nY[,ssd[1:nKeep],drop=FALSE],centers=2,nstart=200)
    vclust2 <- kmres2$cluster
    nScore <- clScore(vclust2,X)

    svdRes2 <- NULL
    svdRes2 <- svdPlot(nY[, ssd[1:nKeep], drop=FALSE], 
                       annot=annot, 
                       labels=lab.and.region, 
                       svdRes=svdRes2, 
                       plAnnots=plAnnots,                     
                       kColors=gender.col, file=NULL)    

    ## Replicate-based

    sRes <- naiveReplicateRUV(Y, cIdx, scIdx, k=20)

    sdY <- apply(sRes$cY, 2, sd)
    ssd <- sort(sdY,decreasing=TRUE,index.return=TRUE)$ix

    kmresRep <- kmeans(sRes$cY[,ssd[1:nKeep],drop=FALSE],centers=2,nstart=200)
    vclustRep <- kmresRep$cluster
    RepScore <- clScore(vclustRep,X)

    svdResRep <- NULL
    svdResRep <- svdPlot(sRes$cY[, ssd[1:nKeep], drop=FALSE], 
                         annot=annot, 
                         labels=lab.and.region, 
                         svdRes=svdResRep, 
                         plAnnots=plAnnots,                     
                         kColors=gender.col, file=NULL)    

    if (!require(spams))
        {
            print('Iterative-based estimators were not run because they require the spams package to be loaded. Spams can be obtained for free from http://spams-devel.gforge.inria.fr/downloads.html')
            IterScore <- IterRandScore <- NA
        }else {
            ## Iterative replicate-based

            cEps <- 1e-6
            maxIter <- 30
            p <- 20

            paramXb <- list()
            paramXb$K <- p
            paramXb$D <- matrix(c(0.),nrow = 0,ncol=0)
            paramXb$batch <- TRUE
            paramXb$iter <- 1

            ## l1
            paramXb$mode <- 'PENALTY'
            paramXb$lambda <- 0.25 

            iRes <- iterativeRUV(Y, cIdx, scIdx, paramXb, k=20, nu.coeff=0,
                                 cEps, maxIter,
                                 Wmethod='rep', wUpdate=11)

            ucY <- iRes$cY

            sdY <- apply(ucY, 2, sd)
            ssd <- sort(sdY,decreasing=TRUE,index.return=TRUE)$ix

            kmresIter <- kmeans(ucY[,ssd[1:nKeep]],centers=2,nstart=200)
            vclustIter <- kmresIter$cluster
            IterScore <- clScore(vclustIter,X)

            svdResIter <- NULL
            svdResIter <- svdPlot(ucY[, ssd[1:nKeep], drop=FALSE], 
                                  annot=annot, 
                                  labels=lab.and.region, 
                                  svdRes=svdResIter, 
                                  plAnnots=plAnnots,                     
                                  kColors=gender.col, file=NULL)    

            ## Iterated ridge

            paramXb <- list()
            paramXb$K <- p
            paramXb$D <- matrix(c(0.),nrow = 0,ncol=0)
            paramXb$batch <- TRUE
            paramXb$iter <- 1
            paramXb$mode <- 'PENALTY' #2
            paramXb$lambda <- 1
            paramXb$lambda2 <- 0

            iRes <- iterativeRUV(Y, cIdx, scIdx=NULL, paramXb, k=nrow(Y), nu.coeff=1e-2/2,
                                 cEps, maxIter,
                                 Wmethod='svd', wUpdate=11)

            nrcY <- iRes$cY

            sdY <- apply(nrcY, 2, sd)
            ssd <- sort(sdY,decreasing=TRUE,index.return=TRUE)$ix

            kmresIter <- kmeans(nrcY[,ssd[1:nKeep]],centers=2,nstart=200)
            vclustIter <- kmresIter$cluster
            IterRandScore <- clScore(vclustIter,X)

            svdResIterRand <- NULL
            svdResIterRand <- svdPlot(nrcY[, ssd[1:nKeep], drop=FALSE], 
                                      annot=annot, 
                                      labels=lab.and.region, 
                                      svdRes=svdResIterRand, 
                                      plAnnots=plAnnots,                     
                                      kColors=gender.col, file=NULL)    
        }
    
    scores <- c(uScore, MCScore, nsScore, nScore, RepScore, IterScore, IterRandScore)
    names(scores) <- c('Uncorrected', 'Centered', 'Naive RUV-2', 'Naive + shrink', 'Replicates', 'Replicates + iter', 'Shrinkage + iter')
    print('Clustering errors after each correction')
    print(scores)
}
