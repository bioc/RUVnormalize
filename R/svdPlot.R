#########################################################################/**
## @RdocFunction svdPlot
##
## @title "Plot the data projected into the space spanned by their first two principal components"
##
## \description{ The function takes as input a gene expression matrix
## and plots the data projected into the space spanned by their first
## two principal components.}
##
## \arguments{
##   \item{Y}{Expression matrix where the rows are the
##   samples and the columns are the genes.}
##   \item{annot}{A matrix
##   containing the annotation to be plotted. Each row must correspond
##   to a sample (row) of argument Y, each column must be a
##   categorical or continuous descriptor for the sample. Optional.}
##   \item{labels}{A vector with one element per sample (row) of
##   argument Y. If this argument is specified, each sample is
##   represented by its label. Otherwise, it is represented by a dot
##   (if no annotation is provided) or by the value of the
##   annotation. Optional.}
##   \item{svdRes}{A list containing the result of svd(Y), possibly
##   restricted to the first few singular values. Optional: if not
##   provided, the function computes the SVD.}
##   \item{plAnnots}{A list specifiying whether each column of the
##   annot argument corresponds to a categorical or continuous
##   factor. Each element of the list is named after a column of
##   annot, and contains a string 'categorical' or 'continuous'. For
##   each element of this list, a plot is produced where the samples
##   are represented by colors corresponding to their annotation. Optional.}
##   \item{kColors}{A vector of colors to be used to represent
##   categorical factors. Optional: a default value is provided. If a
##   categorical factors has more levels than the number of colors
##   provided, colors are not used and the factor is represented in
##   black.}
##   \item{file}{A string giving the path to a pdf file for the plot. Optional.}
## }
##
## \value{
##  A @list containing the result of svd(Y, nu=2, nv=0).
##  }
##
## @synopsis
##
## @examples "../inst/extdata/naiveRandRUV.Rex"
##
##*/########################################################################

svdPlot <- function(Y, annot=NULL, labels=NULL, svdRes=NULL, plAnnots=NULL, kColors=NULL, file=NULL){

    if(is.null(kColors)){
        ## This is dichromat::colorschemes$Categorical.12 plus gray and black
        kColors <- c('#FFBF80', '#FF8000', '#FFFF99', '#FFFF33', '#B2FF8C', '#33FF00', '#A6EDFF',
                     '#1AB2FF', '#CCBFFF', '#664CFF', '#FF99BF', '#E61A33', 'gray80', 'black')
        kColors <- rev(kColors)
    }
      
    col.var <- rep(1,nrow(Y))
    col.vec <- c('black')
  
    CEX <- 2

    if(is.null(svdRes)){
        svdRes <- svd(scale(Y,scale=FALSE), nu=2, nv=0)
    }

    if(is.null(plAnnots)){
        if(!is.null(file)){
            pdf(file=paste(file,'.pdf',sep=''), height=8, width=8)
        }
        par(mar=c(5,5,4,2))
        
        plot(svdRes$u[,1:2],type='n',xlab='PC1',ylab='PC2',cex.axis=CEX,cex.lab=CEX)
        if(is.null(labels)){            
            points(svdRes$u[,1:2],col=col.vec[col.var],pch=19,cex=CEX)
        }else{
            text(svdRes$u[,1:2],col=col.vec[col.var],labels=labels,cex=CEX)
        }
        
        if(!is.null(file)){
            dev.off()
        }
    }else{
        textLab <- names(plAnnots)            
        for(tl in textLab){
            if(!is.null(file)){
                pdf(file=paste(file,'_',tl,'.pdf',sep=''), height=8, width=8)
            }
            par(mar=c(5,5,4,2))
            
            col.var <- as.character(annot[, tl])
            
            if(plAnnots[tl] == 'continuous'){
                cNames <- as.character(sort(unique(annot[, tl])))                
                col.vec <- myPalette(low = "#eeeeff",high = "navy", k=length(cNames))
                names(col.vec) <- cNames
                labs <- sapply(annot[, tl], round, digits=2)
            }
            else{
                labs <- annot[,tl] 
                if(length(unique(col.var))<=length(kColors)){
                    col.vec <- kColors[(1:length(unique(col.var)))]
                    names(col.vec) <- unique(col.var)
                }
                else{
                    col.var <- rep(1,nrow(Y))
                    col.vec <- c('black')
                }
            }

            plot(svdRes$u[,1:2],type='n',xlab='PC1',ylab='PC2',cex.axis=CEX,cex.lab=CEX)
            if(is.null(labels)){
                text(svdRes$u[,1:2],col=col.vec[col.var],labels=labs,cex=CEX)
            }else{
                text(svdRes$u[,1:2],col=col.vec[col.var],labels=labels,cex=CEX)
            }

            if(!is.null(file)){
                dev.off()
            }
        }
    }
    svdRes
}
