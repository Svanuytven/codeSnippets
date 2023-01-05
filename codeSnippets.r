# R analysis functions
## Sebastiaan Vanuytven


## Basic function to convert human to mouse gene names
convertHumanGeneList <- function(x){
  
  require("biomaRt")
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  
  genesV2 = getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = x , mart = human, attributesL = c("mgi_symbol"), martL = mouse, uniqueRows=T)
  
  humanx <- unique(genesV2[, 2])
  
  # Print the first 6 genes found to the screen
  print(head(humanx))
  return(humanx)
}

## Functions to determine the optimal number of PCs for downstreamAnalysis
### Initialize Functions For nPC
PrepDR <- function( # From Seurat
  object,
  genes.use = NULL,
  use.imputed = FALSE,
  assay.type="RNA"
) {
  
  if (length(VariableFeatures(object = object)) == 0 && is.null(x = genes.use)) {
    stop("Variable genes haven't been set. Run MeanVarPlot() or provide a vector
         of genes names in genes.use and retry.")
  }
  if (use.imputed) {
    data.use <- t(x = scale(x = t(x = object@imputed)))
  } else {
    data.use <- GetAssayData(object, assay.type = assay.type,slot = "scale.data")
  }
  genes.use <- if(is.null(genes.use)) VariableFeatures(object = object) else genes.use # Changed
  genes.use <- unique(x = genes.use[genes.use %in% rownames(x = data.use)])
  genes.var <- apply(X = data.use[genes.use, ], MARGIN = 1, FUN = var)
  genes.use <- genes.use[genes.var > 0]
  genes.use <- genes.use[! is.na(x = genes.use)]
  data.use <- data.use[genes.use, ]
  return(data.use)
}

### Estimate nPC
PCA_estimate_nPC<-function(data, whereto, k=10, from.nPC = 2, to.nPC=50, by.nPC=1, maxit=200, seed=617) {
  library(missMDA)
  PC <-seq(from = from.nPC, to = to.nPC, by = by.nPC)
  # Init the error matrices
  error1<-matrix(0, nrow = length(c(1:k)), ncol = length(PC))
  error2<-matrix(0, nrow = length(c(1:k)), ncol = length(PC))
  print(paste0(k,"-fold paritioning..."))
  # K-Fold Partitioning
  dgem.kfold<-dismo::kfold(t(data), k=k)
  # SVD-CV based on https://stats.stackexchange.com/questions/93845/how-to-perform-leave-one-out-cross-validation-for-pca-to-determine-the-number-of
  for(i in c(1:k)) {
    print(paste0("k:",i))
    X.train<-t(data[, dgem.kfold!=i])
    X.test<-t(data[, dgem.kfold==i])
    # Find a few approximate singular values and corresponding singular vectors of a matrix.
    print("Running SVD...")
    # Seurat uses IRLBA to do PCA : https://github.com/satijalab/seurat/blob/cec7cb95c73fd6d605723e9af9a1f96eda5635de/R/dimensional_reduction.R
    pca.results<-irlba::irlba(A = X.train, nv = to.nPC, maxit = maxit) # Otherwise, default maxit=100 do not converge
    gl<-pca.results$v
    for(j in 1:length(PC)) {
      print(paste0("Ndims:",PC[j]))
      P<-gl[,1:PC[j]]%*%t(gl[,1:PC[j]])
      # Naive method
      err1<-X.test %*% (diag(dim(P)[1]) - P)
      # Approximate method
      err2<-X.test %*% (diag(dim(P)[1]) - P + diag(diag(P)))
      error1[i,j]<-sum(err1^2)
      error2[i,j]<-sum(err2^2)
      rm(err1)
      rm(err2)
    }
  }
  errors1<-colSums(error1)
  errors2<-colSums(error2)
  nPC=PC[which(errors2 == min(errors2))]
  saveRDS(nPC,whereto)
  plot(PC,errors1)
  plot(PC,errors2)
  return(nPC)
}

## Plotting of maximal 3 different genes using the RGB colouring
RGBColoring <- function(object, coordinates, genes, thr=0.3, slot='scale.data'){
  if (length(genes) == 1){
    red <- genes[1]
    if(slot =='scale.data'){
      geneMat <- object@assays$RNA@scale.data[red,,drop=FALSE]
    } else if (slot == 'data'){
      geneMat <- object@assays$RNA@data[red,,drop=FALSE]
    }
    
    modCols <- list(red=red)
  } else if (length(genes) == 2) {
    red <- genes[1]
    green <- genes[2]
    if(slot =='scale.data'){
      geneMat <- object@assays$RNA@scale.data[c(red,green),,drop=FALSE]
    } else if (slot == 'data'){
      geneMat <- object@assays$RNA@data[c(red,green),,drop=FALSE]
    }
    modCols <- list(red=red, green=green)
  } else if (length(genes) == 3) {
    red <- genes[1]
    green <- genes[2]
    blue <- genes[3]
    if(slot =='scale.data'){
      geneMat <- object@assays$RNA@scale.data[c(red,green,blue),,drop=FALSE]
    } else if (slot == 'data'){
      geneMat <- object@assays$RNA@data[c(red,green,blue),,drop=FALSE]
    }
    modCols <- list(red=red, green=green, blue=blue)
  } else {
    stop('A minimum of 1 and maximum of 3 the genes can be provided.')
  }
  coordinates <- object@reductions[[coordinates]]@cell.embeddings
  geneMat <- t(apply(geneMat, 1, function(x) (x-min(x))/(max(x)-min(x))))
  offColor <- "#c0c0c030" # Transparent light gray
  modCols <- lapply(modCols, function(x) sapply(x, function(gene) rownames(geneMat)[gene]))
  cellColChan <- sapply(modCols, function(modsCol) apply(as.matrix(geneMat[names(modsCol),]), 1, mean))
  if (length(genes) == 1){
    cellCol <- apply(cellColChan, 1, function(x) rgb(x["red"], 0, 0, alpha=.8))
  } else if (length(genes) == 2) {
    cellCol <- apply(cellColChan, 1, function(x) rgb(x["red"], x["green"], 0, alpha=.8))
  } else if (length(genes) == 3) {
    cellCol <- apply(cellColChan, 1, function(x) rgb(x["red"], x["green"], x["blue"], alpha=.8))
  }
  names(cellCol) <- colnames(geneMat)
  coord_sorted <- names(sort(cellCol[rownames(coordinates)]))
  cellCol[as.vector(which(rowSums(cellColChan) < thr))] <- offColor
  plot(coordinates[coord_sorted,], col=cellCol[coord_sorted], pch=16, axes=FALSE, frame.plot=FALSE, ann=FALSE, cex=0.5)
  if (length(genes) == 1){
    legend("topright", legend =  genes, fill=c('red'), cex=0.5)
  } else if (length(genes) == 2) {
    legend("topright", legend =  genes, fill=c('red', 'green'), cex=0.5)
  } else if (length(genes) == 3) {
    legend("topright", legend =  genes, fill=c('red', 'green', 'blue'), cex=0.5)
  }
}


## Function to generate a fancy UMAP plot
fancyUMAP <- function(dataset) {
  library(Seurat) 
  library(ggplot2)
  library(patchwork)
  
  p1 <- DimPlot(dataset, reduction = "UMAP") + theme_void()
  p2 <- ggplot(data.frame(x = 100, y = 100), aes (x = x, y = y)) +
      geom_point. +
      xlim(c(0, 10)) + ylim(c(0, 10)) + theme_classic() +
      ylab ("UMAP 2") + xlab ("UMAP 1") +
      theme ( axis.text. = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks = element_blank(),
          axis.line = element_line(),
          arrow = arrow(angle = 15, length = unit(.5, "cm"), type = "closed"))
  
  layout <- c(
              area(t = 1, 1 = 2, b = 11, r = 11),
              area(t = 10, 1 = 1, b = 12, r = 2))

  pfinal <- p1 + p2 + plot_layout(design = layout)
  
  return(pfinal)
}

## Function to generate a scatterplot with boxplots to the side
plot_scatterplot_with_boxplots <- function(data, x, y, col) {
  library(tidyverse)
  library(patchwork)
  library(ggpubr)

  scatterplot <- ggplot(data, aes(x, y, color = col)) +
      geom_point() +
      geom_density_2d() +
      theme_classic() +
      theme(axis.text = element_blank(), axis.title = element_blank())
  legend <- as_ggplot(get_legend(scatterplot))
  scatterplot <- scatterplot + theme(legend.position = "none")
  boxplot_left <- ggplot(data, aes(col, y, fill = col)) +
      geom_boxplot() +
      theme_classic() +
      theme(legend.position = "none", axis.ticks.x = element_blank(),
              axis.text.x = element_blank(), axis.title.x = element_blank())
  boxplot_bottom <- ggplot(data, aes(col, x, fill = col)) +
      geom_boxplot() +
      theme_classic() +
      theme(legend.position = "none", axis.ticks.y = element_blank(),
          axis.text.y = element_blank(), axis.title.y = element_blank()) +
      coord_flip()

  top <- wrap_plots(boxplot_left, scatterplot, ncol = 2, widths = c(0.2, 0.8))
  bottom <- wrap_plots(legend, boxplot_bottom, ncol = 2, widths = c(0.22, 0.8))
  final <- wrap_plots(top, bottom, nrow = 2, heights = c(0.8, 0.2))
  
  return(final)
}

## Function to generate a ggplot figuyre with a fixed panel size 
plotFixedGGplot <- function(ggplotObject, colInCm = 2, rowInCm = 4.5) {
  library(ggh4x)
  plot <- ggplotObject +  force_panelsizes(cols = unit(colInCm,"cm"),
                                           rows = unit(rowInCm,"cm"))
  
  return(plot)
}
