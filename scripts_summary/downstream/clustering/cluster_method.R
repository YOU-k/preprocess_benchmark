library(Seurat)
library(scran)
##seurat_pipe 1- Louvain
seurat_pipe_louvain <- function(sce){
  real_num <- length(unique(na.omit(sce$SNG.1ST)))
  as.matrix(counts(sce))->counts(sce)
  as.matrix(logcounts(sce)) -> logcounts(sce)
  as.Seurat(sce)->seu
  seu <- ScaleData(object = seu)
  seu <- RunPCA(object = seu,  npcs = 20, verbose = FALSE, features = rownames(seu))
  seu <- FindNeighbors(seu, reduction = "pca", dims = 1:20)
  seu <- FindClusters(seu, resolution = c(0.15,seq(0.2,2.0,0.2)), algorithm = 1)
  list(seu$originalexp_snn_res.0.15,
       seu$originalexp_snn_res.0.2,
       seu$originalexp_snn_res.0.4,
       seu$originalexp_snn_res.0.6,
       seu$originalexp_snn_res.0.8,
       seu$originalexp_snn_res.1,
       seu$originalexp_snn_res.1.2,seu$originalexp_snn_res.1.4,
       seu$originalexp_snn_res.1.6,seu$originalexp_snn_res.1.8,seu$originalexp_snn_res.2)-> clust
  unlist(lapply(clust, function(x) {length(levels(x))}))-real_num -> diff
  which.min(abs(diff)) -> whichk
  whichk
  as.SingleCellExperiment(seu)->sce
  sce$cluster <- factor(clust[[whichk]])
  return(sce)
}


##seurat_pipe 3- SLM
seurat_pipe_slm <-  function(sce){
  real_num <- length(unique(na.omit(sce$SNG.1ST)))
  as.matrix(counts(sce))->counts(sce)
  as.matrix(logcounts(sce)) -> logcounts(sce)
  as.Seurat(sce)->seu
  seu <- ScaleData(object = seu)
  seu <- RunPCA(object = seu,  npcs = 20, verbose = FALSE, features = rownames(seu))
  seu <- FindNeighbors(seu, reduction = "pca", dims = 1:20)
  seu <- FindClusters(seu, resolution = c(0.15,seq(0.2,2.0,0.2)), algorithm = 3)
  list(seu$originalexp_snn_res.0.15,
       seu$originalexp_snn_res.0.2,
       seu$originalexp_snn_res.0.4,
       seu$originalexp_snn_res.0.6,
       seu$originalexp_snn_res.0.8,
       seu$originalexp_snn_res.1,
       seu$originalexp_snn_res.1.2,seu$originalexp_snn_res.1.4,
       seu$originalexp_snn_res.1.6,seu$originalexp_snn_res.1.8,seu$originalexp_snn_res.2)-> clust
  unlist(lapply(clust, function(x) {length(levels(x))}))-real_num -> diff
  which.min(abs(diff)) -> whichk
  whichk
  as.SingleCellExperiment(seu)->sce
  sce$cluster <- factor(clust[[whichk]])
  return(sce)
}



#scran methods with k=5,10,20,30 , use the one with cluster number closest to real value.
#fast-greedy, louvain, walktrap
scran_walktrap <- function(sce){
  real_num <- length(unique(na.omit(sce$SNG.1ST)))
  scater::runPCA(sce,ncomponents =20)->sce
  set.seed(8)
  k=c(100,50,30,10,5)
  clust <- list()
  for (i in 1:5){
    g <- buildSNNGraph(sce,k=k[i],use.dimred="PCA")
    clust[[i]] <- igraph::cluster_walktrap(g)$membership
  }
  unlist(lapply(clust, max))-real_num -> diff
  which.min(abs(diff)) -> whichk
  whichk
  metadata(sce)$k <- k[whichk]
  sce$cluster <- factor(clust[[whichk]])
  return(sce)
}


scran_fastgreedy <- function(sce){
  real_num <- length(unique(na.omit(sce$SNG.1ST)))
  scater::runPCA(sce,ncomponents =20)->sce
  set.seed(8)
  k=c(100,50,30,10,5)
  clust <- list()
  for (i in 1:5){
    g <- buildSNNGraph(sce,k=k[i],use.dimred="PCA")
    clust[[i]] <- igraph::cluster_fast_greedy(g)$membership
  }
  unlist(lapply(clust, max))-real_num -> diff
  which.min(abs(diff)) -> whichk
  metadata(sce)$k <- k[whichk]
  sce$cluster <- factor(clust[[whichk]])
  return(sce)
}

scran_louvain <- function(sce){
  real_num <- length(unique(na.omit(sce$SNG.1ST)))
  scater::runPCA(sce,ncomponents =20)->sce
  set.seed(8)
  k=c(100,50,30,10,5)
  clust <- list()
  for (i in 1:5){
    g <- buildSNNGraph(sce,k=k[i],use.dimred="PCA")
    clust[[i]] <- igraph::cluster_louvain(g)$membership
  }
  unlist(lapply(clust, max))-real_num -> diff
  which.min(abs(diff)) -> whichk
  metadata(sce)$k <- k[whichk]
  sce$cluster <- factor(clust[[whichk]])
  return(sce)
}

#SC3 estimate cluster numbers but do not use, use given number
library(SC3)

sc3_cluster <- function(sce){
  real_num <- length(unique(na.omit(sce$SNG.1ST)))
  scater::runPCA(sce,ncomponents =20)->sce
  counts(sce) <- as.matrix(counts(sce))
  logcounts(sce) <- as.matrix(logcounts(sce))
  rowData(sce)$feature_symbol <- rownames(sce)
  sce <- sc3(sce, ks = real_num, biology = FALSE, n_cores = 8)
  sce$cluster <- last(colData(sce))
  return(sce)
}

sc3_svm <-function(sce){
  real_num <- length(unique(na.omit(sce$SNG.1ST)))
  scater::runPCA(sce,ncomponents =20)->sce
  counts(sce) <- as.matrix(counts(sce))
  logcounts(sce) <- as.matrix(logcounts(sce))
  rowData(sce)$feature_symbol <- rownames(sce)
  sce <- sc3(sce, ks =(real_num-1):real_num, biology = TRUE, svm_num_cells = 50)
  sce <- sc3_run_svm(sce, ks = (real_num-1):real_num)
  metadata(sce)$sc3$svm_train_inds <- NULL
  sce <- sc3_calc_biology(sce, ks = (real_num-1):real_num)
  paste0("sc3_",real_num,"_clusters")->en
  names(colData(sce))[names(colData(sce))==en] <- "cluster"
  return(sce)
}

#RaceID use given number

library(RaceID)
RaceID <- function(sce){
  real_num <- length(unique(na.omit(sce$SNG.1ST)))
  sc <- SCseq(as.data.frame(as.matrix(counts(sce))))
  sc <- filterdata(sc, mintotal=1, minexpr = 1, minnumber = 1,
                   LBatch = NULL, knn = 10, CGenes = NULL, FGenes = NULL, ccor = 0.4,
                   bmode = "RaceID")
  sc@ndata = sc@expdata
  sc@genes = rownames(sc@ndata)
  sc@counts = rep(1,ncol(sce))
  names(sc@counts) = colnames(sc@ndata)
  sc@cluster$features = sc@genes
  sc <- compdist(sc, metric="pearson", FSelect = FALSE, knn = NULL)
  sc <- clustexp(sc, cln = real_num, clustnr = 30,
                 bootnr = 50, rseed = 17000, FUNcluster = "kmedoids")
  sc <- findoutliers(sc, probthr = 0.001, outminc = 5, outlg = 2,
                     outdistquant = 0.95)
  colData(sce)$cluster = as.factor(sc@cpart)
  return(sce)
}



#evaluate

library(mclust)
ARI_matric = function(sce){
  if(!("cluster" %in% colnames(colData(sce)))){
    return(NA)
  }
  if ("group" %in% colnames(colData(sce))){
    sce<- sce[,!is.na(sce$group)]
    ari_val = adjustedRandIndex(sce$group, sce$cluster)
  }else{
    sce<- sce[,!is.na(sce$SNG.1ST)]
    ari_val = adjustedRandIndex(sce$SNG.1ST, sce$cluster)
  }
  
  return(ari_val)
}

cluster_number = function(sce){
  if(!("cluster" %in% colnames(colData(sce)))){
    return(NA)
  }
  return(length(table(sce$cluster)))
}

cluster_filtered_number = function(sce){
  if(!("cluster" %in% colnames(colData(sce)))){
    return(NA)
  }
  if ("group" %in% colnames(colData(sce))){
    sce<- sce[,!is.na(sce$group)]
  }else{
    sce<- sce[,!is.na(sce$SNG.1ST)]
  }
  return(length(table(sce$cluster)))
}

cal_entropy=function(x){
  freqs <- table(x)/length(x)
  freqs = freqs[freqs>0]
  return(-sum(freqs * log(freqs)))
}


get_cluster_purity=function(sce){
  if(!("cluster" %in% colnames(colData(sce)))){
    return(NA)
  }
  if ("group" %in% colnames(colData(sce))){
    sce<- sce[,!is.na(sce$group)]
    return(mean(unlist(lapply(unique(sce$group),function(x){cal_entropy(sce$cluster[sce$group==x])}))))
  }else{
    sce<- sce[,!is.na(sce$SNG.1ST)]
    return(mean(unlist(lapply(unique(sce$SNG.1ST),function(x){cal_entropy(sce$cluster[sce$SNG.1ST==x])}))))
  }
  
}

get_cluster_accuracy=function(sce){
  if(!("cluster" %in% colnames(colData(sce)))){
    return(NA)
  }
  if ("group" %in% colnames(colData(sce))){
    sce<- sce[,!is.na(sce$group)]
    return(mean(unlist(lapply(unique(sce$cluster),function(x){cal_entropy(sce$group[sce$cluster==x])}))))
  }else{
    sce<- sce[,!is.na(sce$SNG.1ST)]
    return(mean(unlist(lapply(unique(sce$cluster),function(x){cal_entropy(sce$SNG.1ST[sce$cluster==x])}))))
  }
  
}

