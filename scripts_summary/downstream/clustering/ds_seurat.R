#seurat with resolution at 0.1,0.2,0.4,0.6,0.8,1.0
#louvain SLM

write.path <- "/stornext/General/data/user_managed/grpu_mritchie_1/Yue/preprocess/new/SCE/cluster"
read.path <-"/stornext/General/data/user_managed/grpu_mritchie_1/Yue/preprocess/new/SCE/raw"

library(Seurat)
library(tidyverse)
library(CellBench)

design=c("sc3cl","sc5clv2","sc5clv3")
dataname=c("kb","scpipe","zumis","optimus","drop","cellranger","alevin")
paste(apply(expand.grid(design, dataname), 1, paste, collapse="_"), ".rds",sep="") ->files


result0 <- lapply(files, function(file){readRDS(file.path(read.path,file))})

remove_doublet <- function(sce){
  sce$demuxlet_cls[is.na(sce$demuxlet_cls)] <- "no"
  sce <- sce[,!sce$demuxlet_cls=="DBL"]
  return(sce)
}

lapply(result0, remove_doublet) -> result0

rep(dataname,each=length(design))->data

tibble(design=rep(design,length(dataname)), 
       data=data,
       result=result0) -> datasets

#recode data and design
datasets$data <- recode(datasets$data,"scpipe"="scPipe","zumis"="zUMIs","kb"="kallisto bustools",
                        "optimus" ="Optimus","drop"="dropSeqPipe","cellranger"="Cell Ranger",
                        "alevin"="salmon alevin")
datasets$design <- recode(datasets$design,"sc3cl"="10Xv2_3cl","sc5clv2"="10Xv2_5cl","sc5clv3"="10Xv3_5cl",
                          "mus1"="10Xv2_tissue1","mus2"="10Xv2_tissue2")

##seurat_pipe 1- Louvain
seurat_pipe_louvain <- function(sce){
  real_num <- length(unique(na.omit(sce$SNG.1ST)))
  CreateSeuratObject(counts = counts(sce))->seu
  seu@meta.data <-cbind(seu@meta.data,list(colData(sce)))
  seu<- NormalizeData(object = seu, normalization.method = "LogNormalize", scale.factor = 10000)
  seu <- FindVariableFeatures(object = seu, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5, nfeatures = 1500)
  seu <- ScaleData(object = seu)
  seu <- RunPCA(object = seu,  npcs = 20, verbose = FALSE)
  seu <- FindNeighbors(seu, reduction = "pca", dims = 1:20)
  seu <- FindClusters(seu, resolution = c(0.01,0.1,0.2,0.4,0.6,0.8,1.0), algorithm = 1)
  list(seu$RNA_snn_res.0.01,seu$RNA_snn_res.0.1,
       seu$RNA_snn_res.0.2,seu$RNA_snn_res.0.4,seu$RNA_snn_res.0.6,seu$RNA_snn_res.0.8,seu$RNA_snn_res.1)-> clust
  unlist(lapply(clust, function(x) {length(levels(x))}))-real_num -> diff
  which.min(abs(diff)) -> whichk
  as.SingleCellExperiment(seu)->sce
  sce$cluster <- factor(clust[[whichk]])
  return(sce)
}


##seurat_pipe 3- SLM
seurat_pipe_slm <- function(sce){
  real_num <- length(unique(na.omit(sce$SNG.1ST)))
  CreateSeuratObject(counts = counts(sce))->seu
  seu@meta.data <-cbind(seu@meta.data,list(colData(sce)))
  seu<- NormalizeData(object = seu, normalization.method = "LogNormalize", scale.factor = 10000)
  seu <- FindVariableFeatures(object = seu, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5, nfeatures = 1500)
  seu <- ScaleData(object = seu)
  seu <- RunPCA(object = seu,  npcs = 20, verbose = FALSE)
  seu <- FindNeighbors(seu, reduction = "pca", dims = 1:20)
  seu <- FindClusters(seu, resolution = c(0.01,0.1,0.2,0.4,0.6,0.8,1.0), algorithm = 3)
  list(seu$RNA_snn_res.0.01,seu$RNA_snn_res.0.1,
       seu$RNA_snn_res.0.2,seu$RNA_snn_res.0.4,seu$RNA_snn_res.0.6,seu$RNA_snn_res.0.8,seu$RNA_snn_res.1)-> clust
  unlist(lapply(clust, function(x) {length(levels(x))}))-real_num -> diff
  which.min(abs(diff)) -> whichk
  as.SingleCellExperiment(seu)->sce
  sce$cluster <- factor(clust[[whichk]])
  return(sce)
}

cluster_method <- list(
  seurat_pipe_louvain=seurat_pipe_louvain,
  seurat_pipe_slm=seurat_pipe_slm
)

result2 <- datasets %>% mutate(norm_method="seurat") %>% select(design,data,norm_method,result) %>% 
  arrange(design,data) %>% apply_methods(cluster_method)

saveRDS(result2,file.path(write.path,"droplet_singlecell_seurat.rds"))
###evaluation
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

clustering_evaluation <- list(
  ARI=ARI_matric,
  cluster_purity=get_cluster_purity,
  cluster_accuracy=get_cluster_accuracy,
  cluster_number=cluster_number,
  cluster_filtered_number=cluster_filtered_number
)

result3 <- result2 %>% apply_methods(clustering_evaluation)
saveRDS(result3,file.path(write.path,"droplet_singlecell_seu_evals.rds"))


