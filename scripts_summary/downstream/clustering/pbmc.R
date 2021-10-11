library(scran)
library(scater)
library(CellBench)
library(tidyverse)
library(Seurat)
write.path <- "/stornext/HPCScratch/home/you.y/preprocess_update/SCEs/cluster/"
read.path <-"/stornext/HPCScratch/home/you.y/preprocess_update/SCEs/norm/hvg/"
source("/stornext/HPCScratch/home/you.y/preprocess_update/Rscripts/clustering/droplet/cluster_method.R")

cluster_method <- list(
  seurat_louvain=seurat_pipe_louvain,
  seurat_slm=seurat_pipe_slm,
  scran_walktrap = scran_walktrap,
  scran_fastgreedy= scran_fastgreedy,
  scran_louvain=scran_louvain,
  sc3 = sc3_cluster,
  sc3_SVM= sc3_svm,
  RaceID=RaceID
)


clustering_evaluation <- list(
  ARI=ARI_matric,
  cluster_purity=get_cluster_purity,
  cluster_accuracy=get_cluster_accuracy,
  cluster_number=cluster_number,
  cluster_filtered_number=cluster_filtered_number
)

designs=c("pbmc5k")
for (d in designs){
  tmp <- readRDS(paste0(read.path,d,"_hvg.rds"))
  tmp %>% mutate(result=scran_hi) %>% select(-scran_hi) %>% arrange(design,data,norm_method) %>% apply_methods(cluster_method) ->result2
  
  saveRDS(result2,paste0(write.path,d,"_cluster.rds"))
  
  table(result2$norm_method)
  
  ###evaluation
  
  result3 <- result2 %>% apply_methods(clustering_evaluation)
  saveRDS(result3,paste0(write.path,d,"_cluster_eval.rds"))
  
}
