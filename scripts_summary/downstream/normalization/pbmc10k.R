library(scran)
library(scater)
library(CellBench)
library(R.utils)
library(tidyverse)
SCE.path <- "/stornext/HPCScratch/home/you.y/preprocess_update/SCEs/raw/droplet_based/"
write.path <- "/stornext/HPCScratch/home/you.y/preprocess_update/SCEs/norm/"

##norm_methods
source("/stornext/HPCScratch/home/you.y/preprocess_update/Rscripts/normalization/norm_method.R")

norm_method <- list(
  none = raw_count,
  scran = scran_norm,
  DESeq2=DESeq2_norm,
  scone=scone_norm,
  #linnorm=linnorm_norm,
  sctransform=sctransform_norm,
  sctransform_poi=sctransform_poi
)




design=c("pbmc10k")
dataname=c("kb","zumis","scpipe","sa","splici","drop","cellranger")
paste(apply(expand.grid(design, dataname), 1, paste, collapse="_"), ".rds",sep="") ->files


#check
files %in% dir(SCE.path)
readRDS("/stornext/HPCScratch/home/you.y/preprocess_update/pbmc_data/10kpbmc/cd_rm.rds") -> cd_rms
sapply(strsplit(as.character(cd_rms),"-"),function(x) x[1]) -> cd_rms

remove_doublet<- function(sce){
  return(sce[,!colnames(sce) %in% cd_rms])
}
result0 <- lapply(files, function(file){readRDS(file.path(SCE.path,file))})
lapply(result0, remove_doublet) -> result0

rep(dataname,each=length(design))->data

tibble(design=rep(design,length(dataname)), 
       data=data,
       result=result0) -> datasets

result1 <- datasets %>% arrange(design,data) %>% 
  apply_methods(norm_method)


saveRDS(result1,file.path(write.path,"pbmc10k_norm.rds"))


#for evaluation, some need to firstly remove unlabeled cells.
library(cluster)
library(scran)
library(scater)
silhouette_pca = function(sce){
  if("logcounts" %in% assayNames(sce)){
    try_res = try({
      sce = runPCA(sce,exprs_values="logcounts",ncomponents = 20)
    })
    if (class(try_res) == "try-error") {
      return(NA)
    }
  }else{
    return(NA)
  }
  ret_val = NA
  try_res = try({
    if("SNG.1ST" %in% colnames(colData(sce))){
      sce <- sce[,!is.na(sce$SNG.1ST)]
      sil = silhouette(as.numeric(factor(sce$SNG.1ST)),  dist(reducedDim(sce,"PCA")))
      ret_val=(mean(as.data.frame(sil[1:nrow(sil),])[,3]))
    }else if("group" %in% colnames(colData(sce))){
      sce <- sce[,!is.na(sce$group)]
      sil = silhouette(as.numeric(factor(sce$group)),  dist(reducedDim(sce,"PCA")))
      ret_val=(mean(as.data.frame(sil[1:nrow(sil),])[,3]))    
    }else{
      ret_val=NA
    }
  })
  if (class(try_res) == "try-error") {
    return(NA)
  }
  return(ret_val)
}

fac_corr_u = function(sce){
  if("logcounts" %in% assayNames(sce)){
    try_res = try({
      sce = runPCA(sce,exprs_values="logcounts",ncomponents = 20)
    })
    if (class(try_res) == "try-error") {
      return(NA)
    }
  }else{
    return(NA)
  }
  score = 0
  pca_res = reducedDim(sce)
  vec = log10(sce$sum)
  for (i in 1:dim(pca_res)[2]){
    r = summary(lm(pca_res[,i]~vec))$r.squared
    score = c(score,r*attr(pca_res,"percentVar")[i])
  }
  ret_val = sum(score)
  return(ret_val)
}

norm_evaluation <- list(
  silhouette=silhouette_pca,
  unwanted_variation_corr=fac_corr_u
)


result1 %>% arrange(design,data,norm_method) %>% apply_methods(norm_evaluation) -> result2
saveRDS(result2,file.path(write.path,"eval/pbmc10k_eval.rds"))
