library(scran)
library(scater)
library(CellBench)
library(R.utils)
library(tidyverse)
SCE.path <- "/stornext/General/data/user_managed/grpu_mritchie_1/Yue/preprocess/new/SCE/raw"
write.path <- "/stornext/General/data/user_managed/grpu_mritchie_1/Yue/preprocess/new/SCE/norm"

##norm_methods
library(DESeq2)
library(scran)
library(edgeR)
library(Linnorm)
library(scone)
library(sctransform)
library(Seurat)

raw_count = function(sce){
  logcounts(sce) = counts(sce)
  return(sce)
}

scran_norm = function(sce){
  clusters <- quickCluster(sce)
  sce = computeSumFactors(sce,clusters=clusters)
  if(min(sizeFactors(sce))==0){
    cleanSizeFactors(sizeFactors(sce),num.detected = sce$detected) ->sizeFactors(sce)
  }
  sce <- logNormCounts(sce) # goes to `logcounts` by default
  return(sce)
}

DESeq2_norm = function(sce){
  new <- SingleCellExperiment(
    assays = list(counts = matrix(1,ncol=ncol(counts(sce)))))
  int_elementMetadata(new) <- int_elementMetadata(sce)[1,]
  colnames(new) <- colnames(sce)
  rownames(new) <-"add"
  rbind(sce,new) ->sce
  sizeFactors(sce) <- estimateSizeFactorsForMatrix(counts(sce))
  sce <- logNormCounts(sce)
  return(sce)
}


linnorm_norm = function(sce){
  if(table(apply(counts(sce), 1, function(x) sum(x>0)/length(x)) >0.4)[2] <400){
    logcounts(sce) = Linnorm(counts(sce),minNonZeroPortion = 0.2)
  }else if(table(apply(counts(sce), 1, function(x) sum(x>0)/length(x)) >0.75)[2] <400){
    logcounts(sce) = Linnorm(counts(sce),minNonZeroPortion = 0.4)
  }else{
    logcounts(sce) = Linnorm(counts(sce))
  }
  return(sce)
}


sctransform_norm = function(sce){
  seurat <- CreateSeuratObject(counts = counts(sce), project = "seu-sce")
  norm_seuratz<- SCTransform(object = seurat, verbose = FALSE,variable.features.n = 1500)
  #results in top 3000 high variance features
  sce <- sce[rownames(norm_seuratz@assays$SCT@scale.data),]
  logcounts(sce)<-norm_seuratz@assays$SCT@scale.data
  return(sce)
}


scone_norm = function(sce){
  scaling=list(none=identity, # Identity - do nothing
               sum = SUM_FN, # SCONE library wrappers...
               tmm = TMM_FN, 
               uq = UQ_FN,
               fq = FQT_FN,
               deseq = DESEQ_FN)
  results = scone(SconeExperiment(as.matrix(counts(sce))), 
                  scaling=scaling,
                  run=TRUE, k_qc=0, k_ruv=0,
                  return_norm = "in_memory",
                  zero = "postadjust",
                  bpparam = BiocParallel::SerialParam())
  out_norm = get_normalized(results,
                            method = rownames(get_params(results))[1])
  logcounts(sce) = log2(out_norm + 1)
  return(sce)
}

norm_method <- list(
  none = raw_count,
  scran = scran_norm,
  DESeq2=DESeq2_norm,
  Linnorm=linnorm_norm,
  scone=scone_norm,
  sctransform=sctransform_norm
)


design=c("mus1","mus2")
dataname=c("kb","scpipe","zumis","optimus","drop","cellranger","alevin")
paste(apply(expand.grid(design, dataname), 1, paste, collapse="_"), ".rds",sep="") ->files


#check
files %in% dir(SCE.path)

result0 <- lapply(files, function(file){readRDS(file.path(SCE.path,file))})


rep(dataname,each=length(design))->data

tibble(design=rep(design,length(dataname)), 
       data=data,
       result=result0) -> datasets

result1 <- datasets %>% arrange(design,data) %>% 
  apply_methods(norm_method) ->result1
saveRDS(result1,file.path(write.path,"droplet_tissue.rds"))

#for evaluation, some need to firstly remove unlabeled cells.
library(cluster)
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

fac_corr_w = function(sce){
  if("SNG.1ST" %in% colnames(colData(sce))){
    sce <- sce[,!is.na(sce$SNG.1ST)]}
  if("group" %in% colnames(colData(sce))){
    sce <- sce[,!is.na(sce$group)]}
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
  
  if("group" %in% colnames(colData(sce))){
    vec = sce$group
  }else{
    vec = sce$SNG.1ST
  }
  pca_res = reducedDim(sce)
  for (i in 1:dim(pca_res)[2]){
    r = summary(lm(pca_res[,i]~vec))$r.squared
    score = c(score,r*attr(pca_res,"percentVar")[i])
  }
  ret_val = sum(score)
  return(ret_val)
}
norm_evaluation <- list(
  silhouette=silhouette_pca,
  unwanted_variation_corr=fac_corr_u,
  wanted_variation_corr=fac_corr_w
)


result1 %>% apply_methods(norm_evaluation) -> result2
saveRDS(result2,file.path(write.path,"eval/droplet_tissue_eval.rds"))
