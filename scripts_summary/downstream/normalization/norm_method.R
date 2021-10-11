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
  if(min(sizeFactors(sce))<=0){
    cleanSizeFactors(sizeFactors(sce),num.detected = sce$detected) ->sizeFactors(sce)
  }
  sce <- logNormCounts(sce) # goes to `logcounts` by default
  return(sce)
}


DESeq2_norm = function(sce){
  sizeFactors(sce) <- estimateSizeFactorsForMatrix(counts(sce)+1)
  sce <- logNormCounts(sce)
  return(sce)
}


linnorm_norm = function(sce){
  sce <- sce[(apply(as.matrix(counts(sce)),1, function(x){max(x)-min(x)>1})),]
  logcounts(sce) = Linnorm(counts(sce))
  return(sce)
}



scone_norm = function(sce){
  sce <- sce[!(apply(as.matrix(counts(sce)),1, function(x){max(x)-min(x)==0})),]
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


sctransform_norm = function(sce){
  seurat <- CreateSeuratObject(counts = counts(sce), project = "seu-sce")
  norm_seuratz<- SCTransform(object = seurat, verbose = FALSE,variable.features.n = 1500)
  #results in top 3000 high variance features
  sce <- sce[rownames(norm_seuratz@assays$SCT@scale.data),]
  logcounts(sce)<-norm_seuratz@assays$SCT@scale.data
  return(sce)
}

sctransform_poi = function(sce){
  seurat <- CreateSeuratObject(counts = counts(sce), project = "seu-sce")
  norm_seuratz<- SCTransform(object = seurat,method = "glmGamPoi", verbose = FALSE,variable.features.n = 1500)
  #results in top 3000 high variance features
  sce <- sce[rownames(norm_seuratz@assays$SCT@scale.data),]
  logcounts(sce)<-norm_seuratz@assays$SCT@scale.data
  return(sce)
}
