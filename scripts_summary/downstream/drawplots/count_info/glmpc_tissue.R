# glmpca analysis.
library(cluster)
library(glmpca)
library(CellBench)
library(SingleCellExperiment)
library(scran)
library(scater)
library(R.utils)
library(tidyverse)
remove_na<- function(sce){
  if("SNG.1ST" %in% colnames(colData(sce))){
    sce <- sce[,!is.na(sce$SNG.1ST)]}
  
  if("group" %in% colnames(colData(sce))){
    sce <- sce[,!is.na(sce$group)]}
  
  return(sce)
}

# protein
pro_sil <- function(sce) {
  sce <- sce[!is.na(rowData(sce)$gene_biotype),]
  sce <- sce[,!sce$SNG.1ST==""]
  tmp_sce <- sce[rowData(sce)$gene_biotype=="protein_coding",]
  tmp_sce[,apply(counts(tmp_sce),2, sum)>5] -> tmp_sce
  tmp_sce[apply(counts(tmp_sce),1, sum)>1,] -> tmp_sce
  glmpca(counts(tmp_sce), 20) -> pcs
  sil = silhouette(as.numeric(factor(tmp_sce$SNG.1ST)),  dist(pcs$factors))
  ret_val=(mean(as.data.frame(sil[1:nrow(sil),])[,3]))
  return(ret_val)
}
# lncRNA
lnc_sil <- function(sce) {
  sce <- sce[!is.na(rowData(sce)$gene_biotype),]
  sce <- sce[,!sce$SNG.1ST==""]
  tmp_sce <- sce[rowData(sce)$gene_biotype=="lncRNA",]
  tmp_sce[,apply(counts(tmp_sce),2, sum)>5] -> tmp_sce
  tmp_sce[apply(counts(tmp_sce),1, sum)>1,] -> tmp_sce
  glmpca(counts(tmp_sce), 20) -> pcs
  sil = silhouette(as.numeric(factor(tmp_sce$SNG.1ST)),  dist(pcs$factors))
  ret_val=(mean(as.data.frame(sil[1:nrow(sil),])[,3]))
  return(ret_val)
}


# pseudogene
pseudo_sil <- function(sce) {
  sce <- sce[!is.na(rowData(sce)$gene_biotype),]
  sce <- sce[,!sce$SNG.1ST==""]
  tmp_sce <- sce[grepl("pseudogene",rowData(sce)$gene_biotype),]
  
  if (nrow(tmp_sce)>100){
  tmp_sce[,apply(counts(tmp_sce),2, sum)>5] -> tmp_sce
  tmp_sce[apply(counts(tmp_sce),1, sum)>1,] -> tmp_sce
  glmpca(counts(tmp_sce), 20) -> pcs
  sil = silhouette(as.numeric(factor(tmp_sce$SNG.1ST)),  dist(pcs$factors))
  ret_val=(mean(as.data.frame(sil[1:nrow(sil),])[,3]))}else{
    ret_val=NA
  }
  return(ret_val)
}


glmpca_sil <- list(
  pro_sil=pro_sil,
  lnc_sil=lnc_sil,
  pseudo_sil=pseudo_sil)


SCE_drop.path <- "/stornext/HPCScratch/home/you.y/preprocess_update/SCEs/raw/droplet_based/"
write.path <- "/stornext/HPCScratch/home/you.y/preprocess_update/SCEs/glmpc_sil/"


droplet_recode <- function(vec){
  return(recode(vec, "scpipe"="scPipe","zumis"="zUMIs","kb"="kallisto bustools",
                "optimus" ="Optimus","drop"="dropSeqPipe","cellranger"="Cell Ranger",
                "sa"="salmon_SA","splici"="salmon_splici"))
}
dir(SCE_drop.path)
#"mus1","mus2","pbmc5k",
designs=c("pbmc10k")

for (d in designs){
  design=d
  if (d %in% c("mus1","mus2")){
    dataname=c("scpipe","zumis","drop","kb","sa","splici","optimus","cellranger")
  }else{
    dataname=c("scpipe","zumis","drop","kb","sa","splici","cellranger")
  }

  paste(apply(expand.grid(design, dataname), 1, paste, collapse="_"), ".rds",sep="") ->files
  
  files %in% dir(SCE_drop.path)
  files <- files[files %in% dir(SCE_drop.path)]
  result0 <- lapply(files, function(file){readRDS(file.path(SCE_drop.path,file))})
  
  
  lapply(result0, remove_na) -> result0
  rep(dataname,each=length(design))->data
  
  library(tidyverse)
  tibble(design=rep(design,length(dataname)), 
         data=droplet_recode(data),
         result=result0) -> datasets
  datasets$design <- recode(datasets$design,"mus1"="10xv2_lungtissue1","mus2"="10xv2_lungtissue2","sc3cl"="10xv2_3cellline",
                            "sc5cl"="10xv2_5cellline","sc5clv3"="10xv3_5cellline","pbmc5k"="10xv3_pbmc5k",
                            "pbmc10k"="10xv3_pbmc10k")
  
  result1 <- datasets %>% arrange(design,data) %>% 
    apply_methods(glmpca_sil)
  saveRDS(result1,paste0(write.path,d,".rds"))
  
}
