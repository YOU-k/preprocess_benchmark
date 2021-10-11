# droplet - 
# only common features mus1, sc5cl, pbmc5k
#1. normalization 
# 2. facet by biotype. (no filtering, all features are applied)
# 3. calculate bcv
library(scater)
library(edgeR)
library(scran)
library(ggplot2)
library(SingleCellExperiment)
remove_na<- function(sce){
  if("SNG.1ST" %in% colnames(colData(sce))){
    sce <- sce[,!is.na(sce$SNG.1ST)]}
  
  if("group" %in% colnames(colData(sce))){
    sce <- sce[,!is.na(sce$group)]}
  
  return(sce)
}

norm <- function(sce){
  return(logNormCounts(sce,log=FALSE)) # goes to `logcounts` by default)
}

SCE_drop.path <- "/stornext/HPCScratch/home/you.y/preprocess_update/SCEs/raw/droplet_based/"
write.path <- "/stornext/HPCScratch/home/you.y/preprocess_update/results/droplet-based/bcv/"


droplet_recode <- function(vec){
  return(recode(vec, "scpipe"="scPipe","zumis"="zUMIs","kb"="kallisto bustools",
                "optimus" ="Optimus","drop"="dropSeqPipe","cellranger"="Cell Ranger",
                "sa"="salmon_SA","splici"="salmon_splici"))
}
dir(SCE_drop.path)
#
designs=c("mus1","sc5cl","pbmc5k")

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
  tibble(design=rep(design,length(dataname)), 
         data=droplet_recode(data),
         result=result0) -> datasets
  
  Reduce(intersect, lapply(datasets$result[!datasets$data=="Cell Ranger"], rownames)) -> allgenes
  lapply(result0, function(sce){sce[rownames(sce) %in% allgenes,]}) -> result0
  lapply(result0, norm) -> result0
  datasets$result <- result0
  
  library(tidyverse)

  datasets$design <- recode(datasets$design,"mus1"="10xv2_lungtissue1","mus2"="10xv2_lungtissue2","sc3cl"="10xv2_3cellline",
                            "sc5cl"="10xv2_5cellline","sc5clv3"="10xv3_5cellline","pbmc5k"="10xv3_pbmc5k",
                            "pbmc10k"="10xv3_pbmc10k")
  #1. grepl protein coding genes
  bcv_pro <- tibble()
  for (i in 1:nrow(datasets)){
    sce <- result0[[i]]
    sce <- sce[!is.na(rowData(sce)$gene_biotype),]
    tmp_sce <- sce[rowData(sce)$gene_biotype=="protein_coding",]
    tmp_sce[,apply(counts(tmp_sce),2, sum)>0] -> tmp_sce
    if("SNG.1ST" %in% colnames(colData(tmp_sce))){
      celltype <- factor(tmp_sce$SNG.1ST)
    }
    if("group" %in% colnames(colData(tmp_sce))){
      celltype <- factor(tmp_sce$group)
    }
    dd <- model.matrix(~celltype)
    estimateDisp(counts(tmp_sce),design=dd,robust = TRUE) -> disp
    saveRDS(disp,paste0("/stornext/HPCScratch/home/you.y/preprocess_update/SCEs/disp/common/",design,"_",datasets$data[i],"_disp_pro.rds"))
    
    tmp_bcv <- data.frame(method=datasets$data[i],
                          x=log2(rowMeans(normcounts(tmp_sce))),
                          y=disp$trended.dispersion)
    
    bind_rows(tmp_bcv,bcv_pro) -> bcv_pro
  }
  saveRDS(bcv_pro,paste0("/stornext/HPCScratch/home/you.y/preprocess_update/SCEs/bcv/common/",design,"_bcv_pro.rds"))
  
  
  
  #2. grepl lncRNA
  bcv_lnc <- tibble()
  for (i in 1:nrow(datasets)){
    sce <- result0[[i]]
    sce <- sce[!is.na(rowData(sce)$gene_biotype),]
    tmp_sce <- sce[rowData(sce)$gene_biotype=="lncRNA",]
    tmp_sce[,apply(counts(tmp_sce),2, sum)>0] -> tmp_sce
    if("SNG.1ST" %in% colnames(colData(tmp_sce))){
      celltype <- factor(tmp_sce$SNG.1ST)
    }
    if("group" %in% colnames(colData(tmp_sce))){
      celltype <- factor(tmp_sce$group)
    }
    dd <- model.matrix(~celltype)
    estimateDisp(counts(tmp_sce),design=dd,robust = TRUE) -> disp
    saveRDS(disp,paste0("/stornext/HPCScratch/home/you.y/preprocess_update/SCEs/disp/common/",design,"_",datasets$data[i],"_disp_lnc.rds"))
    tmp_bcv <- data.frame(method=datasets$data[i],
                          x=log2(rowMeans(normcounts(tmp_sce))),
                          y=disp$trended.dispersion)
    
    bind_rows(tmp_bcv,bcv_lnc) -> bcv_lnc
  }
  saveRDS(bcv_lnc,paste0("/stornext/HPCScratch/home/you.y/preprocess_update/SCEs/bcv/common/",design,"_bcv_lnc.rds"))
  
  
  #3. grepl pseudogene
  bcv_pseudo <- tibble()
  for (i in 1:nrow(datasets)){
    sce <- result0[[i]]
    sce <- sce[!is.na(rowData(sce)$gene_biotype),]
    tmp_sce <- sce[grepl("pseudogene",rowData(sce)$gene_biotype),]
    
    if (nrow(tmp_sce)>100){
      tmp_sce[,apply(counts(tmp_sce),2, sum)>0] -> tmp_sce
      if("SNG.1ST" %in% colnames(colData(tmp_sce))){
        celltype <- factor(tmp_sce$SNG.1ST)
      }
      if("group" %in% colnames(colData(tmp_sce))){
        celltype <- factor(tmp_sce$group)
      }
      dd <- model.matrix(~celltype)
      estimateDisp(counts(tmp_sce),design=dd,robust = TRUE) -> disp
      saveRDS(disp,paste0("/stornext/HPCScratch/home/you.y/preprocess_update/SCEs/disp/common/",design,"_",datasets$data[i],"_disp_pseudo.rds"))
      tmp_bcv <- data.frame(method=datasets$data[i],
                            x=log2(rowMeans(normcounts(tmp_sce))),
                            y=disp$trended.dispersion)
    }else{tmp_bcv<- tibble()}
    
    bind_rows(tmp_bcv,bcv_pseudo) -> bcv_pseudo
    
    
  }
  saveRDS(bcv_pseudo,paste0("/stornext/HPCScratch/home/you.y/preprocess_update/SCEs/bcv/common/",design,"_bcv_pseudo.rds"))
  
  pdf(paste0(write.path,d,"_bcv_common.pdf"),width = 9,height = 3.5)
  p <- bind_rows(bcv_pro %>% mutate(type="Protein_coding"),
                 bcv_lnc %>% mutate(type="lncRNA"),
                 bcv_pseudo %>% mutate(type="Pseudogenes")) %>% 
    ggplot(aes(x=x,y=y,col=method)) +
    facet_grid(.~type) +
    geom_smooth(span=0.1,method = "loess") +
    scale_color_manual(values=droplet_col)+
    theme(text = element_text(size=14),
          axis.ticks.x=element_blank(),panel.background = element_rect(fill = "white", colour = NA), 
          panel.border = element_rect(fill = NA,colour = "grey20"), panel.grid = element_line(colour = "grey92"), 
          panel.grid.minor = element_line(size = rel(0.5)), 
          strip.background = element_rect(fill = "grey85",colour = "grey20"))+
    labs(y="Biological coefficient of variation",
         title = datasets$design[1],x="Log average normalized counts",col="Preprocess workflow")
  print(p)
  dev.off()
  
  
  
  bcv_pseudo %>% mutate(type="Pseudogenes") %>% 
    ggplot(aes(x=x,y=y,col=method)) +
    facet_grid(.~type) +
    geom_point() +
    scale_color_manual(values=droplet_col)+
    theme(text = element_text(size=14),
          axis.ticks.x=element_blank(),panel.background = element_rect(fill = "white", colour = NA), 
          panel.border = element_rect(fill = NA,colour = "grey20"), panel.grid = element_line(colour = "grey92"), 
          panel.grid.minor = element_line(size = rel(0.5)), 
          strip.background = element_rect(fill = "grey85",colour = "grey20"))+
    labs(y="Biological coefficient of variation",
         title = datasets$design[1],x="Log average normalized counts",col="Preprocess workflow")
  
}


