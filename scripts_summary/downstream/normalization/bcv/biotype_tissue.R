library(edgeR)
library(scran)
library(tidyverse)
library(scater)
library(CellBench)
SCE.path <- "/stornext/General/data/user_managed/grpu_mritchie_1/Yue/preprocess/new/SCE/raw"
write.path <- "/stornext/General/data/user_managed/grpu_mritchie_1/Yue/preprocess/new/SCE/evals"

#droplet single cell
design=c("mus1","mus2")
dataname=c("kb","scpipe","zumis","alevin","optimus","drop","cellranger")

paste(apply(expand.grid(design, dataname), 1, paste, collapse="_"), ".rds",sep="") ->files
#check
files %in% dir(SCE.path)

result0 <- lapply(files, function(file){readRDS(file.path(SCE.path,file))})
remove_na<- function(sce){
  if("SNG.1ST" %in% colnames(colData(sce))){
    sce <- sce[,!is.na(sce$SNG.1ST)]}
  
  if("group" %in% colnames(colData(sce))){
    sce <- sce[,!is.na(sce$group)]}
  
  return(sce)
}

lapply(result0, remove_na) -> result0


rep(dataname,each=length(design))->data


scran_norm = function(sce){
  clusters <- quickCluster(sce)
  sce = computeSumFactors(sce,clusters=clusters)
  if(min(sizeFactors(sce))==0){
    cleanSizeFactors(sizeFactors(sce),num.detected = sce$detected) ->sizeFactors(sce)
  }
  sce <- logNormCounts(sce,log=FALSE) # goes to `logcounts` by default
  return(log2(rowMeans(normcounts(sce))))
}

pro_genesidex <- function(sce){
  rowData(sce)$gene_biotype[is.na(rowData(sce)$gene_biotype)] <- "no"
  return(rowData(sce)$gene_biotype =="protein_coding")
}


lnc_genesidex <- function(sce){
  rowData(sce)$gene_biotype[is.na(rowData(sce)$gene_biotype)] <- "no"
  return(rowData(sce)$gene_biotype =="lncRNA")
}



pseudo_genesidex <- function(sce){
  rowData(sce)$gene_biotype[is.na(rowData(sce)$gene_biotype)] <- "no"
  return(grepl("pseudo",rowData(sce)$gene_biotype))
}


calculate_bcv_pro <- function(sce){
  rowData(sce)$gene_biotype[is.na(rowData(sce)$gene_biotype)] <- "no"
  sce <- sce[rowData(sce)$gene_biotype =="protein_coding",]
  if("SNG.1ST" %in% colnames(colData(sce))){
    celltype <- factor(sce$SNG.1ST)
  }
  if("group" %in% colnames(colData(sce))){
    celltype <- factor(sce$group)
  }
  design <- model.matrix(~celltype)
  estimateDisp(counts(sce),design=design,robust = TRUE) -> disp
  return(disp$trended.dispersion)
}

calculate_bcv_lnc <- function(sce){
  rowData(sce)$gene_biotype[is.na(rowData(sce)$gene_biotype)] <- "no"
  sce <- sce[rowData(sce)$gene_biotype =="lncRNA",]
  if("SNG.1ST" %in% colnames(colData(sce))){
    celltype <- factor(sce$SNG.1ST)
  }
  if("group" %in% colnames(colData(sce))){
    celltype <- factor(sce$group)
  }
  design <- model.matrix(~celltype)
  estimateDisp(counts(sce),design=design,robust = TRUE) -> disp
  return(disp$trended.dispersion)
}


calculate_bcv_pseudo <- function(sce){
  rowData(sce)$gene_biotype[is.na(rowData(sce)$gene_biotype)] <- "no"
  sce <- sce[grepl("pseudo",rowData(sce)$gene_biotype),]
  if("SNG.1ST" %in% colnames(colData(sce))){
    celltype <- factor(sce$SNG.1ST)
  }
  if("group" %in% colnames(colData(sce))){
    celltype <- factor(sce$group)
  }
  design <- model.matrix(~celltype)
  estimateDisp(counts(sce),design=design,robust = TRUE) -> disp
  return(disp$trended.dispersion)
}


trend <- list(avelognorm=scran_norm,
              pro_genesidex= pro_genesidex,
              lnc_genesidex= lnc_genesidex,
              pseudo_genesidex= pseudo_genesidex,
              calculate_bcv_pro=calculate_bcv_pro,
              calculate_bcv_lnc= calculate_bcv_lnc,
              calculate_bcv_pseudo= calculate_bcv_pseudo)


drop_raw <- tibble(
  design=rep(design,length(dataname)),
  data=data,
  result=result0)

drop_raw %>% arrange(design,data) %>% apply_methods(trend) ->  drop_raw_trend

drop_raw_trend %>% filter(!is.task_error(result))-> drop_raw_trend

drop_raw_trend %>% spread(trend, result) -> drop_raw_trend_spread

drop_raw_trend_spread$data <- recode(drop_raw_trend_spread$data,"scpipe"="scPipe","zumis"="zUMIs","kb"="kallisto bustools",
                                     "optimus" ="Optimus","drop"="dropSeqPipe","cellranger"="Cell Ranger",
                                     "alevin"="salmon alevin")
saveRDS(drop_raw_trend_spread,file.path(write.path,"biotypes/drop_raw_trend_tissue.rds"))
savepath <- "/stornext/General/data/user_managed/grpu_mritchie_1/Yue/preprocess/new/results/norm/bcv"

drop_raw_trend_spread <- drop_raw_trend_spread[!drop_raw_trend_spread$data=="zUMIs",]
drawplot <- function(n){
  drop_raw_trend_spread$ave1 <- drop_raw_trend_spread$avelognorm
  
  for (i in 1:nrow(drop_raw_trend_spread)) {
    drop_raw_trend_spread$ave1[[i]] <- 
      drop_raw_trend_spread$avelognorm[[i]][drop_raw_trend_spread$pro_genesidex[[i]]]
  }
  
  
  drop_raw_trend_spread$ave2 <- drop_raw_trend_spread$avelognorm
  
  for (i in 1:nrow(drop_raw_trend_spread)) {
    drop_raw_trend_spread$ave2[[i]] <- 
      drop_raw_trend_spread$avelognorm[[i]][drop_raw_trend_spread$lnc_genesidex[[i]]]
  }
  
  
  drop_raw_trend_spread$ave3 <- drop_raw_trend_spread$avelognorm
  
  for (i in 1:nrow(drop_raw_trend_spread)) {
    drop_raw_trend_spread$ave3[[i]] <- 
      drop_raw_trend_spread$avelognorm[[i]][drop_raw_trend_spread$pseudo_genesidex[[i]]]
  }
  
  
  drop_raw_trend_spread %>% dplyr::filter(design==n) -> trend
  
  trend1 <- trend[!is.task_error(trend$calculate_bcv_pro),]
  trend1 <- trend1[!trend1$calculate_bcv_pro=="NULL",]
  
  pdf(paste0(savepath,"/biotypes/",n,"-protein.pdf"),width=4,height=3)
  p1=ggplot() +scale_color_manual(values=mycol)+theme_bw() +labs(y="Biological coefficient of variation",title = "Protein-coding genes")
  for (i in 1:nrow(trend1)){
    data.frame(avelognorm=unlist(trend1$ave1[[i]]),bcv=sqrt(unlist(trend1$calculate_bcv_pro[[i]])),preprocess=trend1$data[[i]]) -> d
    p1+ geom_smooth(data=d, aes(x=avelognorm,y=bcv,col=preprocess),span=0.01,method = "loess")->p1
  }
  print(p1)
  dev.off()
  
  trend1 <- trend[!is.task_error(trend$calculate_bcv_lnc),]
  trend1 <- trend1[!trend1$calculate_bcv_lnc=="NULL",]
  
  pdf(paste0(savepath,"/biotypes/",n,"-lnc.pdf"),width=4,height=3)
  p1=ggplot() +scale_color_manual(values=mycol)+theme_bw() +labs(y="Biological coefficient of variation",title = "LncRNAs")
  for (i in 1:nrow(trend1)){
    data.frame(avelognorm=unlist(trend1$ave2[[i]]),bcv=sqrt(unlist(trend1$calculate_bcv_lnc[[i]])),preprocess=trend1$data[[i]]) -> d
    p1+ geom_smooth(data=d, aes(x=avelognorm,y=bcv,col=preprocess),span=0.01,method = "loess")->p1
  }
  print(p1)
  dev.off()
  
  pdf(paste0(savepath,"/biotypes/",n,"-lnc-point.pdf"),width=4,height=3)
  p1=ggplot() +scale_color_manual(values=mycol)+theme_bw() +labs(y="Biological coefficient of variation",title = "LncRNAs")
  for (i in 1:nrow(trend1)){
    data.frame(avelognorm=unlist(trend1$ave2[[i]]),bcv=sqrt(unlist(trend1$calculate_bcv_lnc[[i]])),preprocess=trend1$data[[i]]) -> d
    p1+ geom_point(data=d, aes(x=avelognorm,y=bcv,col=preprocess),size=0.1)->p1
  }
  print(p1)
  dev.off()
  
  
  trend1 <- trend[!is.task_error(trend$calculate_bcv_pseudo),]
  trend1 <- trend1[!trend1$calculate_bcv_pseudo=="NULL",]
  
  pdf(paste0(savepath,"/biotypes/",n,"-pseudo.pdf"),width=4,height=3)
  p1=ggplot() +scale_color_manual(values=mycol)+theme_bw() +labs(y="Biological coefficient of variation",title = "Pseudogenes")
  for (i in 1:nrow(trend1)){
    data.frame(avelognorm=unlist(trend1$ave3[[i]]),bcv=sqrt(unlist(trend1$calculate_bcv_pseudo[[i]])),preprocess=trend1$data[[i]]) -> d
    p1+ geom_smooth(data=d, aes(x=avelognorm,y=bcv,col=preprocess),span=0.01,method = "loess")->p1
  }
  print(p1)
  dev.off()
  
  
  pdf(paste0(savepath,"/biotypes/",n,"-pseudo-point.pdf"),width=4,height=3)
  p1=ggplot() +scale_color_manual(values=mycol)+theme_bw() +labs(y="Biological coefficient of variation",title = "Pseudogenes")
  for (i in 1:nrow(trend1)){
    data.frame(avelognorm=unlist(trend1$ave3[[i]]),bcv=sqrt(unlist(trend1$calculate_bcv_pseudo[[i]])),preprocess=trend1$data[[i]]) -> d
    p1+ geom_point(data=d, aes(x=avelognorm,y=bcv,col=preprocess),size=0.1)->p1
  }
  print(p1)
  dev.off()
  
}
drawplot("mus1")
drawplot("mus2")


