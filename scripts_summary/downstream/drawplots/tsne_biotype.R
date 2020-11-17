# check for biotype of each 
library(scran)
scran_norm = function(sce){
  clusters <- quickCluster(sce)
  sce = computeSumFactors(sce,clusters=clusters)
  if(min(sizeFactors(sce))==0){
    cleanSizeFactors(sizeFactors(sce),num.detected = sce$detected) ->sizeFactors(sce)
  }
  sce <- logNormCounts(sce) # goes to `logcounts` by default
  return(sce)
}


SCE.path <- "/stornext/General/data/user_managed/grpu_mritchie_1/Yue/preprocess/new/SCE/raw"
write.path <- "/stornext/General/data/user_managed/grpu_mritchie_1/Yue/preprocess/new/SCE/evals"

#droplet tissue
design=c("mus2")
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

drop_raw <- tibble(
  design=rep(design,length(dataname)),
  data=data,
  result=result0)

norm <- list(scran_norm=scran_norm)
drop_raw  %>% arrange(design,data) %>% apply_methods(norm) ->  drop_raw_norm

savepath <- "/stornext/General/data/user_managed/grpu_mritchie_1/Yue/preprocess/new/results/raw/tsne"
draw_tsne <- function(sce,design,preprocess,biotype){
  as.character(sce$SNG.1ST)->sce$SNG.1ST
  sce$SNG.1ST[sce$SNG.1ST==""] <- "others"
  #tsne
  tr = reducedDim(sce,"TSNE")
  pdf(paste0(savepath,"/",design,"-",preprocess,"-",biotype,".pdf"),width=9,height = 5)
  tsne <- ggplot(data=as.data.frame(tr),aes(x=tr[,1],y=tr[,2],col=sce$SNG.1ST))+
    geom_point(show.legend = T)+
    scale_color_manual(values=tissue_color) +
    labs(x="tSNE-1",y="tSNE-2")+
    theme_bw()+
    theme(text = element_text(size=20),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank())
  print(tsne)
  dev.off()
}

for (i in 3:nrow(drop_raw_norm)) {
  drop_raw_norm$result[[i]]->sce
  rowData(sce)$gene_biotype[is.na(rowData(sce)$gene_biotype)] <- "no"
  
  runPCA(sce,ncomponents =20) ->sce
  runTSNE(sce)->sce
  draw_tsne(sce,drop_raw_norm$design[[i]],drop_raw_norm$data[[i]],"all")
  
  sce[rowData(sce)$gene_biotype=="protein_coding",]->sce3
  runPCA(sce3,ncomponents =20) ->sce3
  runTSNE(sce3)->sce3
  draw_tsne(sce3,drop_raw_norm$design[[i]],drop_raw_norm$data[[i]],"protein_coding")
  
  sce[rowData(sce)$gene_biotype=="lncRNA",]->sce1
  runPCA(sce1,ncomponents =20) ->sce1
  runTSNE(sce1)->sce1
  draw_tsne(sce1,drop_raw_norm$design[[i]],drop_raw_norm$data[[i]],"lncRNA")
  
  sce[grepl("pseudo",rowData(sce)$gene_biotype),]->sce2
  runPCA(sce2,ncomponents =20) ->sce2
  runTSNE(sce2)->sce2
  draw_tsne(sce2,drop_raw_norm$design[[i]],drop_raw_norm$data[[i]],"pseudo")
}





    design=c("sc5clv2")
dataname=c("kb","scpipe","zumis","alevin","optimus","drop","cellranger")

paste(apply(expand.grid(design, dataname), 1, paste, collapse="_"), ".rds",sep="") ->files
#check
files %in% dir(SCE.path)

result0 <- lapply(files, function(file){readRDS(file.path(SCE.path,file))})
remove_doublet <- function(sce){
  if("SNG.1ST" %in% colnames(colData(sce))){
    sce$demuxlet_cls <- as.character(sce$demuxlet_cls)
    sce$demuxlet_cls[is.na(sce$demuxlet_cls)] <- "no"
    sce <- sce[,!sce$demuxlet_cls=="DBL"] }
  
  if("SNG.1ST" %in% colnames(colData(sce))){
    sce <- sce[,!is.na(sce$SNG.1ST)]}
  
  if("group" %in% colnames(colData(sce))){
    sce <- sce[,!is.na(sce$group)]}
  
  return(sce)
}

lapply(result0, remove_doublet) -> result0

rep(dataname,each=length(design))->data

drop_raw <- tibble(
  design=rep(design,length(dataname)),
  data=data,
  result=result0)

norm <- list(scran_norm=scran_norm)
drop_raw  %>% arrange(design,data) %>% apply_methods(norm) ->  drop_raw_norm

draw_tsne <- function(sce,design,preprocess,biotype){
  as.character(sce$SNG.1ST)->sce$SNG.1ST
  #tsne
  tr = reducedDim(sce,"TSNE")
  pdf(paste0(savepath,"/",design,"-",preprocess,"-",biotype,".pdf"),width=7,height = 5)
  tsne <- ggplot(data=as.data.frame(tr),aes(x=tr[,1],y=tr[,2],col=sce$SNG.1ST))+
    geom_point(show.legend = T)+
    scale_color_manual(values=cellcolor) +
    labs(x="tSNE-1",y="tSNE-2")+
    theme_bw()+
    theme(text = element_text(size=20),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank())
  print(tsne)
  dev.off()
}

