#norm by scran
#facet by biotype
#tsne plot
library(scater)
library(edgeR)
library(scran)
library(ggplot2)
library(SingleCellExperiment)
remove_na<- function(sce){
  if("SNG.1ST" %in% colnames(colData(sce))){
    sce <- sce[,!is.na(sce$SNG.1ST)]}
  sce <- sce[,!sce$SNG.1ST==""]
  #sce <- sce[,sce$doublet=="SNG"]
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

SCE_drop.path <- "/stornext/HPCScratch/home/you.y/preprocess_update/SCEs/raw/droplet_based/"
save.path<- "/stornext/HPCScratch/home/you.y/preprocess_update/results/droplet-based/raw_tsne/"


droplet_recode <- function(vec){
  return(recode(vec, "scpipe"="scPipe","zumis"="zUMIs","kb"="kallisto bustools",
                "optimus" ="Optimus","drop"="dropSeqPipe","cellranger"="Cell Ranger",
                "sa"="salmon_SA","splici"="salmon_splici"))
}
dir(SCE_drop.path)
#"mus1","mus2","sc3cl","sc5cl",
designs=c("pbmc5k")
designs=c("pbmc10k")
for (d in designs){
  design=d
  if (d %in% c("mus1","mus2")){
    dataname=c("scpipe","zumis","drop","kb","sa","splici","optimus","cellranger")
    col_use=tissue_col
  }else{
    dataname=c("scpipe","zumis","drop","kb","sa","splici","cellranger")
  }
  if(d=="sc3cl"){col_use=cellcolor}
  if(d=="pbmc5k"){col_use=pbmc5k_col}
  if(d=="pbmc10k"){col_use=pbmc10k_col}
  paste(apply(expand.grid(design, dataname), 1, paste, collapse="_"), ".rds",sep="") ->files
  
  files %in% dir(SCE_drop.path)
  files <- files[files %in% dir(SCE_drop.path)]
  result0 <- lapply(files, function(file){readRDS(file.path(SCE_drop.path,file))})
  
  
  lapply(result0, remove_na) -> result0
  lapply(result0, scran_norm) -> result0
  rep(dataname,each=length(design))->data
  
  library(tidyverse)
  tibble(design=rep(design,length(dataname)), 
         data=droplet_recode(data),
         result=result0) -> datasets
  datasets$design <- recode(datasets$design,"mus1"="10xv2_lungtissue1","mus2"="10xv2_lungtissue2","sc3cl"="10xv2_3cellline",
                            "sc5cl"="10xv2_5cellline","sc5clv3"="10xv3_5cellline","pbmc5k"="10xv3_pbmc5k",
                            "pbmc10k"="10xv3_pbmc10k")
  
  for (i in 1:nrow(datasets)){
    datasets$result[[i]] -> sce
    #protein
    sce <- sce[!is.na(rowData(sce)$gene_biotype),]
    tmp_sce <- sce[rowData(sce)$gene_biotype=="protein_coding",]
    tmp_sce[,apply(counts(tmp_sce),2, sum)>0] -> tmp_sce
    runPCA(tmp_sce,ncomponents =20) ->tmp_sce
    runTSNE(tmp_sce)->tmp_sce
    as.character(tmp_sce$SNG.1ST)->tmp_sce$SNG.1ST
    #tsne
    tr = reducedDim(tmp_sce,"TSNE")
    pdf(paste0(save.path,"/",d,"-",datasets$data[i],"-pro.pdf"),width=3,height = 3)
    tsne <- ggplot(data=as.data.frame(tr),aes(x=tr[,1],y=tr[,2],col=tmp_sce$SNG.1ST))+
      geom_point(show.legend = F)+
      scale_color_manual(values=col_use) +
      labs(x="tSNE-1",y="tSNE-2",title = paste0(d,"-",datasets$data[i],"-protein_coding"))+
      theme_bw()+
      theme(text = element_text(size=12),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"),
            axis.text.x=element_blank(),
            axis.ticks.x=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks.y=element_blank())
    print(tsne)
    dev.off()
    
    
    #lncRNA
    sce <- sce[!is.na(rowData(sce)$gene_biotype),]
    tmp_sce <- sce[rowData(sce)$gene_biotype=="lncRNA",]
    tmp_sce[,apply(counts(tmp_sce),2, sum)>0] -> tmp_sce
    runPCA(tmp_sce,ncomponents =20) ->tmp_sce
    runTSNE(tmp_sce)->tmp_sce
    as.character(tmp_sce$SNG.1ST)->tmp_sce$SNG.1ST
    #tsne
    tr = reducedDim(tmp_sce,"TSNE")
    pdf(paste0(save.path,"/",d,"-",datasets$data[i],"-lnc.pdf"),width=3,height = 3)
    tsne <- ggplot(data=as.data.frame(tr),aes(x=tr[,1],y=tr[,2],col=tmp_sce$SNG.1ST))+
      geom_point(show.legend = F)+
      scale_color_manual(values=col_use) +
      labs(x="tSNE-1",y="tSNE-2",title = paste0(d,"-",datasets$data[i],"-lncRNA"))+
      theme_bw()+
      theme(text = element_text(size=12),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"),
            axis.text.x=element_blank(),
            axis.ticks.x=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks.y=element_blank())
    print(tsne)
    dev.off()
    
    
    #pseudogene
    
    sce <- sce[!is.na(rowData(sce)$gene_biotype),]
    tmp_sce <- sce[grepl("pseudogene",rowData(sce)$gene_biotype),]
    tmp_sce[,apply(counts(tmp_sce),2, sum)>0] -> tmp_sce
    if (nrow(tmp_sce)>100){
    runPCA(tmp_sce,ncomponents =20) ->tmp_sce
    runTSNE(tmp_sce)->tmp_sce
    as.character(tmp_sce$SNG.1ST)->tmp_sce$SNG.1ST
    #tsne
    tr = reducedDim(tmp_sce,"TSNE")
    pdf(paste0(save.path,"/",d,"-",datasets$data[i],"-pseudo.pdf"),width=3,height = 3)
    tsne <- ggplot(data=as.data.frame(tr),aes(x=tr[,1],y=tr[,2],col=tmp_sce$SNG.1ST))+
      geom_point(show.legend = F)+
      scale_color_manual(values=col_use) +
      labs(x="tSNE-1",y="tSNE-2",title = paste0(d,"-",datasets$data[i],"-pseudogene"))+
      theme_bw()+
      theme(text = element_text(size=12),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"),
            axis.text.x=element_blank(),
            axis.ticks.x=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks.y=element_blank())
    print(tsne)
    dev.off()
    }
  }

}
