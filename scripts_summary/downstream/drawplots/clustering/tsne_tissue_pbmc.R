sce.path <- "/stornext/HPCScratch/home/you.y/preprocess_update/SCEs/cluster/"
savepath <- "/stornext/HPCScratch/home/you.y/preprocess_update/results/tsneplots/"
library(scater)
library(scran)
source("/stornext/HPCScratch/home/you.y/preprocess_update/Rscripts/color.R")

setwd(sce.path)

library(ggpubr)
library(CellBench)
library(tidyverse)
#tissue
dir()
readRDS(file.path(sce.path,"mus2_cluster.rds")) ->droplet_tissue
droplet_tissue %>% filter(!is.task_error(result)) ->droplet_tissue

for (i in 1:nrow(droplet_tissue)){
  droplet_tissue[i,]->info
  droplet_tissue$result[[i]]->sce
  sce = runPCA(sce,exprs_values="logcounts",ncomponents = 20)
  runTSNE(sce)->sce
  tr = reducedDim(sce,"TSNE")
  as.character(sce$SNG.1ST) -> sce$SNG.1ST
  sce$SNG.1ST[is.na(sce$SNG.1ST)] <- " "
  tsne1 <- ggplot(data=as.data.frame(tr),aes(x=tr[,1],y=tr[,2],col=sce$SNG.1ST))+
    geom_point(show.legend = F)+
    scale_color_manual(values = tissue_col) +
    labs(x="tSNE-1",y="tSNE-2")+
    theme_bw()+
    theme(text = element_text(size=15),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          axis.text.x=element_blank(),
          legend.position = "none",
          axis.ticks.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank()) +labs(title = paste0(info$data,"-",info$norm_method,"-",info$cluster_method))
  as.character(sce$cluster) ->sce$cluster
  tsne2 <- ggplot(data=as.data.frame(tr),aes(x=tr[,1],y=tr[,2],col=sce$cluster))+
    geom_point(show.legend = F)+
    labs(x="tSNE-1",y="tSNE-2")+
    theme_bw()+
    theme(text = element_text(size=15),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          legend.position = "none",
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank()) +labs(title = paste0(info$data,"-",info$norm_method,"-",info$cluster_method))
  pdf(paste0(savepath,"tissue/",info$design,info$data,"-",info$norm_method,"-",info$cluster_method,".pdf"),width=10,height = 5)
  print(ggarrange(tsne1,tsne2,ncol=2,widths = c(5,5)))
  dev.off()
}




#pbmc5k
readRDS(file.path(sce.path,"pbmc5k_cluster.rds")) ->droplet_pbmc
droplet_pbmc %>% filter(!is.task_error(result)) ->droplet_pbmc

for (i in 1:nrow(droplet_pbmc)){
  droplet_pbmc[i,]->info
  droplet_pbmc$result[[i]]->sce
  sce = runPCA(sce,exprs_values="logcounts",ncomponents = 20)
  runTSNE(sce)->sce
  tr = reducedDim(sce,"TSNE")
  as.character(sce$SNG.1ST) -> sce$SNG.1ST
  sce$SNG.1ST[is.na(sce$SNG.1ST)] <- " "
  tsne1 <- ggplot(data=as.data.frame(tr),aes(x=tr[,1],y=tr[,2],col=sce$SNG.1ST))+
    geom_point(show.legend = F)+
    scale_color_manual(values = pbmc5k_col) +
    labs(x="tSNE-1",y="tSNE-2")+
    theme_bw()+
    theme(text = element_text(size=15),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          axis.text.x=element_blank(),
          legend.position = "none",
          axis.ticks.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank()) +labs(title = paste0(info$data,"-",info$norm_method,"-",info$cluster_method))
  as.character(sce$cluster) ->sce$cluster
  tsne2 <- ggplot(data=as.data.frame(tr),aes(x=tr[,1],y=tr[,2],col=sce$cluster))+
    geom_point(show.legend = F)+
    labs(x="tSNE-1",y="tSNE-2")+
    theme_bw()+
    theme(text = element_text(size=15),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          legend.position = "none",
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank()) +labs(title = paste0(info$data,"-",info$norm_method,"-",info$cluster_method))
  pdf(paste0(savepath,"pbmc/",info$design,info$data,"-",info$norm_method,"-",info$cluster_method,".pdf"),width=10,height = 5)
  print(ggarrange(tsne1,tsne2,ncol=2,widths = c(5,5)))
  dev.off()
}

