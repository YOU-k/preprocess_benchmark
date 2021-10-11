sce.path <- "/stornext/HPCScratch/home/you.y/preprocess_update/SCEs/norm/"
savepath <- "/stornext/HPCScratch/home/you.y/preprocess_update/results/pcaplots/"

readRDS(file.path(sce.path,"plate_rnamix_v2.rds")) ->plate_mixture
library(scran)
library(tidyverse)
library(CellBench)
library(scater)

plate_mixture %>% dplyr::filter(!is.task_error(result))-> plate_rnamix

for (i in 1:nrow(plate_rnamix)){
  plate_rnamix$result[[i]]->sce
  sce = runPCA(sce,exprs_values="logcounts",ncomponents = 2)
  col <- rgb(sce$H1975_prop, sce$HCC827_prop, sce$H2228_prop,alpha=0.9)
  #pca
  dr = reducedDim(sce,"PCA")
  pdf(paste0(savepath,"rnamix/",plate_rnamix$design[[i]],plate_rnamix$data[[i]],"-",plate_rnamix$norm_method[[i]],".pdf"),width=5,height = 5)
  pca <- ggplot(data=as.data.frame(dr),aes(x=PC1,y=PC2,col=sce$group))+
    geom_point(show.legend = F)+
    scale_size_continuous(range = c(1,4))+
    scale_color_manual(guide=FALSE,values = unique(col), limits = unique(sce$group))+
    labs(x=paste0("PC1 (",format(attr(dr,"percentVar")[1],digits=3),"%)"),y=paste0("PC2 (",format(attr(dr,"percentVar")[2],digits=3),"%)"))+
    theme_bw()+theme(text = element_text(size=20),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"),
                     axis.text.x=element_blank(),
                     axis.ticks.x=element_blank(),
                     axis.text.y=element_blank(),
                     axis.ticks.y=element_blank())+
    ggtitle(paste0(plate_rnamix$data[[i]],"-",plate_rnamix$norm_method[[i]]))
  print(pca)
  dev.off()
}


#plate_singlecell
readRDS(file.path(sce.path,"plate_singlecell_v2.rds")) ->plate_singlecell
plate_singlecell  %>% dplyr::filter(!is.task_error(result))-> plate_singlecell
cellcolor <- c(rgb(1,0,0,alpha = 0.9),rgb(0,1,0,alpha = 0.9),rgb(0,0,1,alpha = 0.9),"purple","orange","grey")
names(cellcolor) <- c("H1975","HCC827","H2228","H838","A549"," ")
for (i in 1:nrow(plate_singlecell)){
  plate_singlecell$result[[i]]->sce
  as.character(sce$SNG.1ST) -> sce$SNG.1ST
  sce$SNG.1ST[is.na(sce$SNG.1ST)] <- " "
  sce = runPCA(sce,exprs_values="logcounts",ncomponents = 2)
  #pca
  dr = reducedDim(sce,"PCA")
  pdf(paste0(savepath,"simplemix/",plate_singlecell$design[[i]],plate_singlecell$data[[i]],"-",plate_singlecell$norm_method[[i]],".pdf"),width=5,height = 5)
  pca <- ggplot(data=as.data.frame(dr),aes(x=PC1,y=PC2,col=sce$SNG.1ST))+
    geom_point(show.legend = F)+
    scale_size_continuous(range = c(1,4))+
    scale_color_manual(values = cellcolor)+
    labs(x=paste0("PC1 (",format(attr(dr,"percentVar")[1],digits=3),"%)"),y=paste0("PC2 (",format(attr(dr,"percentVar")[2],digits=3),"%)"))+
    theme_bw()+theme(text = element_text(size=20),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"),
                     axis.text.x=element_blank(),
                     axis.ticks.x=element_blank(),
                     axis.text.y=element_blank(),
                     axis.ticks.y=element_blank())+
    ggtitle(paste0(plate_singlecell$data[[i]],"-",plate_singlecell$norm_method[[i]]))
  print(pca)
  dev.off()
}


