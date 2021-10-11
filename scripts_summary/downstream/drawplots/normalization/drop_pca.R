# for droplet based pca plots
sce.path <- "/stornext/HPCScratch/home/you.y/preprocess_update/SCEs/norm/"
savepath <- "/stornext/HPCScratch/home/you.y/preprocess_update/results/pcaplots/"

library(scran)
library(tidyverse)
library(CellBench)
library(scater)

dir(sce.path)
#sc3cl
designs=c("sc3cl","sc5cl","sc5clv3")
for (d in designs) {
  readRDS(paste0(sce.path,d,"_norm.rds")) -> droplet_singlecell
  droplet_singlecell  %>% dplyr::filter(!is.task_error(result))-> droplet_singlecell
  
  cellcolor <- c(rgb(1,0,0,alpha = 0.9),rgb(0,1,0,alpha = 0.9),rgb(0,0,1,alpha = 0.9),"purple","orange","grey")
  names(cellcolor) <- c("H1975","HCC827","H2228","H8383","A549"," ")
  
  for (i in 1:nrow(droplet_singlecell)){
    droplet_singlecell$result[[i]]->sce
    as.character(sce$SNG.1ST)->sce$SNG.1ST
    sce$SNG.1ST[(is.na(sce$SNG.1ST))] <- " "
    sce = runPCA(sce,exprs_values="logcounts",ncomponents = 2)
    #pca
    dr = reducedDim(sce,"PCA")
    pdf(paste0(savepath,"10x_cellline/",droplet_singlecell$design[[i]],droplet_singlecell$data[[i]],"-",droplet_singlecell$norm_method[[i]],".pdf"),width=5,height = 5)
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
      ggtitle(paste0(droplet_singlecell$data[[i]],"-",droplet_singlecell$norm_method[[i]]))
    print(pca)
    dev.off()
  }
}


#mus1
designs=c("mus1","mus2")
for (d in designs) {
  readRDS(paste0(sce.path,d,"_norm.rds")) -> mus
  mus  %>% dplyr::filter(!is.task_error(result))-> mus

for (i in 1:nrow(mus)){
  mus$result[[i]]->sce
  #sce[,!is.na(sce$SNG.1ST)] ->sce
  #sce[,!sce$SNG.1ST==""]->sce
  as.character(sce$SNG.1ST)->sce$SNG.1ST

  sce = runPCA(sce,exprs_values="logcounts",ncomponents = 20)
  #pca
  dr = reducedDim(sce,"PCA")
  pdf(paste0(savepath,"mus/",mus$design[[i]],mus$data[[i]],"-",mus$norm_method[[i]],".pdf"),width=10,height = 5)
  pca <- ggplot(data=as.data.frame(dr),aes(x=PC1,y=PC2,col=sce$SNG.1ST))+
    geom_point(show.legend = T)+
    scale_size_continuous(range = c(1,4))+
    scale_color_manual(values = tissue_col)+
    labs(x=paste0("PC1 (",format(attr(dr,"percentVar")[1],digits=3),"%)"),y=paste0("PC2 (",format(attr(dr,"percentVar")[2],digits=3),"%)"),col="Labels")+
    theme_bw()+theme(text = element_text(size=20),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"),
                     axis.text.x=element_blank(),
                     axis.ticks.x=element_blank(),
                     axis.text.y=element_blank(),
                     axis.ticks.y=element_blank())+
    ggtitle(paste0(mus$data[[i]],"-",mus$norm_method[[i]]))
  print(pca)
  dev.off()
}
}

##pbmc
designs=c("pbmc5k","pbmc10k")
for (d in designs) {
  readRDS(paste0(sce.path,d,"_norm.rds")) -> pbmc
  pbmc  %>% dplyr::filter(!is.task_error(result))-> pbmc
  if (d=="pbmc5k") {col=pbmc5k_col} else{col=pbmc10k_col}
  for (i in 1:nrow(pbmc)){
    pbmc$result[[i]]->sce
    #sce[,!is.na(sce$SNG.1ST)] ->sce
    #sce[,!sce$SNG.1ST==""]->sce
    as.character(sce$SNG.1ST)->sce$SNG.1ST
    
    sce = runPCA(sce,exprs_values="logcounts",ncomponents = 20)
    #pca
    dr = reducedDim(sce,"PCA")
    pdf(paste0(savepath,"pbmc/",pbmc$design[[i]],pbmc$data[[i]],"-",pbmc$norm_method[[i]],".pdf"),width=7.5,height = 5)
    pca <- ggplot(data=as.data.frame(dr),aes(x=PC1,y=PC2,col=sce$SNG.1ST))+
      geom_point(show.legend = T)+
      scale_size_continuous(range = c(1,4))+
      scale_color_manual(values = col)+
      labs(x=paste0("PC1 (",format(attr(dr,"percentVar")[1],digits=3),"%)"),
           y=paste0("PC2 (",format(attr(dr,"percentVar")[2],digits=3),"%)"),
           col="Labels")+
      theme_bw()+theme(text = element_text(size=20),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"),
                       axis.text.x=element_blank(),
                       axis.ticks.x=element_blank(),
                       axis.text.y=element_blank(),
                       axis.ticks.y=element_blank())+
      ggtitle(paste0(pbmc$data[[i]],"-",pbmc$norm_method[[i]]))
    print(pca)
    dev.off()
  }
}
