sce.path <- "/stornext/General/data/user_managed/grpu_mritchie_1/Yue/preprocess/new/SCE/cluster"
savepath <- "/stornext/General/data/user_managed/grpu_mritchie_1/Yue/preprocess/new/results/cluster/plots/"
library(scater)
library(scran)
#plate_singlecell
readRDS(file.path(sce.path,"plate_singlecell.rds")) ->plate_singlecell
plate_singlecell %>% filter(!is.task_error(result)) %>% select(-hivar_method) ->plate_singlecell

cellcolor <- c(rgb(1,0,0,alpha = 0.9),rgb(0,1,0,alpha = 0.9),rgb(0,0,1,alpha = 0.9),"purple","orange","grey")
names(cellcolor) <- c("H1975","HCC827","H2228","H8383","A549"," ")

for (i in 1:nrow(plate_singlecell)){
  plate_singlecell$result[[i]]->sce
  sce = runPCA(sce,exprs_values="logcounts",ncomponents = 2)
  #pca
  dr = reducedDim(sce,"PCA")
  pdf(paste0(savepath,"simplemix/cluster/",plate_singlecell$design[[i]],plate_singlecell$data[[i]],"-",plate_singlecell$norm_method[[i]],plate_singlecell$cluster_method[[i]],".pdf"),width=5,height = 5)
  as.factor(sce$cluster) ->sce$cluster
  pca2 <- ggplot(data=as.data.frame(dr),aes(x=PC1,y=PC2,col=sce$cluster))+
    geom_point(show.legend = F)+
    scale_size_continuous(range = c(1,4))+
    labs(x=paste0("PC1 (",format(attr(dr,"percentVar")[1],digits=3),"%)"),y=paste0("PC2 (",format(attr(dr,"percentVar")[2],digits=3),"%)"))+
    theme_bw()+theme(text = element_text(size=20),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"),
                     axis.text.x=element_blank(),
                     axis.ticks.x=element_blank(),
                     axis.text.y=element_blank(),
                     axis.ticks.y=element_blank())+
    ggtitle(paste0(plate_singlecell$cluster_method[[i]],"-",plate_singlecell$norm_method[[i]]))
  print(pca2)
  dev.off()
  
  pdf(paste0(savepath,"simplemix/cell/",plate_singlecell$design[[i]],plate_singlecell$data[[i]],"-",plate_singlecell$norm_method[[i]],plate_singlecell$cluster_method[[i]],".pdf"),width=5,height = 5)
  pca2 <- ggplot(data=as.data.frame(dr),aes(x=PC1,y=PC2,col=sce$SNG.1ST))+
    geom_point(show.legend = F)+
    scale_size_continuous(range = c(1,4))+
    scale_color_manual(values = cellcolor)+
    labs(x=paste0("PC1 (",format(attr(dr,"percentVar")[1],digits=3),"%)"),y=paste0("PC2 (",format(attr(dr,"percentVar")[2],digits=3),"%)"))+
    theme_bw()+theme(text = element_text(size=20),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"),
                     axis.text.x=element_blank(),
                     axis.ticks.x=element_blank(),
                     axis.text.y=element_blank(),
                     axis.ticks.y=element_blank())+
    ggtitle(paste0(plate_singlecell$cluster_method[[i]],"-",plate_singlecell$norm_method[[i]]))
  print(pca2)
  dev.off()
}


#plate_singlecell_seu
readRDS(file.path(sce.path,"plate_singlecell_seurat.rds")) ->plate_singlecell
plate_singlecell %>% filter(!is.task_error(result)) ->plate_singlecell

cellcolor <- c(rgb(1,0,0,alpha = 0.9),rgb(0,1,0,alpha = 0.9),rgb(0,0,1,alpha = 0.9),"purple","orange","grey")
names(cellcolor) <- c("H1975","HCC827","H2228","H8383","A549"," ")

for (i in 1:nrow(plate_singlecell)){
  plate_singlecell$result[[i]]->sce
  sce = runPCA(sce,exprs_values="logcounts",ncomponents = 2)
  #pca
  dr = reducedDim(sce,"PCA")
  pdf(paste0(savepath,"simplemix/cluster/",plate_singlecell$design[[i]],plate_singlecell$data[[i]],"-",plate_singlecell$cluster_method[[i]],".pdf"),width=5,height = 5)
  as.factor(sce$cluster) ->sce$cluster
  pca2 <- ggplot(data=as.data.frame(dr),aes(x=PC1,y=PC2,col=sce$cluster))+
    geom_point(show.legend = F)+
    scale_size_continuous(range = c(1,4))+
    labs(x=paste0("PC1 (",format(attr(dr,"percentVar")[1],digits=3),"%)"),y=paste0("PC2 (",format(attr(dr,"percentVar")[2],digits=3),"%)"))+
    theme_bw()+theme(text = element_text(size=20),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"),
                     axis.text.x=element_blank(),
                     axis.ticks.x=element_blank(),
                     axis.text.y=element_blank(),
                     axis.ticks.y=element_blank())+
    ggtitle(plate_singlecell$cluster_method[[i]])
  print(pca2)
  dev.off()
  
  pdf(paste0(savepath,"simplemix/cell/",plate_singlecell$design[[i]],plate_singlecell$data[[i]],"-",plate_singlecell$cluster_method[[i]],".pdf"),width=5,height = 5)
  pca2 <- ggplot(data=as.data.frame(dr),aes(x=PC1,y=PC2,col=sce$SNG.1ST))+
    geom_point(show.legend = F)+
    scale_size_continuous(range = c(1,4))+
    scale_color_manual(values = cellcolor)+
    labs(x=paste0("PC1 (",format(attr(dr,"percentVar")[1],digits=3),"%)"),y=paste0("PC2 (",format(attr(dr,"percentVar")[2],digits=3),"%)"))+
    theme_bw()+theme(text = element_text(size=20),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"),
                     axis.text.x=element_blank(),
                     axis.ticks.x=element_blank(),
                     axis.text.y=element_blank(),
                     axis.ticks.y=element_blank())+
    ggtitle(plate_singlecell$cluster_method[[i]])
  print(pca2)
  dev.off()
}


#droplet_singlecell
readRDS(file.path(sce.path,"droplet_singlecell.rds")) ->droplet_singlecell
droplet_singlecell %>% filter(!is.task_error(result)) %>% select(-hivar_method) ->droplet_singlecell

cellcolor <- c(rgb(1,0,0,alpha = 0.9),rgb(0,1,0,alpha = 0.9),rgb(0,0,1,alpha = 0.9),"purple","orange","grey")
names(cellcolor) <- c("H1975","HCC827","H2228","H8383","A549"," ")

for (i in 1:nrow(droplet_singlecell)){
  droplet_singlecell$result[[i]]->sce
  sce = runPCA(sce,exprs_values="logcounts",ncomponents = 2)
  #pca
  dr = reducedDim(sce,"PCA")
  pdf(paste0(savepath,"droplet/cluster/",droplet_singlecell$design[[i]],droplet_singlecell$data[[i]],"-",droplet_singlecell$norm_method[[i]],droplet_singlecell$cluster_method[[i]],".pdf"),width=5,height = 5)
  as.factor(sce$cluster) ->sce$cluster
  pca2 <- ggplot(data=as.data.frame(dr),aes(x=PC1,y=PC2,col=sce$cluster))+
    geom_point(show.legend = F)+
    scale_size_continuous(range = c(1,4))+
    labs(x=paste0("PC1 (",format(attr(dr,"percentVar")[1],digits=3),"%)"),y=paste0("PC2 (",format(attr(dr,"percentVar")[2],digits=3),"%)"))+
    theme_bw()+theme(text = element_text(size=20),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"),
                     axis.text.x=element_blank(),
                     axis.ticks.x=element_blank(),
                     axis.text.y=element_blank(),
                     axis.ticks.y=element_blank())+
    ggtitle(paste0(droplet_singlecell$cluster_method[[i]],"-",droplet_singlecell$norm_method[[i]]))
  print(pca2)
  dev.off()
  
  pdf(paste0(savepath,"droplet/cell/",droplet_singlecell$design[[i]],droplet_singlecell$data[[i]],"-",droplet_singlecell$norm_method[[i]],droplet_singlecell$cluster_method[[i]],".pdf"),width=5,height = 5)
  pca2 <- ggplot(data=as.data.frame(dr),aes(x=PC1,y=PC2,col=sce$SNG.1ST))+
    geom_point(show.legend = F)+
    scale_size_continuous(range = c(1,4))+
    scale_color_manual(values = cellcolor)+
    labs(x=paste0("PC1 (",format(attr(dr,"percentVar")[1],digits=3),"%)"),y=paste0("PC2 (",format(attr(dr,"percentVar")[2],digits=3),"%)"))+
    theme_bw()+theme(text = element_text(size=20),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"),
                     axis.text.x=element_blank(),
                     axis.ticks.x=element_blank(),
                     axis.text.y=element_blank(),
                     axis.ticks.y=element_blank())+
    ggtitle(paste0(droplet_singlecell$cluster_method[[i]],"-",droplet_singlecell$norm_method[[i]]))
  print(pca2)
  dev.off()
}
