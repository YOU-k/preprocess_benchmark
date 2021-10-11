sce.path <- "/stornext/HPCScratch/home/you.y/preprocess_update/SCEs/cluster/"
savepath <- "/stornext/HPCScratch/home/you.y/preprocess_update/results/tsneplots/"
library(scater)
library(scran)

#plate_rnamix
setwd(sce.path)
readRDS("plate_rnamix_all.rds") -> plate_mixture
plate_mixture<- plate_mixture[!is.task_error(plate_mixture$result),]
library(ggpubr)

for (i in 1:nrow(plate_mixture)){
  plate_mixture[i,]->info
  plate_mixture$result[[i]]->sce
  sce = runPCA(sce,exprs_values="logcounts",ncomponents = 20)
  runTSNE(sce)->sce
  col <- rgb(sce$H1975_prop, sce$HCC827_prop, sce$H2228_prop,alpha=0.9)
  tr = reducedDim(sce,"TSNE")
  tsne1 <- ggplot(data=as.data.frame(tr),aes(x=tr[,1],y=tr[,2],col=sce$group))+
    geom_point(show.legend = T)+
    scale_color_manual(guide=FALSE,values = unique(col), limits = unique(sce$group)) +
    labs(x="tSNE-1",y="tSNE-2")+
    theme_bw()+
    theme(text = element_text(size=15),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank()) +labs(title = paste0(info$data,"-",info$norm_method,"-",info$cluster_method))
  as.character(sce$cluster) ->sce$cluster
  tsne2 <- ggplot(data=as.data.frame(tr),aes(x=tr[,1],y=tr[,2],col=sce$cluster))+
    geom_point(show.legend = T)+
    scale_color_brewer(palette = "Set2") +
    labs(x="tSNE-1",y="tSNE-2")+
    theme_bw()+
    theme(text = element_text(size=15),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          legend.position = "none",
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank()) +labs(title = paste0(info$data,"-",info$norm_method,"-",info$cluster_method))
  pdf(paste0(savepath,"/rnamix/",info$data,"-",info$norm_method,"-",info$cluster_method,".pdf"),width=10,height = 5)
  print(ggarrange(tsne1,tsne2,ncol=2,widths = c(5,5)))
  dev.off()
}



#plate_singlecell
readRDS(file.path(sce.path,"plate_cellline_all.rds")) ->plate_singlecell
plate_singlecell %>% filter(!is.task_error(result)) ->plate_singlecell

cellcolor <- c(rgb(1,0,0,alpha = 0.9),rgb(0,1,0,alpha = 0.9),rgb(0,0,1,alpha = 0.9),"purple","orange","grey")
names(cellcolor) <- c("H1975","HCC827","H2228","H838","A549"," ")

for (i in 1:nrow(plate_singlecell)){
  plate_singlecell[i,]->info
  plate_singlecell$result[[i]]->sce
  sce = runPCA(sce,exprs_values="logcounts",ncomponents = 20)
  runTSNE(sce)->sce
  tr = reducedDim(sce,"TSNE")
  as.character(sce$SNG.1ST) -> sce$SNG.1ST
  sce$SNG.1ST[is.na(sce$SNG.1ST)] <- " "
  tsne1 <- ggplot(data=as.data.frame(tr),aes(x=tr[,1],y=tr[,2],col=sce$SNG.1ST))+
    geom_point(show.legend = T)+
    scale_color_manual(values = cellcolor) +
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
    geom_point(show.legend = T)+
    scale_color_brewer(palette = "Set2") +
    labs(x="tSNE-1",y="tSNE-2")+
    theme_bw()+
    theme(text = element_text(size=15),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          legend.position = "none",
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank()) +labs(title = paste0(info$data,"-",info$norm_method,"-",info$cluster_method))
  pdf(paste0(savepath,"cellline_p/",info$design,info$data,"-",info$norm_method,"-",info$cluster_method,".pdf"),width=10,height = 5)
  print(ggarrange(tsne1,tsne2,ncol=2,widths = c(5,5)))
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
