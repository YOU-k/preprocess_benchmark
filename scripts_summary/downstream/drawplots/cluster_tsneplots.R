# plate_mixture
getwd()
readRDS("cluster/plate_mixture.rds") -> plate_mixture
readRDS("cluster/plate_mixture_sctransform.rds") -> plate_mixture_sctransform
readRDS("cluster/plate_mixture_seurat.rds") -> plate_mixture_seurat

savepath<- "/stornext/General/data/user_managed/grpu_mritchie_1/Yue/preprocess/new/results/cluster/plots"
bind_rows(plate_mixture %>% filter(!norm_method=="sctransform") %>% select(-hivar_method) ,
          plate_mixture_sctransform, plate_mixture_seurat)  -> plate_mixture
plate_mixture <- plate_mixture[plate_mixture$design=="rnamix",]

plate_mixture<- plate_mixture[!is.task_error(plate_mixture$result),]
library(ggpubr)
for (i in 1:nrow(plate_mixture)){
  plate_mixture[i,]->info
  plate_mixture$result[[i]]->sce
  sce = runPCA(sce,exprs_values="logcounts",ncomponents = 2)
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
#we only plot plate_3cl and plate_5cl 

readRDS("cluster/plate_singlecell.rds") -> plate_singlecell
readRDS("cluster/plate_singlecell_sctransform.rds") -> plate_singlecell_sctransform
readRDS("cluster/plate_singlecell_seurat.rds") -> plate_singlecell_seurat

bind_rows(plate_singlecell %>% filter(!norm_method=="sctransform") %>% select(-hivar_method) ,
          plate_singlecell_sctransform, plate_singlecell_seurat)  -> plate_singlecell

plate_singlecell<- plate_singlecell[!is.task_error(plate_singlecell$result),]

plate_singlecell3 <- plate_singlecell[plate_singlecell$design=="plate3cl",]


cellcolor <- c(rgb(1,0,0,alpha = 0.9),rgb(0,1,0,alpha = 0.9),rgb(0,0,1,alpha = 0.9),"purple","orange","grey")
names(cellcolor) <- c("H1975","HCC827","H2228","H838","A549"," ")

for (i in 1:nrow(plate_singlecell3)){
  plate_singlecell3[i,]->info
  plate_singlecell3$result[[i]]->sce
  sce = runPCA(sce,exprs_values="logcounts",ncomponents = 2)
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
  pdf(paste0(savepath,"/plate_3cl/",info$data,"-",info$norm_method,"-",info$cluster_method,".pdf"),width=10,height = 5)
  print(ggarrange(tsne1,tsne2,ncol=2,widths = c(5,5)))
  dev.off()
}





plate_singlecell5 <- plate_singlecell[plate_singlecell$design=="plate5cl1",]

for (i in 1:nrow(plate_singlecell5)){
  plate_singlecell5[i,]->info
  plate_singlecell5$result[[i]]->sce
  sce = runPCA(sce,exprs_values="logcounts",ncomponents = 2)
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
  pdf(paste0(savepath,"/plate_5cl1/",info$data,"-",info$norm_method,"-",info$cluster_method,".pdf"),width=10,height = 5)
  print(ggarrange(tsne1,tsne2,ncol=2,widths = c(5,5)))
  dev.off()
}
#droplet_singlecell
# we only plate sc3cl and scv35cl


readRDS("cluster/droplet_singlecell.rds") -> droplet_singlecell
readRDS("cluster/droplet_singlecell_sctransform.rds") -> droplet_singlecell_sctransform
readRDS("cluster/droplet_singlecell_seurat.rds") -> droplet_singlecell_seurat

bind_rows(droplet_singlecell %>% filter(!norm_method=="sctransform") %>% select(-hivar_method) ,
          droplet_singlecell_sctransform, droplet_singlecell_seurat)  -> plate_singlecell

droplet_singlecell<- droplet_singlecell[!is.task_error(droplet_singlecell$result),]

droplet_singlecell3 <- droplet_singlecell[droplet_singlecell$design=="10Xv2_3cl",]

for (i in 1:nrow(droplet_singlecell3)){
  droplet_singlecell3[i,]->info
  droplet_singlecell3$result[[i]]->sce
  sce = runPCA(sce,exprs_values="logcounts",ncomponents = 2)
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
  pdf(paste0(savepath,"/sc3cl/",info$data,"-",info$norm_method,"-",info$cluster_method,".pdf"),width=10,height = 5)
  print(ggarrange(tsne1,tsne2,ncol=2,widths = c(5,5)))
  dev.off()
}


droplet_singlecell5 <- droplet_singlecell[droplet_singlecell$design=="10Xv3_5cl",]

for (i in 1:nrow(droplet_singlecell5)){
  droplet_singlecell5[i,]->info
  droplet_singlecell5$result[[i]]->sce
  sce = runPCA(sce,exprs_values="logcounts",ncomponents = 2)
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
  pdf(paste0(savepath,"/sc5clv3/",info$data,"-",info$norm_method,"-",info$cluster_method,".pdf"),width=10,height = 5)
  print(ggarrange(tsne1,tsne2,ncol=2,widths = c(5,5)))
  dev.off()
}
#droplet tissue 
# we only plate tissue2


readRDS("cluster/droplet_tissue.rds") -> droplet_tissue
readRDS("cluster/droplet_tissue_seurat.rds") -> droplet_tissue_seurat
readRDS("cluster/droplet_tissue_seurat_cellranger.rds") -> droplet_tissue_cellranger
readRDS("cluster/droplet_tissue_sctransform.rds") -> droplet_tissue_sctransform


bind_rows(droplet_tissue  %>% select(-hivar_method) ,
          droplet_tissue_sctransform, droplet_tissue_seurat %>% filter(!data=="Cell Ranger"),droplet_tissue_cellranger )  -> droplet_tissue

droplet_tissue<- droplet_tissue[!is.task_error(droplet_tissue$result),]


droplet_tissue$data <- recode(droplet_tissue$data,"scpipe"="scPipe","zumis"="zUMIs","kb"="kallisto bustools",
                       "optimus" ="Optimus","drop"="dropSeqPipe","cellranger"="Cell Ranger",
                       "alevin"="salmon alevin")

droplet_tissue$design <- recode(droplet_tissue$design,"mus1"="10Xv2_tissue1","mus2"="10Xv2_tissue2")


droplet_tissue$result[[1]] ->sce
tissue_color <- c("dark grey","#D3E6B7","#7C71D8" ,"#7A99D7", "#D5E4E2", "#8643E6", "#DF7C44", "#C7BADB", "#7ACADE" ,"#CC6D7F",
                  "#75E29C" ,"#E2BE70" ,"#D595D8" ,"#7BE5D3")
#"#E056A7"
names(tissue_color) <- unique(as.character(sce$SNG.1ST))

tissue_color[1] <- adjustcolor("dark grey", alpha.f = 0.2)


droplet_tissue <- droplet_tissue[droplet_tissue$design=="10Xv2_tissue2",]

library(ggpubr)
#cluser color
library(RColorBrewer)
nb.cols <- 23
mycolors <- colorRampPalette(brewer.pal(8, "Set1"))(nb.cols)

for (i in 1:nrow(droplet_tissue)){
  droplet_tissue[i,]->info
  droplet_tissue$result[[i]]->sce
  as.character(sce$SNG.1ST)->sce$SNG.1ST
  sce$SNG.1ST[(is.na(sce$SNG.1ST))] <- "NA"
  sce$SNG.1ST[sce$SNG.1ST==""] <- "others"
  sce = runPCA(sce,exprs_values="logcounts",ncomponents = 20)
  runTSNE(sce)->sce
  tr = reducedDim(sce,"TSNE")
  tsne1 <- ggplot(data=as.data.frame(tr),aes(x=tr[,1],y=tr[,2],col=sce$SNG.1ST))+
    geom_point(show.legend = T)+
    scale_color_manual(values = tissue_color) +
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
    scale_color_manual(values=mycolors) +
    labs(x="tSNE-1",y="tSNE-2")+
    theme_bw()+
    theme(text = element_text(size=15),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          legend.position = "none",
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank()) +labs(title = paste0(info$data,"-",info$norm_method,"-",info$cluster_method))
  pdf(paste0(savepath,"/tissue2/",info$data,"-",info$norm_method,"-",info$cluster_method,".pdf"),width=10,height = 5)
  print(ggarrange(tsne1,tsne2,ncol=2,widths = c(5,5)))
  dev.off()
}
