SCE_drop.path <- "/stornext/HPCScratch/home/you.y/preprocess_update/SCEs/raw/droplet_based/"
write.path <- "/stornext/HPCScratch/home/you.y/preprocess_update/results/droplet-based/"
get_barcodes <- function(sce){
  colnames(sce)[sce$keep]
}
get_overlap_counts <- function(sce){
  sce[rownames(sce) %in% allgenes,colnames(sce) %in% allbc]  -> sce
  sce[match(allgenes,rownames(sce)),match(allbc,colnames(sce))] ->sce
  return(as.matrix(counts(sce)))
}
droplet_recode <- function(vec){
  return(recode(vec, "scpipe"="scPipe","zumis"="zUMIs","kb"="kallisto bustools",
                "optimus" ="Optimus","drop"="dropSeqPipe","cellranger"="Cell Ranger",
                "sa"="salmon_SA","splici"="salmon_splici"))
}



dir(SCE_drop.path)
#"mus1","mus2","sc3cl","sc5clv3"
designs=c("pbmc5k","pbmc10k","sc5cl")
designs=c("sc5clv3")
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
  
  
  rep(dataname,each=length(design))->data
  
  library(tidyverse)
  tibble(design=rep(design,length(dataname)), 
         data=droplet_recode(data),
         result=result0) -> datasets
  datasets$design <- recode(datasets$design,"mus1"="10xv2_lungtissue1","mus2"="10xv2_lungtissue2","sc3cl"="10xv2_3cellline",
                            "sc5cl"="10xv2_5cellline","sc5clv3"="10xv3_5cellline","pbmc5k"="10xv3_pbmc5k",
                            "pbmc10k"="10xv3_pbmc10k")
  
  
  library(SingleCellExperiment)
  #total and detected
  lapply(datasets$result,function(sce){sce$total}) -> total
  lapply(datasets$result,function(sce){sce$detected}) ->detected
  lapply(datasets$result,function(sce){!is.na(sce$SNG.1ST)}) ->label
  
  
  data.frame(total=log10(unlist(total)),detected=log10(unlist(detected)),
             Preprocess=rep(datasets$data,unlist(lapply(total, length))),
             label=unlist(label)) -> df
 

  
  pdf(paste0(write.path,"total_detected/",d,"_total_counts.pdf"),width = 6,height =3)
  p_t <- ggplot(df) +geom_violin(aes(x=reorder(Preprocess,total),y=total,col=Preprocess)) +
    scale_color_manual(values=c(droplet_col))+
    labs(x="",y="Total counts per cell (log10)",col="Preprocessing methods",title = unique(datasets$design)) +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),text = element_text(size=14),
          axis.ticks.x=element_blank(),panel.background = element_rect(fill = "white", colour = NA), 
          panel.border = element_rect(fill = NA,colour = "grey20"), panel.grid = element_line(colour = "grey92"), 
          panel.grid.minor = element_line(size = rel(0.5)), 
          strip.background = element_rect(fill = "grey85",colour = "grey20"))
  print(p_t)
  dev.off()
  
  pdf(paste0(write.path,"total_detected/",d,"_detected.pdf"),width = 6.2,height =3)
  p_d <-ggplot(df) + geom_violin(aes(x=reorder(Preprocess,detected),y=detected,col=Preprocess)) +
    scale_color_manual(values=droplet_col) +
    labs(x="",y="Number of features detected \nper cell (log10)",col="Preprocessing methods",title = unique(datasets$design)) +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),text = element_text(size=14),
          axis.ticks.x=element_blank(),panel.background = element_rect(fill = "white", colour = NA), 
          panel.border = element_rect(fill = NA,colour = "grey20"), panel.grid = element_line(colour = "grey92"), 
          panel.grid.minor = element_line(size = rel(0.5)), 
          strip.background = element_rect(fill = "grey85",colour = "grey20"))
  print(p_d)
  dev.off()
  
  #cellnum
  datasets$cellnum <- lapply(datasets$result, ncol)
  pdf(paste0(write.path,"cellnum/",d,"_cellnum.pdf"),width = 4,height =3)
  p_cn <- datasets %>% 
    mutate(cellnum=as.numeric(cellnum)) %>% 
    ggplot() + geom_col(aes(x=reorder(data,cellnum),y=cellnum,fill=data))+
    scale_fill_manual(values=droplet_col) +
    labs(y="Number of cells", x="",fill="Preprocess methods",title = unique(datasets$design))+
    coord_flip()+
    theme(text = element_text(size=14),
          axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          panel.background = element_rect(fill = "white", colour = NA), 
          panel.border = element_rect(fill = NA,colour = "grey20"), panel.grid = element_line(colour = "grey92"), 
          panel.grid.minor = element_line(size = rel(0.5)), 
          strip.background = element_rect(fill = "grey85",colour = "grey20"))
  print(p_cn)
  dev.off()
  
  # upsetR
  lapply(datasets$result, get_barcodes) -> bc_all
  library(UpSetR)
  names(bc_all) <- datasets$data
  pdf(paste0(write.path,"upsetr/",d,"_upsetr.pdf"),width = 7,height =6.5)
  p_upset <- upset(fromList(bc_all),nsets = nrow(datasets),order.by="freq",nintersects = 7, 
        sets.bar.color = "#56B4E9",matrix.color = "blue",text.scale = 2,)
  print(p_upset)
  dev.off()
  
  #cor
  Reduce(intersect, lapply(datasets$result, colnames)) -> allbc
  Reduce(intersect, lapply(datasets$result, rownames)) -> allgenes
  
  lapply(datasets$result, get_overlap_counts) -> counts_all
  names(counts_all) <- datasets$data
  
  calculate_cor <- function(pp){
    cor = lapply((1+pp):length(counts_all), 
                 function (x){diag(cor(counts_all[[pp]],counts_all[[x]],method = "pearson"))})
    return(unlist(cor))}
  
  cors=list()
  
  for (n in 1:(length(counts_all)-1)){
    cors[[n]] <- data.frame(preprocess1 = rep(names(counts_all)[(n+1):length(counts_all)],each=length(allbc)),
                            preprocess2 = rep(names(counts_all)[n], length(allbc)*(length(counts_all)-n)),
                            correlation = calculate_cor(n))
  }
  
  Reduce(rbind,cors) ->cors_use
  
  saveRDS(cors_use,paste0("/stornext/HPCScratch/home/you.y/preprocess_update/SCEs/raw/summary/",d,"_cor_use.rds"))
  

  full_join(cors_use  %>% group_by(preprocess1,preprocess2) %>% summarise(median=median(correlation)), cors_use) -> tmp_cor
  factor(tmp_cor$preprocess2,levels = datasets$data) -> tmp_cor$preprocess2
  factor(tmp_cor$preprocess1,levels = datasets$data) -> tmp_cor$preprocess1
  pdf(paste0(write.path,"cor/",d,"_cor.pdf"),width=8, height =8)
  p_cor<-  ggplot(tmp_cor) + geom_boxplot(aes(y=correlation)) + facet_grid(preprocess1~preprocess2) +
    geom_text(aes(0,0.8,label=round(median,2)),size=4,check_overlap = TRUE,col="deeppink3") + 
    labs(y="Pearson correlation coefficient",x="",title=unique(datasets$design))+
      theme(axis.title.x=element_blank(),
            axis.text.x=element_blank(),text = element_text(size=14),
            axis.ticks.x=element_blank(),panel.background = element_rect(fill = "white", colour = NA), 
            panel.border = element_rect(fill = NA,colour = "grey20"), panel.grid = element_line(colour = "grey92"), 
            panel.grid.minor = element_line(size = rel(0.5)), 
            strip.background = element_rect(fill = "grey85",colour = "grey20"))
  print(p_cor)
  dev.off()
  
  
  }




#this is specially for mus1
col=c("lightblue","darkblue")
names(col) <- c(TRUE,FALSE)
pdf(paste0(write.path,"total_detected/mus1_before_filter_total_counts.pdf"),width = 6,height =3)
ggplot(df) + geom_violin(aes(x=reorder(Preprocess,total),y=total,col=Preprocess)) +
  geom_jitter(aes(x=reorder(Preprocess,total),y=total,col=label),alpha=0.1)+
  scale_color_manual(values=c(droplet_col,col))+
  labs(x="",y="Total counts per cell (log10)",col="Preprocessing methods")+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),text = element_text(size=14),
        axis.ticks.x=element_blank(),panel.background = element_rect(fill = "white", colour = NA), 
        panel.border = element_rect(fill = NA,colour = "grey20"), panel.grid = element_line(colour = "grey92"), 
        panel.grid.minor = element_line(size = rel(0.5)), 
        strip.background = element_rect(fill = "grey85",colour = "grey20"))
dev.off()


