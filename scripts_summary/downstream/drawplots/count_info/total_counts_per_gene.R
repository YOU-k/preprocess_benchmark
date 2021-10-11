SCE_drop.path <- "/stornext/HPCScratch/home/you.y/preprocess_update/SCEs/raw/droplet_based/"
write.path <- "/stornext/HPCScratch/home/you.y/preprocess_update/results/droplet-based/density"

library(ggpubr)
cal_total_counts_per_gene <- function(sce){
  return(apply(as.matrix(counts(sce)),1,sum))
}

droplet_recode <- function(vec){
  return(recode(vec, "scpipe"="scPipe","zumis"="zUMIs","kb"="kallisto bustools",
                "optimus" ="Optimus","drop"="dropSeqPipe","cellranger"="Cell Ranger",
                "sa"="salmon_SA","splici"="salmon_splici"))
}



dir(SCE_drop.path)
#"mus1","mus2",
designs=c("sc3cl","sc5clv3","sc5cl","pbmc5k","pbmc10k")
designs=c("sc5cl")
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
  
  
  
  library(SingleCellExperiment)
  
  lapply(datasets$result,cal_total_counts_per_gene) -> tpg
  lapply(datasets$result,function(sce){rowData(sce)$gene_biotype}) ->biotype
  
  data.frame(tpg=log10(unlist(tpg)),biotype=unlist(biotype),
             Preprocess=rep(datasets$data,unlist(lapply(tpg, length)))) -> df
  
  df$biotype2 <- df$biotype
  df$biotype2[grepl("pseudo",df$biotype)] <- "pseudogene"
  df$biotype2[!df$biotype2 %in% c("lncRNA","pseudogene","protein_coding")] <- "others"
  df$biotype2<- factor(df$biotype2,levels=c("protein_coding","lncRNA","pseudogene","others"))
  pdf(paste0(write.path,"/",design,"_facet.pdf"),width = 9,height = 3)
  p1<- ggdensity(df, x = "tpg",xlab = "Total counts per gene (log10)",ylab="Density",
             rug = FALSE, facet.by = "biotype2",
            color = "Preprocess",nrow=1, legend = "right",
            palette = droplet_col,size = 0.8)
  print(p1)
  dev.off()
  
  pdf(paste0(write.path,"/",design,".pdf"),width = 4.5,height = 3)
  p <-ggdensity(df, x = "tpg",xlab = "Total counts per gene (log10)",ylab="Density",
            rug = FALSE,
            color = "Preprocess", legend = "right",
            palette = droplet_col,size = 0.8)
  print(p)
  dev.off()
  
  pdf(paste0(write.path,"/",design,"_by_method.pdf"),width = 6,height = 4.5)
  p <-ggdensity(df, x = "tpg",xlab = "Total counts per gene (log10)",ylab="Density",
                rug = FALSE,facet.by = "Preprocess",
                color = "Preprocess", legend = "right",
                palette = droplet_col,size = 0.8)
  print(p)
  dev.off()
  
  Reduce(intersect,lapply(datasets$result,rownames)) -> cg
  
  datasets$result_filter <- lapply(datasets$result,function(sce){sce[cg,]})
  lapply(datasets$result_filter,cal_total_counts_per_gene) -> tpg2
  lapply(datasets$result_filter,function(sce){rowData(sce)$gene_biotype}) ->biotype2
  
  data.frame(tpg=log10(unlist(tpg2)),biotype=unlist(biotype2),
             Preprocess=rep(datasets$data,unlist(lapply(tpg2, length)))) -> df2
  
  pdf(paste0(write.path,"/",design,"_filtered.pdf"),width = 4.5,height = 3)
  p2 <-ggdensity(df2, x = "tpg",xlab = "Total counts per gene (log10)",ylab="Density",
                rug = FALSE,
                color = "Preprocess", legend = "right",
                palette = droplet_col,size = 0.8)
  print(p2)
  dev.off()
  
  pdf(paste0(write.path,"/",design,"_filtered_by_method.pdf"),width = 6,height = 4.5)
  p <-ggdensity(df2, x = "tpg",xlab = "Total counts per gene (log10)",ylab="Density",
                rug = FALSE,facet.by = "Preprocess",
                color = "Preprocess", legend = "right",
                palette = droplet_col,size = 0.8)
  print(p)
  dev.off()
  
  df2$biotype2 <- df2$biotype
  df2$biotype2[grepl("pseudo",df2$biotype)] <- "pseudogene"
  df2$biotype2[!df2$biotype2 %in% c("lncRNA","pseudogene","protein_coding")] <- "others"
  df2$biotype2<- factor(df2$biotype2,levels=c("protein_coding","lncRNA","pseudogene","others"))
  pdf(paste0(write.path,"/",design,"_facet_filtered.pdf"),width = 9,height = 3)
  p3<- ggdensity(df2, x = "tpg",xlab = "Total counts per gene (log10)",ylab="Density",
                 rug = FALSE, facet.by = "biotype2",
                 color = "Preprocess",nrow=1, legend = "right",
                 palette = droplet_col,size = 0.8)
  print(p3)
  dev.off()
  
}




#for plate_3cl
plate_recode <- function(vec){
  return(recode(vec, "scpipe"="scPipe","zumis"="zUMIs","kb"="kallisto bustools",
                "celseq"="celseq2"))
}
design=d="plate3cl"
dataname=c("scpipe","zumis","kb","celseq","scruff")

paste(apply(expand.grid(design, dataname), 1, paste, collapse="_"), ".rds",sep="") ->files
SCE_plate.path ="/stornext/HPCScratch/home/you.y/preprocess_update/SCEs/raw/plate_based/"
files %in% dir(SCE_plate.path)
files <- files[files %in% dir(SCE_plate.path)]
result0 <- lapply(files, function(file){readRDS(file.path(SCE_plate.path,file))})

rep(dataname,each=length(design))->data

library(tidyverse)
tibble(design=rep(design,length(dataname)), 
       data=plate_recode(data),
       result=result0) -> datasets



library(SingleCellExperiment)

lapply(datasets$result,cal_total_counts_per_gene) -> tpg
lapply(datasets$result,function(sce){rowData(sce)$gene_biotype}) ->biotype

data.frame(tpg=log10(unlist(tpg)),biotype=unlist(biotype),
           Preprocess=rep(datasets$data,unlist(lapply(tpg, length)))) -> df

df$biotype2 <- df$biotype
df$biotype2[grepl("pseudo",df$biotype)] <- "pseudogene"
df$biotype2[!df$biotype2 %in% c("lncRNA","pseudogene","protein_coding")] <- "others"
df$biotype2<- factor(df$biotype2,levels=c("protein_coding","lncRNA","pseudogene","others"))
pdf(paste0(write.path,"/",design,"_facet.pdf"),width = 9,height = 3)
p1<- ggdensity(df, x = "tpg",xlab = "Total counts per gene (log10)",ylab="Density",
               rug = FALSE, facet.by = "biotype2",
               color = "Preprocess",nrow=1, legend = "right",
               palette = plate_col_k,size = 0.8)
print(p1)
dev.off()

pdf(paste0(write.path,"/",design,".pdf"),width = 4.5,height = 3)
p <-ggdensity(df, x = "tpg",xlab = "Total counts per gene (log10)",ylab="Density",
              rug = FALSE,
              color = "Preprocess", legend = "right",
              palette = plate_col_k,size = 0.8)
print(p)
dev.off()

pdf(paste0(write.path,"/",design,"_by_method.pdf"),width = 6,height = 4.5)
p <-ggdensity(df, x = "tpg",xlab = "Total counts per gene (log10)",ylab="Density",
              rug = FALSE,facet.by = "Preprocess",
              color = "Preprocess", legend = "right",
              palette = plate_col_k,size = 0.8)
print(p)
dev.off()
