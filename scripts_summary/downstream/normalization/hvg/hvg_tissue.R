library(scran)
library(scater)
library(CellBench)
library(tidyverse)

read.path<- "/stornext/General/data/user_managed/grpu_mritchie_1/Yue/preprocess/new/SCE/norm"

scran_high_var = function(sce,topn=1500){
  if (nrow(sce)>1500){
    dec <- modelGeneVar(sce)
    # Ordering by most interesting genes for inspection.
    dec[order(dec$bio, decreasing=TRUE),] 
    chosen <- getTopHVGs(dec, n=topn)
    sce <- sce[chosen,]
  }
  return(sce)
}

get_hvg_curve = function(sce,topn=1500){
  if (nrow(sce)>1500){
    dec <- modelGeneVar(sce)
    plot(dec$mean, dec$total, pch=16, xlab="Mean of log-expression",
         ylab="Variance of log-expression")
    data.frame(curve(metadata(dec)$trend(x), col="dodgerblue",add = TRUE))->db
  }
  return(db)
}
hivar_method = list(scran_hi = scran_high_var)

droplet_tissue <- readRDS(file.path(read.path,"droplet_tissue.rds"))
droplet_tissue %>% dplyr::filter(design=="mus2") ->droplet_tissue
droplet_tissue_c <- tibble()
for (d in unique(droplet_tissue$design)){
  droplet_tissue %>% dplyr::filter(design==d) -> tmp1
  Reduce(intersect,lapply(tmp1$result,colnames)) -> commoncells
  lapply(tmp1$result,function(sce){sce[,commoncells]}) ->tmp2
  tmp1$result <- tmp2
  bind_rows(droplet_tissue_c,tmp1) ->droplet_tissue_c
}


#droplet_tissue

droplet_tissue_c$design <-
  recode(droplet_tissue_c$design,"sc3cl"="10Xv2_3cl","sc5clv2"="10Xv2_5cl","sc5clv3"="10Xv3_5cl",
         "mus1"="10Xv2_tissue1","mus2"="10Xv2_tissue2")
droplet_tissue_c$data <- recode(droplet_tissue_c$data, "scpipe"="scPipe","zumis"="zUMIs","kb"="kallisto bustools",
                              "optimus" ="Optimus","drop"="dropSeqPipe","cellranger"="Cell Ranger",
                              "alevin"="salmon alevin")

gethvg <- list(hvg=get_hvg_curve)
droplet_tissue %>% 
  dplyr::filter(design=="mus2") %>% 
  dplyr::filter(norm_method=="scone") %>% 
  arrange(design,data,norm_method) %>% apply_methods(gethvg) -> hvgdb
bind_rows(hvgdb$result) -> hvgp
hvgp$preprocess <- rep(plothvg$data,each=101)
ggplot(hvgp) + geom_line(aes(x=x,y=y,col=preprocess),size=1) + 
  labs(x="Mean of log-expression",y="Variance of log-expression",col="Preprocessing workflows") +
  scale_color_manual(values=mycol) +theme_bw()

ggsave("/stornext/General/data/user_managed/grpu_mritchie_1/Yue/preprocess/new/results/norm/hvged_tissue_c_all.pdf",
       width = 5,height = 2.5)

droplet_tissue %>% arrange(design,data,norm_method) %>% apply_methods(hivar_method) -> droplet_tissue_hvged
droplet_tissue_hvged$hivar_method <- NULL
saveRDS(droplet_tissue_hvged,file.path(read.path,"hvged/droplet_tissue_hvged_commoncells.rds"))



#droplet_tissue_only protein coding genes
select_protein_coding_genes <- function(sce){
  sce<- sce[!is.na(rowData(sce)$gene_biotype),]
  sce <- sce[rowData(sce)$gene_biotype=="protein_coding",]
  return(sce)
}
lapply(droplet_tissue_c$result,select_protein_coding_genes) -> dtp
droplet_tissue_c -> droplet_tissue_cp
droplet_tissue_cp$result <- dtp
rm(dtp)
droplet_tissue_cp %>% arrange(design,data,norm_method) %>% apply_methods(hivar_method) -> droplet_tissue_hvged
droplet_tissue_hvged$hivar_method <- NULL
saveRDS(droplet_tissue_hvged,file.path(read.path,"hvged/droplet_tissue_hvged_commoncells_protein.rds"))


