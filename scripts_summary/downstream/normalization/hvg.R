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

hivar_method = list(scran_hi = scran_high_var)

#plate_mixture
plate_mixture <- readRDS(file.path(read.path,"plate_mixture.rds"))
plate_mixture$data <- recode(plate_mixture$data,"celseq"="celseq2", "scpipe"="scPipe","zumis"="zUMIs")
plate_mixture %>% arrange(design,data,norm_method) %>% apply_methods(hivar_method) -> plate_mixture_hvged
saveRDS(plate_mixture_hvged,file.path(read.path,"hvged/plate_mixture_hvged.rds"))

#plate_singlecell
plate_singlecell <- readRDS(file.path(read.path,"plate_singlecell.rds"))
plate_singlecell$data <- recode(plate_singlecell$data,"celseq"="celseq2", "scpipe"="scPipe","zumis"="zUMIs")

plate_singlecell %>% arrange(design,data,norm_method) %>% apply_methods(hivar_method) -> plate_singlecell_hvged
saveRDS(plate_singlecell_hvged,file.path(read.path,"hvged/plate_singlecell_hvged.rds"))


#droplet_singlecell 

droplet_singlecell <- readRDS(file.path(read.path,"droplet_singlecell.rds"))
droplet_singlecell$design <-
  recode(droplet_singlecell$design,"sc3cl"="10Xv2_3cl","sc5clv2"="10Xv2_5cl","sc5clv3"="10Xv3_5cl",
         "mus1"="10Xv2_tissue1","mus2"="10Xv2_tissue2")
droplet_singlecell$data <- recode(droplet_singlecell$data, "scpipe"="scPipe","zumis"="zUMIs","kb"="kallisto bustools",
                                  "optimus" ="Optimus","drop"="dropSeqPipe","cellranger"="Cell Ranger",
                                  "alevin"="salmon alevin")
droplet_singlecell %>% arrange(design,data,norm_method) %>% apply_methods(hivar_method) -> droplet_singlecell_hvged
saveRDS(droplet_singlecell_hvged,file.path(read.path,"hvged/droplet_singlecell_hvged.rds"))


#droplet_tissue
droplet_tissue <- readRDS(file.path(read.path,"droplet_tissue.rds"))
droplet_tissue$design <-
  recode(droplet_tissue$design,"sc3cl"="10Xv2_3cl","sc5clv2"="10Xv2_5cl","sc5clv3"="10Xv3_5cl",
         "mus1"="10Xv2_tissue1","mus2"="10Xv2_tissue2")
droplet_tissue$data <- recode(droplet_tissue$data, "scpipe"="scPipe","zumis"="zUMIs","kb"="kallisto bustools",
                              "optimus" ="Optimus","drop"="dropSeqPipe","cellranger"="Cell Ranger",
                              "alevin"="salmon alevin")
droplet_tissue %>% arrange(design,data,norm_method) %>% apply_methods(hivar_method) -> droplet_tissue_hvged
saveRDS(droplet_tissue_hvged,file.path(read.path,"hvged/droplet_tissue_hvged.rds"))





#droplet_tissue cellranger
droplet_tissue <- readRDS(file.path(read.path,"droplet_tissue_cellranger.rds"))
droplet_tissue$design <-
  recode(droplet_tissue$design,"sc3cl"="10Xv2_3cl","sc5clv2"="10Xv2_5cl","sc5clv3"="10Xv3_5cl",
         "mus1"="10Xv2_tissue1","mus2"="10Xv2_tissue2")
droplet_tissue$data <- recode(droplet_tissue$data, "scpipe"="scPipe","zumis"="zUMIs","kb"="kallisto bustools",
                              "optimus" ="Optimus","drop"="dropSeqPipe","cellranger"="Cell Ranger",
                              "alevin"="salmon alevin")
droplet_tissue %>% arrange(design,data,norm_method) %>% apply_methods(hivar_method) -> droplet_tissue_hvged
saveRDS(droplet_tissue_hvged,file.path(read.path,"hvged/droplet_tissue_hvged_cellranger.rds"))
