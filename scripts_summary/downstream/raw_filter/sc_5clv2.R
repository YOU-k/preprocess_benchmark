library(scran)
library(scater)
library(CellBench)
library(R.utils)
library(tidyverse)

setwd("/stornext/HPCScratch/home/you.y/preprocess_update/SCEs/")
data.path <- "/stornext/General/data/user_managed/grpu_mritchie_1/Yue/preprocess/data"
write.path <- "/stornext/HPCScratch/home/you.y/preprocess_update/SCEs/raw/droplet_based/"


barcode_info <- read.table("/stornext/General/data/user_managed/grpu_mritchie_1/SCmixology/scbench_5cellline_10x/scPipe_output/Lib90.best",header = TRUE)

sapply(strsplit(as.character(barcode_info$BEST),"-"),function(x) x[1]) -> barcode_info$doublet
sapply(strsplit(as.character(barcode_info$BARCODE),"-1"),"[",1)->barcode_info$barcode

library(AnnotationHub)
ah <- AnnotationHub()
EnsDb.Hsapiens.v98 <- query(ah, c("EnsDb", "Homo Sapiens", 98))[[1]]

gene_filter <- function(sce){
  keep1 = (apply(counts(sce), 1, function(x) mean(x[x>0])) > 1)  # average count larger than 1
  keep2 = (rowSums(counts(sce)>0) > 10)  
  table(keep1&keep2)
  sce = sce[(keep1 & keep2), ]
  return(sce)
}


calculateQCMet <- function(sce){
  location <- mapIds(
    x = EnsDb.Hsapiens.v98, 
    # NOTE: Need to remove gene version number prior to lookup.
    keys = rownames(sce),
    keytype = "GENEID",
    column = "SEQNAME")
  rowData(sce)$CHR =location
  
  symbol <- mapIds(
    x = EnsDb.Hsapiens.v98, 
    # NOTE: Need to remove gene version number prior to lookup.
    keys = rownames(sce),
    keytype = "GENEID",
    column = "SYMBOL")
  rowData(sce)$symbol =symbol
  
  biotype<- mapIds(
    x = EnsDb.Hsapiens.v98, 
    # NOTE: Need to remove gene version number prior to lookup.
    keys = rownames(sce),
    keytype = "GENEID",
    column = "GENEBIOTYPE")
  rowData(sce)$gene_biotype =biotype
  
  is_mito <- rowData(sce)$CHR == "MT"
  sce <- addPerCellQC(sce, subsets = list(Mito=which(is_mito)))
  return(sce)
}

scater_filter <- function(sce){
  libsize_drop <- isOutlier(
    metric = sce$sum, 
    nmads = 3,
    type = "lower", 
    log = TRUE)
  feature_drop <- isOutlier(
    metric = sce$detected,
    nmads = 3, 
    type = "lower", 
    log = TRUE)
  mito_drop <- isOutlier(
    metric = sce$subsets_Mito_percent, 
    nmads = 3, 
    type = "higher")
  f <- !(libsize_drop | feature_drop | mito_drop)
  return(f)
}



##scpipe 
library(DropletUtils)
library(scPipe)
sc1_p <- "/stornext/HPCScratch/home/you.y/preprocess_update/raw_results/scpipe/sc_5cl/"
sc1 <- create_sce_by_dir(sc1_p, organism = "hsapiens_gene_ensembl", gene_id_type="ensembl_gene_id")
sc1_barcode <- read.csv("/stornext/HPCScratch/home/you.y/preprocess_update/raw_results/scpipe/sc_5cl/barcode_anno.csv",header = TRUE)
colData(sc1) = cbind(colData(sc1), barcode=sc1_barcode$barcode_sequence)

colnames(sc1) <- sc1$barcode
sc_empty <- emptyDrops(counts(sc1), lower = 100)
sc1 <- sc1[,which(sc_empty$FDR<=0.001)]

calculateQCMet(sc1) ->sc1
sc1$SNG.1ST <- barcode_info$SNG.1ST[match(sc1$barcode, barcode_info$barcode)]
sc1$doublet <- barcode_info$doublet[match(sc1$barcode, barcode_info$barcode)]

sc1 <- sc1[rowSums(counts(sc1))>0,]
sc1_f<- scater_filter(sc1)
sc1$keep <-sc1_f

gene_filter(sc1[,sc1$keep])-> sc11

saveRDS(sc11,file.path(write.path,"sc5cl_scpipe.rds"))

#zumis

library(SingleCellExperiment)
zumis_path <- "/stornext/HPCScratch/home/you.y/preprocess_update/raw_results/zumis/sc_5cl/zUMIs_output/expression/"
zumis_sc1_data <- readRDS(file.path(zumis_path,"sc_5cl.dgecounts.rds"))

zumis_sc1 <- SingleCellExperiment(
  assays = list(counts = zumis_sc1_data$umicount$inex$all)
)

zumis_empty <- emptyDrops(counts(zumis_sc1), lower = 100)
zumis_sc1 <- zumis_sc1[,which(zumis_empty$FDR<=0.001)]

zumis_sc1$SNG.1ST <- barcode_info$SNG.1ST[match(colnames(zumis_sc1), barcode_info$barcode)]
zumis_sc1$doublet <- barcode_info$doublet[match(colnames(zumis_sc1), barcode_info$barcode)]

zumis_sc1 <- zumis_sc1[rowSums(counts(zumis_sc1))>0,]
calculateQCMet(zumis_sc1) ->zumis_sc1

zumis_f<- scater_filter(zumis_sc1)

zumis_sc1$keep <-zumis_f
gene_filter(zumis_sc1[,zumis_sc1$keep])-> zumis_sc11

saveRDS(zumis_sc11,file.path(write.path,"sc5cl_zumis.rds"))

#kallisto 

library(Matrix)
kb_path <-"/stornext/General/data/user_managed/grpu_mritchie_1/Yue/preprocess_update/raw_results/kb/sc_5cl/genecount/"
kb_mat <- t(readMM(file=file.path(kb_path,"genes.mtx")))
kb_well <- read.table(file.path(kb_path,"genes.barcodes.txt"))
kb_genes <- read.table(file.path(kb_path,"genes.genes.txt"))
rownames(kb_mat) <- kb_genes[,1]
colnames(kb_mat) <- kb_well[,1]

kb_empty <- emptyDrops(kb_mat, lower = 100)
kb_sc1<- SingleCellExperiment(
  assays = list(counts = as.matrix(kb_mat[,which(kb_empty$FDR<=0.001)]))
)

rownames(kb_sc1) <- unlist(strsplit(as.character(rownames(kb_sc1[])),".",fix=TRUE))[seq(1,nrow(kb_sc1)*2,2)]

calculateQCMet(kb_sc1) ->kb_sc1

kb_sc1$SNG.1ST <- barcode_info$SNG.1ST[match(colnames(kb_sc1), barcode_info$barcode)]
kb_sc1$doublet <- barcode_info$doublet[match(colnames(kb_sc1), barcode_info$barcode)]

kb_sc1 <- kb_sc1[rowSums(counts(kb_sc1))>0,]


kb_f<- scater_filter(kb_sc1)
kb_sc1$keep <-kb_f

gene_filter(kb_sc1[,kb_sc1$keep])-> kb_sc11

saveRDS(kb_sc11,file.path(write.path,"sc5cl_kb.rds"))




#salmon_alevin 

library(fishpond)
library(tximport)
sal_ale_path<-"/stornext/General/data/user_managed/grpu_mritchie_1/Yue/preprocess_update/raw_results/SA/sc_5cl/output/alevin/"

alevin <- tximport(file.path(sal_ale_path,"quants_mat.gz"),type="alevin")
alevin_empty <- emptyDrops(alevin$counts,lower = 100)
alevin_filtered  <- alevin$counts[,which(alevin_empty$FDR<=0.001)]

alevin_sc <- SingleCellExperiment(
  assays = list(counts = as.matrix(alevin_filtered))
)

rownames(alevin_sc) <- unlist(strsplit(as.character(rownames(alevin_sc[])),".",fix=TRUE))[seq(1,nrow(alevin_sc)*2,2)]
barcode_info$SNG.1ST[match(colnames(alevin_sc), barcode_info$barcode)]-> alevin_sc$SNG.1ST
barcode_info$doublet[match(colnames(alevin_sc), barcode_info$barcode)]-> alevin_sc$doublet

calculateQCMet(alevin_sc) ->alevin_sc
alevin_sc <- alevin_sc[rowSums(counts(alevin_sc))>0,]

alevin_f<- scater_filter(alevin_sc)

alevin_sc$keep <-alevin_f
gene_filter(alevin_sc[,alevin_sc$keep])-> alevin_sc1


saveRDS(alevin_sc1,file.path(write.path,"sc5cl_sa.rds"))


#splici
source("/stornext/General/data/user_managed/grpu_mritchie_1/Yue/preprocess_update/raw_results/splici/load_fry.R")
splici_path <-"/stornext/General/data/user_managed/grpu_mritchie_1/Yue/preprocess_update/tools/salmon/quants/sc_5cl/quant/"
af_all = load_fry(splici_path)
af_empty <- emptyDrops(counts(af_all))
#saveRDS(af_empty,file.path(write.path,"sc3cl_splici_empty.rds"))
af  <- af_all[,which(af_empty$FDR<=0.001)]
#af  <- af[,!is.na(af_empty$FDR)]
#af <- af_all[,defaultDrops(counts(af_all))]

barcode_info$SNG.1ST[match(colnames(af), barcode_info$barcode)]-> af$SNG.1ST
barcode_info$doublet[match(colnames(af), barcode_info$barcode)]-> af$doublet

calculateQCMet(af) ->af
af <- af[rowSums(counts(af))>0,]
af_f<- scater_filter(af)

af$keep <-af_f

gene_filter(af[,af$keep])-> af1
saveRDS(af1,file.path(write.path,"sc5cl_splici.rds"))


#dropseqpipe

dropseq_path <- "/stornext/General/data/user_managed/grpu_mritchie_1/Yue/preprocess/dsp/sc_5cl_v2/results/summary/umi/"
drop_mat <- readMM(file.path(dropseq_path,"matrix.mtx"))
drop_well <- read.table(file.path(dropseq_path,"barcodes.tsv"))
drop_genes <- read.table(file.path(dropseq_path,"genes.tsv"))

drop_genes <- mapIds(
  x = EnsDb.Hsapiens.v98, 
  # NOTE: Need to remove gene version number prior to lookup.
  keys = as.character(drop_genes[,1]),
  keytype = "SYMBOL",
  column = "GENEID")
names(drop_genes) <- NULL
rownames(drop_mat) <- drop_genes
colnames(drop_mat) <- as.character(drop_well[,1])

drop_empty <- emptyDrops(drop_mat, lower = 100)

drop_sc1<- SingleCellExperiment(
  assays = list(counts = as.matrix(drop_mat[,which(drop_empty$FDR<=0.001)]))
)
colnames(drop_sc1) <- unlist(strsplit(as.character(colnames(drop_sc1[])),"_",fix=TRUE))[seq(5,ncol(drop_sc1)*5,5)]
barcode_info$SNG.1ST[match(colnames(drop_sc1), barcode_info$barcode)]-> drop_sc1$SNG.1ST
barcode_info$doublet[match(colnames(drop_sc1), barcode_info$barcode)]-> drop_sc1$doublet


calculateQCMet(drop_sc1) ->drop_sc1

drop_sc1 <- drop_sc1[rowSums(counts(drop_sc1))>0,]
drop_f<- scater_filter(drop_sc1)

drop_sc1$keep <-drop_f
gene_filter(drop_sc1[,drop_sc1$keep])-> drop_sc11

saveRDS(drop_sc11,file.path(write.path,"sc5cl_drop.rds"))


#cellranger
cellr_path <-"/stornext/HPCScratch/home/you.y/preprocess_update/raw_results/cellranger/sc_5cl/Lib90_combined/outs/raw_feature_bc_matrix/"
cellrangr_scall <- read10xCounts(cellr_path,col.names=TRUE)
colnames(cellrangr_scall)<-sapply(strsplit(colnames(cellrangr_scall),"-"),"[[",1)
cellr_empty <- emptyDrops(counts(cellrangr_scall), lower = 100)
cellrangr_sc <- cellrangr_scall[,which(cellr_empty$FDR<0.001)]
#cellrangr_sc<- cellrangr_scall[,defaultDrops(counts(cellrangr_scall))]

barcode_info$SNG.1ST[match(colnames(cellrangr_sc), barcode_info$barcode)]-> cellrangr_sc$SNG.1ST
barcode_info$doublet[match(colnames(cellrangr_sc), barcode_info$barcode)]-> cellrangr_sc$doublet


calculateQCMet(cellrangr_sc) ->cellrangr_sc

cellrangr_sc <- cellrangr_sc[rowSums(counts(cellrangr_sc))>0,]
cellrangr_f<- scater_filter(cellrangr_sc)

cellrangr_sc$keep <-cellrangr_f
gene_filter(cellrangr_sc[,cellrangr_sc$keep])-> cellrangr_sc1


saveRDS(cellrangr_sc1,file.path(write.path,"sc5cl_cellranger.rds"))





###all together

bind_rows(data.frame(data="scpipe",total=sc11$total,detect=sc11$detected, label=sc11$SNG.1ST), 
          data.frame(data="zumis",total=zumis_sc11$total,detect=zumis_sc11$detected,label=zumis_sc11$SNG.1ST),
          data.frame(data="drop",total=drop_sc11$total,detect=drop_sc11$detected,label=drop_sc11$SNG.1ST),
          data.frame(data="cellranger",total=cellrangr_sc1$total,detect=cellrangr_sc1$detected,label=cellrangr_sc1$SNG.1ST),
          data.frame(data="splici",total=af1$total,detect=af1$detected,label=af1$SNG.1ST),
          #data.frame(data="optimus",total=optimus_sce1$total,detect=optimus_sce1$detected,label=optimus_sce1$SNG.1ST),
          data.frame(data="sa",total=alevin_sc1$total,detect=alevin_sc1$detected,label=alevin_sc1$SNG.1ST),
          data.frame(data="kb",total=kb_sc11$total,detect=kb_sc11$detected,label=kb_sc11$SNG.1ST)) -> tmp

tmp$data <- recode(tmp$data , "scpipe"="scPipe","zumis"="zUMIs","kb"="kallisto bustools",
                   "optimus" ="Optimus","drop"="dropSeqPipe","cellranger"="Cell Ranger",
                   "sa"="salmon_SA","splici"="salmon_splici")

tmp$log10total <- log10(tmp$total)
tmp$log10detect <- log10(tmp$detect)
ggplot(tmp,aes(x=data,y=log10total,col=data)) +geom_violin() 
tmp$labelled <- !is.na(tmp$label)
col=c("lightblue","darkblue")
names(col) <- c(TRUE,FALSE)
library(ggpubr)
ggplot(tmp) +
  geom_violin(aes(x=data,y=log10detect,col=data)) +
  #geom_jitter(aes(x=data,y=log10detect,col=labelled),alpha=0.1)+
  scale_color_manual(values=c(droplet_col,col))+
  theme_bw()
ggplot(tmp) +
  geom_violin(aes(x=data,y=log10total,col=data)) +
  #geom_jitter(aes(x=data,y=log10detect,col=labelled),alpha=0.1)+
  scale_color_manual(values=c(droplet_col,col))+
  theme_bw()
ggviolin(tmp, x = "data", y = "log10detect", fill = "white",col="data",
         palette = c(droplet_col,col), legend="right" ,add = "jitter", add.params = list(color = "labelled",alpha=0.2))






