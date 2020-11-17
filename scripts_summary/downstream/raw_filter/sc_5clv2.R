library(scran)
library(scater)
library(CellBench)
library(R.utils)
library(tidyverse)

setwd("/stornext/General/data/user_managed/grpu_mritchie_1/Yue/preprocess/new")
data.path <- "/stornext/General/data/user_managed/grpu_mritchie_1/Yue/preprocess/data"
write.path <- "/stornext/General/data/user_managed/grpu_mritchie_1/Yue/preprocess/new/SCE/raw"

barcode_info <- read.table("/stornext/General/data/user_managed/grpu_mritchie_1/SCmixology/scbench_5cellline_10x/scPipe_output/Lib90.best",header = TRUE)

sapply(strsplit(as.character(barcode_info$BEST),"-"),function(x) x[1]) -> barcode_info$doublet
sapply(strsplit(as.character(barcode_info$BARCODE),"-1"),"[",1)->barcode_info$BARCODE
empty.path <- "/stornext/General/data/user_managed/grpu_mritchie_1/Yue/preprocess/SCE/empty"

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

plotpath <- "/stornext/General/data/user_managed/grpu_mritchie_1/Yue/preprocess/new/results/raw/"
library(edgeR)

##scpipe 

library(scPipe)
sc1_p <- "/stornext/General/data/user_managed/grpu_mritchie_1/Yue/preprocess/scpipe/sc_5clv2_2/"
sc1 <- create_sce_by_dir(sc1_p, organism = "hsapiens_gene_ensembl", gene_id_type="ensembl_gene_id")
sc1_barcode <- read.csv("/stornext/General/data/user_managed/grpu_mritchie_1/Yue/preprocess/scpipe/sc_5clv2_2/barcode_anno.csv",header = TRUE)
colData(sc1) = cbind(colData(sc1), barcode=sc1_barcode$barcode_sequence)

colnames(sc1) <- sc1$barcode

sc_empty <- readRDS(file.path(empty.path,"scpipe_5cl_emptied.rds"))
sc1 <- sc1[,which(sc_empty$FDR<=0.001)]
calculateQCMet(sc1) ->sc1
sc1$SNG.1ST <- barcode_info$SNG.1ST[match(sc1$barcode, barcode_info$BARCODE)]
sc1$demuxlet_cls <- barcode_info$doublet[match(sc1$barcode, barcode_info$BARCODE)]

#zero
sc1 <- sc1[rowSums(counts(sc1))>0,]
data.frame(total_counts_per_gene =log10(rowSums(counts(sc1))) ,
           pct.zeros= rowSums(counts(sc1)==0)/ncol(sc1),
           preprocess="scPipe") -> zero.d.scpipe

#biotype
library(tidyverse)
as.data.frame(rowData(sc1)) %>% dplyr::group_by(gene_biotype) %>% summarise(count=n()) %>% mutate(preprocess="scPipe") ->gene_biotype_scpipe


#filter by remove outlier
sc2_p <- "/stornext/General/data/user_managed/grpu_mritchie_1/Yue/preprocess/scpipe/sc_5cl_v2/"
sc2 <- create_sce_by_dir(sc2_p, organism = "hsapiens_gene_ensembl", gene_id_type="ensembl_gene_id")
sc2_barcode <- read.csv("/stornext/General/data/user_managed/grpu_mritchie_1/Yue/preprocess/scpipe/sc_5cl_v2/barcode_anno.csv",header = TRUE)
colData(sc2) = cbind(colData(sc2), barcode=sc2_barcode$barcode_sequence)

sc2 <- calculate_QC_metrics(sc2)
calculateQCMet(sc2) ->sc2
sc2 <- sc2[,colSums(counts(sc2))>0]
scp_sc1_qc <- detect_outlier(sc2, type= "low", comp=2)
table(QC_metrics(scp_sc1_qc)$outliers)


scp_sc1_qc$keep <- FALSE
scp_sc1_qc$keep[scp_sc1_qc$outliers=="FALSE"] <- TRUE

scp_sc1_qc$SNG.1ST <- barcode_info$SNG.1ST[match(scp_sc1_qc$barcode, barcode_info$BARCODE)]
scp_sc1_qc$demuxlet_cls <- barcode_info$doublet[match(scp_sc1_qc$barcode, barcode_info$BARCODE)]

#plotfc(scp_sc1_qc,"sc_5clv3/fc/scpipe_rmout_fc")

data.frame(total_counts_per_cell = log10(scp_sc1_qc$total),
           kept = scp_sc1_qc$keep,
           Mito_percent_per_cell= scp_sc1_qc$subsets_Mito_percent,
           detected = scp_sc1_qc$detected,
           preprocess="scPipe_rm") -> filter.scpipe.rm



#use scater
scater_filter(sc1) -> sc1_f
table(sc1_f)
sc1$keep <-sc1_f


#plotfc(sc1,"sc_5clv3/fc/scpipe_scater_fc")


data.frame(total_counts_per_cell = log10(sc1$total),
           kept = sc1$keep,
           Mito_percent_per_cell= sc1$subsets_Mito_percent,
           detected = sc1$detected,
           preprocess="scPipe") -> filter.scpipe

saveRDS(gene_filter(sc1[,sc1$keep]),file.path(write.path,"sc5clv2_scpipe.rds"))

#zumis

library(SingleCellExperiment)
library(DropletUtils)
zumis_path <- "/stornext/General/data/user_managed/grpu_mritchie_1/Yue/preprocess/zumis2/sc_5cl/zUMIs_output/expression"
zumis_sc1_data <- readRDS(file.path(zumis_path,"sc_5cl.dgecounts.rds"))

zumis_sc1 <- SingleCellExperiment(
  assays = list(counts = zumis_sc1_data$umicount$inex$all)
)


zumis_empty <- readRDS(file.path(empty.path,"zumis_5cl_emptied.rds"))
#defaultDrops(zumis_sc1_data$umicount$inex$all) -> c.keep
zumis_sc1 <- zumis_sc1[,which(zumis_empty$FDR<=0.001)]

zumis_sc1$SNG.1ST <- barcode_info$SNG.1ST[match(colnames(zumis_sc1), barcode_info$BARCODE)]
zumis_sc1$demuxlet_cls <- barcode_info$doublet[match(colnames(zumis_sc1), barcode_info$BARCODE)]

zumis_sc1 <- zumis_sc1[rowSums(counts(zumis_sc1))>0,]
data.frame(total_counts_per_gene =log10(rowSums(counts(zumis_sc1))) ,
           pct.zeros= rowSums(counts(zumis_sc1)==0)/ncol(zumis_sc1),
           preprocess="zUMIs") -> zero.d.zumis
calculateQCMet(zumis_sc1) ->zumis_sc1


zumis_f<- scater_filter(zumis_sc1)

as.data.frame(rowData(zumis_sc1)) %>% dplyr::group_by(gene_biotype) %>% summarise(count=n()) %>% mutate(preprocess="zUMIs") ->gene_biotype_zumis


zumis_sc1$keep <-zumis_f

data.frame(total_counts_per_cell = log10(zumis_sc1$total),
           kept = zumis_sc1$keep,
           Mito_percent_per_cell= zumis_sc1$subsets_Mito_percent,
           detected = zumis_sc1$detected,
           preprocess="zUMIs") -> filter.zumis


#plotfc(zumis_sc1,"sc_5clv3/fc/zumis_scater_fc")


#automatic
zumis_path2 <- "/stornext/General/data/user_managed/grpu_mritchie_1/Yue/preprocess/zumis/sc_5cl/zUMIs_output/expression"
zumis_sc2_data <- readRDS(file.path(zumis_path2,"scbench_5cellline_10x.dgecounts.rds"))

zumis_sc2 <- SingleCellExperiment(
  assays = list(counts = zumis_sc2_data$umicount$inex$all)
)

zumis_sc2 <- zumis_sc2[rowSums(counts(zumis_sc2))>0,]

calculateQCMet(zumis_sc2) ->zumis_sc2

zumis_f2<- scater_filter(zumis_sc2)
zumis_sc2$SNG.1ST <- barcode_info$SNG.1ST[match(colnames(zumis_sc2), barcode_info$BARCODE)]
zumis_sc2$demuxlet_cls <- barcode_info$doublet[match(colnames(zumis_sc2), barcode_info$BARCODE)]

zumis_sc2$keep <-zumis_f2
#plotfc(zumis_sc1,"sc_5clv3/fc/zumis_auto_scater_fc")


data.frame(total_counts_per_cell = log10(zumis_sc2$total),
           kept = zumis_sc2$keep,
           Mito_percent_per_cell= zumis_sc2$subsets_Mito_percent,
           detected = zumis_sc2$detected,
           preprocess="zUMIs_auto") -> filter.zumis.auto

saveRDS(gene_filter(zumis_sc1[,zumis_sc1$keep]),file.path(write.path,"sc5clv2_zumis.rds"))



#kallisto 

library(Matrix)

kb_path <-"/stornext/General/data/user_managed/grpu_mritchie_1/Yue/preprocess/kallisto_bus/sc_5cl_v2/genecount/"
kb_mat <- t(readMM(file=file.path(kb_path,"genes.mtx")))
kb_well <- read.table(file.path(kb_path,"genes.barcodes.txt"))
kb_genes <- read.table(file.path(kb_path,"genes.genes.txt"))
rownames(kb_mat) <- kb_genes[,1]
colnames(kb_mat) <- kb_well[,1]


kb_empty <- readRDS(file.path(empty.path,"kb_5cl_v2emptied.rds"))

kb_sc1<- SingleCellExperiment(
  assays = list(counts = as.matrix(kb_mat[,which(kb_empty$FDR<=0.001)]))
)

rownames(kb_sc1) <- unlist(strsplit(as.character(rownames(kb_sc1[])),".",fix=TRUE))[seq(1,nrow(kb_sc1)*2,2)]

calculateQCMet(kb_sc1) ->kb_sc1

kb_sc1$SNG.1ST <- barcode_info$SNG.1ST[match(colnames(kb_sc1), barcode_info$BARCODE)]
kb_sc1$demuxlet_cls <- barcode_info$doublet[match(colnames(kb_sc1), barcode_info$BARCODE)]


kb_sc1 <- kb_sc1[rowSums(counts(kb_sc1))>0,]

data.frame(total_counts_per_gene =log10(rowSums(counts(kb_sc1))) ,
           pct.zeros= rowSums(counts(kb_sc1)==0)/ncol(kb_sc1),
           preprocess="kallisto bustools") -> zero.d.kb


as.data.frame(rowData(kb_sc1)) %>% dplyr::group_by(gene_biotype) %>% summarise(count=n()) %>% mutate(preprocess="kallisto bustools") ->gene_biotype_kb

kb_f<- scater_filter(kb_sc1)
kb_sc1$keep <-kb_f
kb_sc1 <- kb_sc1[rowSums(counts(kb_sc1))>0,]
#plotfc(kb_sc1,"sc_5clv3/fc/kb_scater_fc")

data.frame(total_counts_per_cell = log10(kb_sc1$total),
           kept = kb_sc1$keep,
           Mito_percent_per_cell= kb_sc1$subsets_Mito_percent,
           detected = kb_sc1$detected,
           preprocess="kallisto bustools") -> filter.kb

saveRDS(gene_filter(kb_sc1[,kb_sc1$keep]),file.path(write.path,"sc5clv2_kb.rds"))




#salmon_alevin 

library(fishpond)
alevin_empty <- readRDS(file.path(empty.path,"alevin_5cl_v2.rds"))
library(tximport)
sal_ale_path<-"/stornext/General/data/user_managed/grpu_mritchie_1/Yue/preprocess/sal_ale/sc_5cl_v2/alevin_output/alevin/"
alevin <- tximport(file.path(sal_ale_path,"quants_mat.gz"),type="alevin")
alevin_filtered  <- alevin$counts[,which(alevin_empty$FDR<=0.001)]

alevin_sc <- SingleCellExperiment(
  assays = list(counts = as.matrix(alevin_filtered))
)
rownames(alevin_sc) <- unlist(strsplit(as.character(rownames(alevin_sc[])),".",fix=TRUE))[seq(1,nrow(alevin_sc)*2,2)]
barcode_info$SNG.1ST[match(colnames(alevin_sc), barcode_info$BARCODE)]-> alevin_sc$SNG.1ST
alevin_sc$demuxlet_cls <- barcode_info$doublet[match(colnames(alevin_sc), barcode_info$BARCODE)]


calculateQCMet(alevin_sc) ->alevin_sc
alevin_sc <- alevin_sc[rowSums(counts(alevin_sc))>0,]
#alevin_sc <- alevin_sc[,alevin_sc$total >100,]
data.frame(total_counts_per_gene =log10(rowSums(counts(alevin_sc))) ,
           pct.zeros= rowSums(counts(alevin_sc)==0)/ncol(alevin_sc),
           preprocess="salmon alevin") -> zero.d.alevin


as.data.frame(rowData(alevin_sc)) %>% dplyr::group_by(gene_biotype) %>% summarise(count=n()) %>% mutate(preprocess="salmon alevin") ->gene_biotype_alevin


alevin_f<- scater_filter(alevin_sc)
#alevin_sc <- alevin_sc[rowSums(counts(alevin_sc))>0,]

alevin_sc$keep <-alevin_f
#plotfc(alevin_sc,"sc_5clv3/fc/alevin_scater_fc")


data.frame(total_counts_per_cell = log10(alevin_sc$total),
           kept = alevin_sc$keep,
           Mito_percent_per_cell= alevin_sc$subsets_Mito_percent,
           detected = alevin_sc$detected,
           preprocess="salmon alevin") -> filter.alevin

saveRDS(gene_filter(alevin_sc[,alevin_sc$keep]),file.path(write.path,"sc5clv2_alevin.rds"))

#dropseqpipe
#######
dropseq_path <- "/stornext/General/data/user_managed/grpu_mritchie_1/Yue/preprocess/dsp/sc_5cl_v2/results/summary/umi"
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
drop_empty <- readRDS(file.path(empty.path,"dropseqpipe_5cl_emptied.rds"))

drop_mat <- drop_mat[,which(drop_empty$FDR<=0.001)]
drop_sc1<- SingleCellExperiment(
  assays = list(counts = as.matrix(drop_mat))
)
colnames(drop_sc1) <- unlist(strsplit(as.character(colnames(drop_sc1[])),"_",fix=TRUE))[seq(5,ncol(drop_sc1)*5,5)]
barcode_info$SNG.1ST[match(colnames(drop_sc1), barcode_info$BARCODE)]-> drop_sc1$SNG.1ST
drop_sc1$demuxlet_cls <- barcode_info$doublet[match(colnames(drop_sc1), barcode_info$BARCODE)]


calculateQCMet(drop_sc1) ->drop_sc1

data.frame(total_counts_per_gene =log10(rowSums(counts(drop_sc1))) ,
           pct.zeros= rowSums(counts(drop_sc1)==0)/ncol(drop_sc1),
           preprocess="dropSeqPipe") -> zero.d.drop

drop_sc1 <- drop_sc1[rowSums(counts(drop_sc1))>0,]
as.data.frame(rowData(drop_sc1)) %>% dplyr::group_by(gene_biotype) %>% summarise(count=n()) %>% mutate(preprocess="dropSeqPipe") ->gene_biotype_drop


drop_f<- scater_filter(drop_sc1)

drop_sc1$keep <-drop_f
#plotfc(drop_sc1,"sc_5clv3/fc/drop_scater_fc")


data.frame(total_counts_per_cell = log10(drop_sc1$total),
           kept = drop_sc1$keep,
           Mito_percent_per_cell= drop_sc1$subsets_Mito_percent,
           detected = drop_sc1$detected,
           preprocess="dropSeqPipe") -> filter.drop

saveRDS(gene_filter(drop_sc1[,drop_sc1$keep]),file.path(write.path,"sc5clv2_drop.rds"))



#cellranger
cellranger_empty <- readRDS(file.path(empty.path,"cellranger_5cl_v2_emptied.rds"))
cellr_path <-"/stornext/General/data/user_managed/grpu_mritchie_1/Yue/preprocess/cellranger/sc_5cl/Lib90_combined/outs/raw_feature_bc_matrix/"
cellrangr_sc <- read10xCounts(cellr_path,col.names=TRUE)
colnames(cellrangr_sc)<-sapply(strsplit(colnames(cellrangr_sc),"-"),"[[",1)
cellrangr_sc <- cellrangr_sc[,which(cellranger_empty$FDR<=0.001)]

barcode_info$SNG.1ST[match(colnames(cellrangr_sc), barcode_info$BARCODE)]-> cellrangr_sc$SNG.1ST
cellrangr_sc$demuxlet_cls <- barcode_info$doublet[match(colnames(cellrangr_sc), barcode_info$BARCODE)]


calculateQCMet(cellrangr_sc) ->cellrangr_sc

data.frame(total_counts_per_gene =log10(rowSums(counts(cellrangr_sc))) ,
           pct.zeros= rowSums(counts(cellrangr_sc)==0)/ncol(cellrangr_sc),
           preprocess="Cell Ranger") -> zero.d.cellranger


cellrangr_sc <- cellrangr_sc[rowSums(counts(cellrangr_sc))>0,]
as.data.frame(rowData(cellrangr_sc)) %>% dplyr::group_by(gene_biotype) %>% summarise(count=n()) %>% mutate(preprocess="Cell Ranger") ->gene_biotype_cellranger


cellrangr_f<- scater_filter(cellrangr_sc)

cellrangr_sc$keep <-cellrangr_f
#plotfc(cellrangr_sc,"sc_5clv3/fc/cellranger_scater_fc")

data.frame(total_counts_per_cell = log10(cellrangr_sc$total),
           kept = cellrangr_sc$keep,
           Mito_percent_per_cell= cellrangr_sc$subsets_Mito_percent,
           detected = cellrangr_sc$detected,
           preprocess="Cell Ranger") -> filter.cellranger

saveRDS(gene_filter(cellrangr_sc[,cellrangr_sc$keep]),file.path(write.path,"sc5clv2_cellranger.rds"))

#optimus
optimus_path <- "/stornext/General/data/user_managed/grpu_mritchie_1/Yue/preprocess/optimus/sc5clv2/"
readRDS(file.path(optimus_path,"optimus.rds"))->optimus_sc

barcode_info$SNG.1ST[match(colnames(optimus_sc), barcode_info$BARCODE)]-> optimus_sc$SNG.1ST
optimus_sc$demuxlet_cls <- barcode_info$doublet[match(colnames(optimus_sc), barcode_info$BARCODE)]

calculateQCMet(optimus_sc) ->optimus_sc

data.frame(total_counts_per_gene =log10(rowSums(counts(optimus_sc))) ,
           pct.zeros= rowSums(counts(optimus_sc)==0)/ncol(optimus_sc),
           preprocess="Optimus") -> zero.d.optimus


optimus_sc <- optimus_sc[rowSums(counts(optimus_sc))>0,]
as.data.frame(rowData(optimus_sc)) %>% dplyr::group_by(gene_biotype) %>% summarise(count=n()) %>% mutate(preprocess="Optimus") ->gene_biotype_optimus


optimus_f<- scater_filter(optimus_sc)

optimus_sc$keep <-optimus_f
#plotfc(optimus_sc,"sc_5clv3/fc/optimus_scater_fc")

data.frame(total_counts_per_cell = log10(optimus_sc$total),
           kept = optimus_sc$keep,
           Mito_percent_per_cell= optimus_sc$subsets_Mito_percent,
           detected = optimus_sc$detected,
           preprocess="Optimus") -> filter.optimus

saveRDS(gene_filter(optimus_sc[,optimus_sc$keep]),file.path(write.path,"sc5clv2_optimus.rds"))





