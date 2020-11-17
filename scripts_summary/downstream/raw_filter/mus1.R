library(scran)
library(scater)
library(CellBench)
library(R.utils)
library(tidyverse)

setwd("/stornext/General/data/user_managed/grpu_mritchie_1/Yue/preprocess/new")
data.path <- "/stornext/General/data/user_managed/grpu_mritchie_1/Yue/preprocess/data"
write.path <- "/stornext/General/data/user_managed/grpu_mritchie_1/Yue/preprocess/new/SCE/raw"

barcode_info <- read.csv("/stornext/General/data/user_managed/grpu_mritchie_1/Yue/preprocess/mouse/annotations_droplet.csv",header = TRUE)
barcode_info[which(barcode_info$channel=="10X_P7_8"),] -> barcode_info
barcode_info$barcode <- gsub("^.*_", "", unlist(barcode_info$cell))
barcode_info %>% dplyr::select(barcode, cell_ontology_class) -> barcode_info
barcode_info$SNG.1ST <- barcode_info$cell_ontology_class

empty.path <- "/stornext/General/data/user_managed/grpu_mritchie_1/Yue/preprocess/SCE/empty"

library(AnnotationHub)
ah <- AnnotationHub()
EnsDb.Mus.v99 <- query(ah, c("EnsDb", "Mus musculus", 99))[[1]]

gene_filter <- function(sce){
  keep1 = (apply(counts(sce), 1, function(x) mean(x[x>0])) > 1)  # average count larger than 1
  keep2 = (rowSums(counts(sce)>0) > 10)  
  table(keep1&keep2)
  sce = sce[(keep1 & keep2), ]
  return(sce)
}


calculateQCMet <- function(sce){
  location <- mapIds(
    x = EnsDb.Mus.v99, 
    # NOTE: Need to remove gene version number prior to lookup.
    keys = rownames(sce),
    keytype = "GENEID",
    column = "SEQNAME")
  rowData(sce)$CHR =location
  
  symbol <- mapIds(
    x = EnsDb.Mus.v99, 
    # NOTE: Need to remove gene version number prior to lookup.
    keys = rownames(sce),
    keytype = "GENEID",
    column = "SYMBOL")
  rowData(sce)$symbol =symbol
  
  biotype<- mapIds(
    x = EnsDb.Mus.v99, 
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
plotfc <- function(sce, name){
  keep <- sce$keep
  lost <- calculateAverage(counts(sce)[, !keep])
  kept <- calculateAverage(counts(sce)[, keep])
  logged <- cpm(cbind(lost, kept), log = TRUE, prior.count = 2)
  
  logFC <- logged[, 1] - logged[, 2]
  logfc_symbol <- rowData(sce)$symbol
  abundance <- rowMeans(logged)
  
  mito_set <- rowData(sce)$symbol[which(rowData(sce)$CHR == "MT")]
  ribo_set <- rowData(sce)$symbol[grep("^Rp(s|l)", rowData(sce)$symbol)]
  
  is_mito <- rowData(sce)$symbol %in% mito_set
  is_ribo <- rowData(sce)$symbol %in% ribo_set
  png(paste0(plotpath,name,".png"))
  
  plot(
    abundance,
    logFC,
    xlab = "Average count",
    ylab = "Log-FC (lost/kept)",
    pch = 16)
  points(
    abundance[is_mito],
    logFC[is_mito],
    col = "dodgerblue",
    pch = 16,
    cex = 1)
  points(
    abundance[is_ribo],
    logFC[is_ribo],
    col = "orange",
    pch = 16,
    cex = 1)
  abline(h = c(-1, 1), col = "red", lty = 2)
  
  dev.off()
  
  text = data.frame(gene_id=names(logFC[(c(order(logFC)[1:20],order(logFC)[(length(logFC)-20):length(logFC)]))]),
                    symbol = logfc_symbol[(c(order(logFC)[1:20],order(logFC)[(length(logFC)-20):length(logFC)]))]
  )
  
  write.table(text,paste0(plotpath,name,".txt"),quote = FALSE,row.names = FALSE,sep="\t")
  
}



plotfclabel <- function(sce, name){
  sce <- sce[,sce$keep]
  keep <- sce$no_label
  lost <- calculateAverage(counts(sce)[, !keep])
  kept <- calculateAverage(counts(sce)[, keep])
  logged <- cpm(cbind(lost, kept), log = TRUE, prior.count = 2)
  
  logFC <- logged[, 1] - logged[, 2]
  logfc_symbol <- rowData(sce)$symbol
  abundance <- rowMeans(logged)
  
  mito_set <- rowData(sce)$symbol[which(rowData(sce)$CHR == "MT")]
  ribo_set <- rowData(sce)$symbol[grep("^Rp(s|l)", rowData(sce)$symbol)]
  
  is_mito <- rowData(sce)$symbol %in% mito_set
  is_ribo <- rowData(sce)$symbol %in% ribo_set
  png(paste0(plotpath,name,"_label.png"))
  
  plot(
    abundance,
    logFC,
    xlab = "Average count",
    ylab = "Log-FC (lost/kept)",
    pch = 16)
  points(
    abundance[is_mito],
    logFC[is_mito],
    col = "dodgerblue",
    pch = 16,
    cex = 1)
  points(
    abundance[is_ribo],
    logFC[is_ribo],
    col = "orange",
    pch = 16,
    cex = 1)
  abline(h = c(-1, 1), col = "red", lty = 2)
  
  dev.off()
  
  text = data.frame(gene_id=names(logFC[(c(order(logFC)[1:20],order(logFC)[(length(logFC)-20):length(logFC)]))]),
                    symbol = logfc_symbol[(c(order(logFC)[1:20],order(logFC)[(length(logFC)-20):length(logFC)]))]
  )
  
  write.table(text,paste0(plotpath,name,"_label.txt"),quote = FALSE,row.names = FALSE,sep="\t")
  
}


##scpipe 

library(scPipe)
sc1_p <- "/stornext/General/data/user_managed/grpu_mritchie_1/Yue/preprocess/scpipe/mus1_2/"
sc1 <- create_sce_by_dir(sc1_p, organism = "mmusculus_gene_ensembl", gene_id_type="ensembl_gene_id")
sc1_barcode <- read.csv("/stornext/General/data/user_managed/grpu_mritchie_1/Yue/preprocess/scpipe/mus1_2/barcode_anno.csv",header = TRUE)
colData(sc1) = cbind(colData(sc1), barcode=sc1_barcode$barcode_sequence)

colnames(sc1) <- sc1$barcode

sc_empty <- readRDS(file.path(empty.path,"scpipe_mus1_emptied.rds"))
sc1 <- sc1[,which(sc_empty$FDR<=0.001)]
calculateQCMet(sc1) ->sc1
sc1$SNG.1ST <- barcode_info$cell_ontology_class[match(sc1$barcode, barcode_info$barcode)]
#zero
sc1 <- sc1[rowSums(counts(sc1))>0,]
data.frame(total_counts_per_gene =log10(rowSums(counts(sc1))) ,
           pct.zeros= rowSums(counts(sc1)==0)/ncol(sc1),
           preprocess="scPipe") -> zero.d.scpipe

#biotype
library(tidyverse)
as.data.frame(rowData(sc1)) %>% dplyr::group_by(gene_biotype) %>% summarise(count=n()) %>% mutate(preprocess="scPipe") ->gene_biotype_scpipe


#filter by remove outlier
sc2_p <- "/stornext/General/data/user_managed/grpu_mritchie_1/Yue/preprocess/scpipe/mus1/"
sc2 <- create_sce_by_dir(sc2_p, organism = "mmusculus_gene_ensembl", gene_id_type="ensembl_gene_id")
sc2_barcode <- read.csv("/stornext/General/data/user_managed/grpu_mritchie_1/Yue/preprocess/scpipe/mus1/barcode_anno.csv",header = TRUE)
colData(sc2) = cbind(colData(sc2), barcode=sc2_barcode$barcode_sequence)

sc2 <- calculate_QC_metrics(sc2)
calculateQCMet(sc2) ->sc2
scp_sc1_qc <- detect_outlier(sc2, type= "low", comp=2)
table(QC_metrics(scp_sc1_qc)$outliers)


scp_sc1_qc$keep <- FALSE
scp_sc1_qc$keep[scp_sc1_qc$outliers=="FALSE"] <- TRUE
scp_sc1_qc$SNG.1ST <- barcode_info$cell_ontology_class[match(scp_sc1_qc$barcode, barcode_info$barcode)]
scp_sc1_qc$no_label <- !is.na(scp_sc1_qc$SNG.1ST)
#plotfc(scp_sc1_qc,"mus1/fc/scpipe_rmout_fc")
data.frame(total_counts_per_cell = log10(scp_sc1_qc$total),
           kept = scp_sc1_qc$keep,
           label=scp_sc1_qc$no_label,
           Mito_percent_per_cell= scp_sc1_qc$subsets_Mito_percent,
           detected = scp_sc1_qc$detected,
           preprocess="scPipe_rm") -> filter.scpipe.rm



#use scater
scater_filter(sc1) -> sc1_f
table(sc1_f)
sc1$keep <-sc1_f



#all cells kept by scater
#plotfc(sc1,"mus1/fc/scpipe_scater_fc")
sc1$no_label <- !is.na(sc1$SNG.1ST)
#plotfclabel(sc1,"mus1/fc/scpipe_fc_label")

data.frame(total_counts_per_cell = log10(sc1$total),
           kept = sc1$keep,
           label=sc1$no_label,
           Mito_percent_per_cell= sc1$subsets_Mito_percent,
           detected = sc1$detected,
           preprocess="scPipe") -> filter.scpipe

saveRDS(gene_filter(sc1[,sc1$keep]),file.path(write.path,"mus1_scpipe.rds"))

#zumis

library(SingleCellExperiment)
zumis_path <- "/stornext/General/data/user_managed/grpu_mritchie_1/Yue/preprocess/zumis2/mus1/zUMIs_output/expression"
zumis_sc1_data <- readRDS(file.path(zumis_path,"mouse_lung1.dgecounts.rds"))

zumis_sc1 <- SingleCellExperiment(
  assays = list(counts = zumis_sc1_data$umicount$inex$all)
)

zumis_empty <- readRDS(file.path(empty.path,"zumis_mus1_emptied.rds"))
zumis_sc1 <- zumis_sc1[,which(zumis_empty$FDR<=0.001)]
zumis_sc1$SNG.1ST <- barcode_info$cell_ontology_class[match(colnames(zumis_sc1), barcode_info$barcode)]

zumis_sc1 <- zumis_sc1[rowSums(counts(zumis_sc1))>0,]
data.frame(total_counts_per_gene =log10(rowSums(counts(zumis_sc1))) ,
           pct.zeros= rowSums(counts(zumis_sc1)==0)/ncol(zumis_sc1),
           preprocess="zUMIs") -> zero.d.zumis
calculateQCMet(zumis_sc1) ->zumis_sc1


zumis_f<- scater_filter(zumis_sc1)

as.data.frame(rowData(zumis_sc1)) %>% dplyr::group_by(gene_biotype) %>% summarise(count=n()) %>% mutate(preprocess="zUMIs") ->gene_biotype_zumis


zumis_sc1$keep <-zumis_f
zumis_sc1$no_label <- !is.na(zumis_sc1$SNG.1ST)
data.frame(total_counts_per_cell = log10(zumis_sc1$total),
           kept = zumis_sc1$keep,
           label=zumis_sc1$no_label,
           Mito_percent_per_cell= zumis_sc1$subsets_Mito_percent,
           detected = zumis_sc1$detected,
           preprocess="zUMIs") -> filter.zumis

#plotfc(zumis_sc1,"mus1/fc/zumis_scater_fc")
#plotfclabel(zumis_sc1,"mus1/fc/zumis_fc_label")

#automatic
zumis_path2 <- "/stornext/General/data/user_managed/grpu_mritchie_1/Yue/preprocess/zumis/mouse/mus1/zUMIs_output/expression"
zumis_sc2_data <- readRDS(file.path(zumis_path2,"mouse_lung1.dgecounts.rds"))

zumis_sc2 <- SingleCellExperiment(
  assays = list(counts = zumis_sc2_data$umicount$inex$all)
)

zumis_sc2 <- zumis_sc2[rowSums(counts(zumis_sc2))>0,]

calculateQCMet(zumis_sc2) ->zumis_sc2

zumis_f2<- scater_filter(zumis_sc2)
zumis_sc2$SNG.1ST <- barcode_info$cell_ontology_class[match(colnames(zumis_sc2), barcode_info$barcode)]
zumis_sc2$keep <-zumis_f2
#plotfc(zumis_sc1,"mus1/fc/zumis_auto_scater_fc")

zumis_sc2$no_label <- !is.na(zumis_sc2$SNG.1ST)

data.frame(total_counts_per_cell = log10(zumis_sc2$total),
           kept = zumis_sc2$keep,
           label= zumis_sc2$no_label,
           Mito_percent_per_cell= zumis_sc2$subsets_Mito_percent,
           detected = zumis_sc2$detected,
           preprocess="zUMIs_auto") -> filter.zumis.auto

saveRDS(gene_filter(zumis_sc1[,zumis_sc1$keep]),file.path(write.path,"mus1_zumis.rds"))

#exon

zumis_sc3 <- SingleCellExperiment(
  assays = list(counts = zumis_sc1_data$umicount$exon$all)
)

zumis_empty2 <- readRDS(file.path(empty.path,"zumis_mus1_emptied2.rds"))
zumis_sc3 <- zumis_sc3[,which(zumis_empty2$FDR<=0.001)]

zumis_sc3 <- zumis_sc3[rowSums(counts(zumis_sc3))>0,]

calculateQCMet(zumis_sc3) ->zumis_sc3


as.data.frame(rowData(zumis_sc3)) %>% dplyr::group_by(gene_biotype) %>% summarise(count=n()) %>% mutate(preprocess="zUMIs_exon") ->gene_biotype_zumis_exon


#kallisto 

library(Matrix)

kb_empty <- readRDS(file.path(empty.path,"kb_mus1_emptied.rds"))
kb_path <-"/stornext/General/data/user_managed/grpu_mritchie_1/Yue/preprocess/kallisto_bus/mus1/genecount/"
kb_mat <- t(readMM(file=file.path(kb_path,"genes.mtx")))
kb_well <- read.table(file.path(kb_path,"genes.barcodes.txt"))
kb_genes <- read.table(file.path(kb_path,"genes.genes.txt"))
rownames(kb_mat) <- kb_genes[,1]
colnames(kb_mat) <- kb_well[,1]
kb_mat <- kb_mat[,which(kb_empty$FDR<=0.001)]

kb_sc1 <- SingleCellExperiment(
  assays = list(counts = as.matrix(kb_mat))
)

colnames(kb_sc1) <- colnames(kb_mat)
rownames(kb_sc1) <- unlist(strsplit(as.character(rownames(kb_sc1[])),".",fix=TRUE))[seq(1,nrow(kb_sc1)*2,2)]

calculateQCMet(kb_sc1) ->kb_sc1

kb_sc1$SNG.1ST <- barcode_info$cell_ontology_class[match(colnames(kb_sc1), barcode_info$barcode)]

data.frame(total_counts_per_gene =log10(rowSums(counts(kb_sc1))) ,
           pct.zeros= rowSums(counts(kb_sc1)==0)/ncol(kb_sc1),
           preprocess="kallisto bustools") -> zero.d.kb

kb_sc1 <- kb_sc1[rowSums(counts(kb_sc1))>0,]

as.data.frame(rowData(kb_sc1)) %>% dplyr::group_by(gene_biotype) %>% summarise(count=n()) %>% mutate(preprocess="kallisto bustools") ->gene_biotype_kb

kb_f<- scater_filter(kb_sc1)
kb_sc1$keep <-kb_f
kb_sc1 <- kb_sc1[rowSums(counts(kb_sc1))>0,]
#plotfc(kb_sc1,"mus1/fc/kb_scater_fc")

kb_sc1$no_label <- !is.na(kb_sc1$SNG.1ST)
#plotfclabel(kb_sc1,"mus1/fc/kb_fc")
data.frame(total_counts_per_cell = log10(kb_sc1$total),
           kept = kb_sc1$keep,
           label = kb_sc1$no_label,
           Mito_percent_per_cell= kb_sc1$subsets_Mito_percent,
           detected = kb_sc1$detected,
           preprocess="kallisto bustools") -> filter.kb

saveRDS(gene_filter(kb_sc1[,kb_sc1$keep]),file.path(write.path,"mus1_kb.rds"))




#salmon_alevin 

library(fishpond)
alevin_empty <- readRDS(file.path(empty.path,"alevin_mouse1.rds"))
library(tximport)
sal_ale_path<-"/stornext/General/data/user_managed/grpu_mritchie_1/Yue/preprocess/sal_ale/mouse/mus1/alevin_output/alevin/"
alevin <- tximport(file.path(sal_ale_path,"quants_mat.gz"),type="alevin")
alevin_filtered  <- alevin$counts[,which(alevin_empty$FDR<=0.001)]

alevin_sc <- SingleCellExperiment(
  assays = list(counts = as.matrix(alevin_filtered))
)
rownames(alevin_sc) <- unlist(strsplit(as.character(rownames(alevin_sc[])),".",fix=TRUE))[seq(1,nrow(alevin_sc)*2,2)]
barcode_info$SNG.1ST[match(colnames(alevin_sc), barcode_info$barcode)]-> alevin_sc$SNG.1ST

calculateQCMet(alevin_sc) ->alevin_sc
alevin_sc <- alevin_sc[rowSums(counts(alevin_sc))>0,]
#alevin_sc <- alevin_sc[,alevin_sc$total >100,]
data.frame(total_counts_per_gene =log10(rowSums(counts(alevin_sc))) ,
           pct.zeros= rowSums(counts(alevin_sc)==0)/ncol(alevin_sc),
           preprocess="salmon alevin") -> zero.d.alevin

alevin_sc$no_label <- !is.na(alevin_sc$SNG.1ST)

as.data.frame(rowData(alevin_sc)) %>% dplyr::group_by(gene_biotype) %>% summarise(count=n()) %>% mutate(preprocess="salmon alevin") ->gene_biotype_alevin


alevin_f<- scater_filter(alevin_sc)
#alevin_sc <- alevin_sc[rowSums(counts(alevin_sc))>0,]

alevin_sc$keep <-alevin_f
#plotfc(alevin_sc,"mus1/fc/alevin_scater_fc")
#plotfclabel(alevin_sc,"mus1/fc/alevin_fc")

data.frame(total_counts_per_cell = log10(alevin_sc$total),
           kept = alevin_sc$keep,
           label = alevin_sc$no_label,
           Mito_percent_per_cell= alevin_sc$subsets_Mito_percent,
           detected = alevin_sc$detected,
           preprocess="salmon alevin") -> filter.alevin

saveRDS(gene_filter(alevin_sc[,alevin_sc$keep]),file.path(write.path,"mus1_alevin.rds"))

#dropseqpipe

dropseq_path <- "/stornext/General/data/user_managed/grpu_mritchie_1/Yue/preprocess/dsp/mus1/results/summary/umi"
drop_mat <- readMM(file.path(dropseq_path,"matrix.mtx"))
drop_well <- read.table(file.path(dropseq_path,"barcodes.tsv"))
drop_genes <- read.table(file.path(dropseq_path,"genes.tsv"))

drop_genes <- mapIds(
  x = EnsDb.Mus.v99, 
  # NOTE: Need to remove gene version number prior to lookup.
  keys = as.character(drop_genes[,1]),
  keytype = "SYMBOL",
  column = "GENEID")
names(drop_genes) <- NULL
rownames(drop_mat) <- drop_genes
colnames(drop_mat) <- as.character(drop_well[,1])
drop_empty <- readRDS(file.path(empty.path,"dropseqpipe_mus1_emptied.rds"))

drop_mat <- drop_mat[,which(drop_empty$FDR<=0.001)]
drop_sc1<- SingleCellExperiment(
  assays = list(counts = as.matrix(drop_mat))
)
colnames(drop_sc1) <- unlist(strsplit(as.character(colnames(drop_sc1[])),"_",fix=TRUE))[seq(4,ncol(drop_sc1)*4,4)]
barcode_info$cell_ontology_class[match(colnames(drop_sc1), barcode_info$barcode)]-> drop_sc1$SNG.1ST


calculateQCMet(drop_sc1) ->drop_sc1

data.frame(total_counts_per_gene =log10(rowSums(counts(drop_sc1))) ,
           pct.zeros= rowSums(counts(drop_sc1)==0)/ncol(drop_sc1),
           preprocess="dropSeqPipe") -> zero.d.drop

drop_sc1$no_label <- !is.na(drop_sc1$SNG.1ST)
drop_sc1 <- drop_sc1[rowSums(counts(drop_sc1))>0,]
as.data.frame(rowData(drop_sc1)) %>% dplyr::group_by(gene_biotype) %>% summarise(count=n()) %>% mutate(preprocess="dropSeqPipe") ->gene_biotype_drop


drop_f<- scater_filter(drop_sc1)

drop_sc1$keep <-drop_f
#plotfc(drop_sc1,"mus1/fc/drop_scater_fc")
#plotfclabel(drop_sc1,"mus1/fc/drop_fc")

data.frame(total_counts_per_cell = log10(drop_sc1$total),
           kept = drop_sc1$keep,
           label = drop_sc1$no_label,
           Mito_percent_per_cell= drop_sc1$subsets_Mito_percent,
           detected = drop_sc1$detected,
           preprocess="dropSeqPipe") -> filter.drop
saveRDS(gene_filter(drop_sc1[,drop_sc1$keep]),file.path(write.path,"mus1_drop.rds"))



#cellranger
cellranger_empty <- readRDS(file.path(empty.path,"cellranger_mouse1_emptied.rds"))
cellr_path <-"/stornext/General/data/user_managed/grpu_mritchie_1/Yue/preprocess/cellranger/mus1/bamtofastq/outs/raw_feature_bc_matrix/"
cellrangr_sc <- read10xCounts(cellr_path,col.names=TRUE)
colnames(cellrangr_sc)<-sapply(strsplit(colnames(cellrangr_sc),"-"),"[[",1)
cellrangr_sc <- cellrangr_sc[,which(cellranger_empty$FDR<=0.001)]

barcode_info$cell_ontology_class[match(colnames(cellrangr_sc), barcode_info$barcode)]-> cellrangr_sc$SNG.1ST


calculateQCMet(cellrangr_sc) ->cellrangr_sc

data.frame(total_counts_per_gene =log10(rowSums(counts(cellrangr_sc))) ,
           pct.zeros= rowSums(counts(cellrangr_sc)==0)/ncol(cellrangr_sc),
           preprocess="Cell Ranger") -> zero.d.cellranger

cellrangr_sc$no_label <- !is.na(cellrangr_sc$SNG.1ST)
cellrangr_sc <- cellrangr_sc[rowSums(counts(cellrangr_sc))>0,]
as.data.frame(rowData(cellrangr_sc)) %>% dplyr::group_by(gene_biotype) %>% summarise(count=n()) %>% mutate(preprocess="Cell Ranger") ->gene_biotype_cellranger


cellrangr_f<- scater_filter(cellrangr_sc)

cellrangr_sc$keep <-cellrangr_f
#plotfc(cellrangr_sc,"mus1/fc/cellranger_scater_fc")
#plotfclabel(cellrangr_sc,"mus1/fc/cellranger_fc")

data.frame(total_counts_per_cell = log10(cellrangr_sc$total),
           kept = cellrangr_sc$keep,
           label = cellrangr_sc$no_label,
           Mito_percent_per_cell= cellrangr_sc$subsets_Mito_percent,
           detected = cellrangr_sc$detected,
           preprocess="Cell Ranger") -> filter.cellranger

saveRDS(gene_filter(cellrangr_sc[,cellrangr_sc$keep]),file.path(write.path,"mus1_cellranger.rds"))

#optimus
optimus_path <- "/stornext/General/data/user_managed/grpu_mritchie_1/Yue/preprocess/optimus/mus1/"
readRDS(file.path(optimus_path,"optimus.rds"))->optimus_sc
barcode_info$cell_ontology_class[match(colnames(optimus_sc), barcode_info$barcode)]-> optimus_sc$SNG.1ST

calculateQCMet(optimus_sc) ->optimus_sc

data.frame(total_counts_per_gene =log10(rowSums(counts(optimus_sc))) ,
           pct.zeros= rowSums(counts(optimus_sc)==0)/ncol(optimus_sc),
           preprocess="Optimus") -> zero.d.optimus

optimus_sc$no_label <- !is.na(optimus_sc$SNG.1ST)
optimus_sc <- optimus_sc[rowSums(counts(optimus_sc))>0,]
as.data.frame(rowData(optimus_sc)) %>% dplyr::group_by(gene_biotype) %>% summarise(count=n()) %>% mutate(preprocess="Optimus") ->gene_biotype_optimus


optimus_f<- scater_filter(optimus_sc)

optimus_sc$keep <-optimus_f
#plotfc(optimus_sc,"mus1/fc/optimus_scater_fc")
#plotfclabel(optimus_sc,"mus1/fc/optimus_fc")

data.frame(total_counts_per_cell = log10(optimus_sc$total),
           kept = optimus_sc$keep,
           label = optimus_sc$no_label,
           Mito_percent_per_cell= optimus_sc$subsets_Mito_percent,
           detected = optimus_sc$detected,
           preprocess="Optimus") -> filter.optimus

saveRDS(gene_filter(optimus_sc[,optimus_sc$keep]),file.path(write.path,"mus1_optimus.rds"))






###all together

#zero
bind_rows(zero.d.scpipe,zero.d.kb,zero.d.zumis,zero.d.alevin,zero.d.cellranger,zero.d.drop,zero.d.optimus) ->zero.d
saveRDS(zero.d,"SCE/raw/zeros/mus1_zero.rds")
png("results/raw/mus1/zero.png",units = "in",res=300,width = 6,height = 4)
ggplot(zero.d[zero.d$total_counts_per_gene>0,],aes(x=total_counts_per_gene,y=pct.zeros)) +facet_grid(.~preprocess)+
  geom_point(aes(col=preprocess),size=0.1) +scale_color_manual(values = mycol) +geom_smooth(col="black") +
  theme(text = element_text(size=20))+theme_bw()
dev.off()


#total counts density
pdf("results/raw/mus1/total_counts_density.pdf")
ggplot(zero.d) +geom_density(aes(x=total_counts_per_gene,color=preprocess),adjust = 1.2,size=1.1) +scale_color_manual(values = mycol) +
  theme(text = element_text(size=20),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())
dev.off()


#biotype
pdf("results/raw/mus1/biotype.pdf",width = 10)
bind_rows(gene_biotype_kb,gene_biotype_scpipe,gene_biotype_zumis,gene_biotype_alevin,gene_biotype_cellranger,gene_biotype_drop,
          gene_biotype_optimus) ->gene_biotype
saveRDS(gene_biotype,"SCE/raw/biotypes/mus1_gene_biotype.rds")


full_join(gene_biotype,gene_biotype %>% group_by(preprocess) %>% 
            summarise(n=sum(count)),by="preprocess") %>% mutate(proportion=count/n) %>% 
  group_by(preprocess) %>% top_n(10,proportion) -> gene_biotype

n <- length(unique(gene_biotype$gene_biotype))
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
col=sample(col_vector, n) 
col[7] <- "orange"
saveRDS(col,"SCE/raw/biotypes/col.rds")

col <- c( "#7570B3",col, "#66A61E" )
col <- c( "#FFD92F" ,col)
ggplot(gene_biotype[!is.na(gene_biotype$gene_biotype),]) + geom_col(aes(x=preprocess,y=proportion,fill=gene_biotype)) + 
  scale_fill_manual(values=col) +
  theme(text = element_text(size=20)) +theme_bw()
dev.off()

pdf("results/raw/mus1/biotype_count.pdf",width = 10)
ggplot(gene_biotype[!is.na(gene_biotype$gene_biotype),]) + geom_col(aes(x=preprocess,y=count,fill=gene_biotype)) + 
  scale_fill_manual(values=col) +
  theme(text = element_text(size=20)) +theme_bw()
dev.off()



#glmpca
library(glmpca)

plotglmpca <- function(sce,name){
  sce <- sce[,sce$keep=="TRUE"]
  sce <- sce[rowSums(counts(sce))>0,]
  res<-glmpca(counts(sce), 2)
  factors<-res$factors
  
  color <- as.character(sce$SNG.1ST)
  color[which(color=="H1975")] <-"red"
  color[which(color=="H2228")] <- "blue"
  color[which(color=="HCC827")] <- "green"
  color[is.na(color)] <-"grey"
  
  color[!colnames(sce) %in% allbc & color=="green"] <- "darkgreen"
  color[!colnames(sce) %in% allbc & color=="blue"] <- "darkblue"
  color[!colnames(sce) %in% allbc & color=="red"] <- "darkred"
  png(paste0(plotpath,"/",name,"all.png"),res=250,width = 5,height = 5, units = "in")
  plot(factors[,1],factors[,2],col= color,pch=19, cex = .8)
  dev.off()
  
  color <- as.character(sce$demuxlet_cls)
  color[which(color=="SNG")] <-"red"
  color[which(color=="DBL")] <- "green"
  color[is.na(color)] <-"grey"
  
  png(paste0(plotpath,"/",name,"db.png"),res=250,width = 5,height = 5, units = "in")
  plot(factors[,1],factors[,2],col= color,pch=19, cex = .8)
  dev.off()
  
}


plotglmpca(zumis_sc1,"plate_3cl/glm/zumis_glm")
plotglmpca(sc1,"plate_3cl/glm/scpipe_glm")
plotglmpca(scp_sc1_qc,"plate_3cl/glm/scpipe_rm_glm")
plotglmpca(celseq_sc1,"plate_3cl/glm/celseq_glm")
plotglmpca(scruff_sc1,"plate_3cl/glm/scruff_glm")



#filter 

bind_rows(filter.kb,filter.scpipe,filter.scpipe.rm,filter.zumis,filter.zumis.auto,filter.alevin,
          filter.cellranger,filter.optimus) ->filter.d
saveRDS(filter.d,"SCE/raw/filter/mus1_filter.rds")

library(ggsci)
#mito

pdf("results/raw/mus1/mito_filter.pdf")
ggplot(data=filter.d,aes(x=preprocess,y=Mito_percent_per_cell)) +
  geom_violin()  +
  ggbeeswarm::geom_quasirandom(size = 0.5, aes(colour = kept)) +
  coord_flip()+
  scale_y_log10() +scale_color_jco()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
        text=element_text(size=15)) +
  #scale_fill_manual(values = mycol)+
  theme_bw() + labs(x="",y="percentage of mitochondrial genes per cell",fill="preprocess")
dev.off()

#total
png("results/raw/mus1/total_filter.png")
ggplot(data=filter.d,aes(x=preprocess,y=total_counts_per_cell)) +
  geom_violin()  +
  facet_grid(.~kept) +
  ggbeeswarm::geom_quasirandom(size = 0.05, aes(colour = label)) +
  coord_flip()+scale_color_jco()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
        text=element_text(size=15)) +
  #scale_fill_manual(values = mycol)+
  theme_bw() + labs(x="",y="Total counts per cell",fill="preprocess")
dev.off()

#detected
png("results/raw/mus1/detect_filter.png")
ggplot(data=filter.d,aes(x=preprocess,y=detected)) +
  geom_violin()  +
  facet_grid(.~kept) +
  ggbeeswarm::geom_quasirandom(size = 0.05, aes(colour = label)) +
  coord_flip()+scale_color_jco()+scale_y_log10()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
        text=element_text(size=15)) +
  #scale_fill_manual(values = mycol)+
  theme_bw() + labs(x="",y="Number of detected genes per cell",fill="preprocess")
dev.off()

png("results/raw/mus1/kept_label_detected.png")
ggplot(data=filter.d[filter.d$kept,],aes(x=preprocess,y=detected)) +
  geom_violin() +
  #geom_point() + 
  facet_grid(.~label) +
  ggbeeswarm::geom_quasirandom(size = 0.1, aes(colour = label)) +
  coord_flip()+
  scale_y_log10() +scale_color_jco()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
        text=element_text(size=15)) +
  #scale_fill_manual(values = mycol)+
  theme_bw() + labs(x="",y="Number of detected genes per cell",fill="preprocess")
dev.off()
png("results/raw/mus1/kept_label_total.png")
ggplot(data=filter.d[filter.d$kept,],aes(x=preprocess,y=total_counts_per_cell)) +
  geom_violin() +
  #geom_point() + 
  facet_grid(.~label) +
  ggbeeswarm::geom_quasirandom(size = 0.1, aes(colour = label)) +
  coord_flip()+
  scale_y_log10() +scale_color_jco()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
        text=element_text(size=15)) +
  #scale_fill_manual(values = mycol)+
  theme_bw() + labs(x="",y="Total counts per cell",fill="preprocess")
dev.off()



#for bcs kept in all
intersect(intersect(intersect(colnames(sc1),colnames(zumis_sc1)),colnames(kb_sc1)),colnames(alevin_sc)) -> allbc
bind_rows(
  as.data.frame(rowData(sc1[,colnames(sc1)%in% allbc])) %>% dplyr::group_by(gene_biotype) %>% summarise(count=n()) %>% mutate(preprocess="scPipe"),
  as.data.frame(rowData(zumis_sc1[,colnames(zumis_sc1)%in% allbc])) %>% dplyr::group_by(gene_biotype) %>% summarise(count=n()) %>% mutate(preprocess="zUMIs"),
  as.data.frame(rowData(kb_sc1[,colnames(kb_sc1)%in% allbc])) %>% dplyr::group_by(gene_biotype) %>% summarise(count=n()) %>% mutate(preprocess="kallisto bustools"),
  as.data.frame(rowData(alevin_sc[,colnames(alevin_sc)%in% allbc])) %>% dplyr::group_by(gene_biotype) %>% summarise(count=n()) %>% mutate(preprocess="salmon alevin"),
)  ->gt
full_join(gt,gt %>% group_by(preprocess) %>% 
            summarise(n=sum(count)),by="preprocess") %>% mutate(proportion=count/n) %>% 
  group_by(preprocess) %>% top_n(10,proportion)->gt
ggplot(gt[!is.na(gt$gene_biotype),]) + geom_col(aes(x=preprocess,y=proportion,fill=gene_biotype)) + 
  scale_fill_manual(values=col) +
  theme(text = element_text(size=20)) +theme_bw()




#star

intersect(intersect(intersect(colnames(cellrangr_sc),colnames(zumis_sc1)),colnames(optimus_sc)),colnames(drop_sc1)) -> allbc
bind_rows(
  as.data.frame(rowData(cellrangr_sc[,colnames(cellrangr_sc)%in% allbc])) %>% dplyr::group_by(gene_biotype) %>% summarise(count=n()) %>% mutate(preprocess="Cell Ranger"),
  as.data.frame(rowData(zumis_sc1[,colnames(zumis_sc1)%in% allbc])) %>% dplyr::group_by(gene_biotype) %>% summarise(count=n()) %>% mutate(preprocess="zUMIs"),
  as.data.frame(rowData(optimus_sc[,colnames(optimus_sc)%in% allbc])) %>% dplyr::group_by(gene_biotype) %>% summarise(count=n()) %>% mutate(preprocess="Optimus"),
  as.data.frame(rowData(drop_sc1[,colnames(drop_sc1)%in% allbc])) %>% dplyr::group_by(gene_biotype) %>% summarise(count=n()) %>% mutate(preprocess="dropSeqPipe"),
)  ->gt
full_join(gt,gt %>% group_by(preprocess) %>% 
            summarise(n=sum(count)),by="preprocess") %>% mutate(proportion=count/n) %>% 
  group_by(preprocess) %>% top_n(10,proportion)->gt
ggplot(gt[!is.na(gt$gene_biotype),]) + geom_col(aes(x=preprocess,y=proportion,fill=gene_biotype)) + 
  scale_fill_manual(values=col) +
  theme(text = element_text(size=20)) +theme_bw()
ggplot(gt[!is.na(gt$gene_biotype),]) + geom_col(aes(x=preprocess,y=count,fill=gene_biotype)) + 
  scale_fill_manual(values=col) +
  theme(text = element_text(size=20)) +theme_bw()

bind_rows(data.frame(total_counts_per_gene =log10(rowSums(counts(sc1[,colnames(sc1)%in% allbc]))) ,
                     pct.zeros= rowSums(counts(sc1[,colnames(sc1)%in% allbc])==0)/ncol(sc1[,colnames(sc1)%in% allbc]),
                     preprocess="scPipe"),
          data.frame(total_counts_per_gene =log10(rowSums(counts(zumis_sc1[,colnames(zumis_sc1)%in% allbc]))) ,
                     pct.zeros= rowSums(counts(zumis_sc1[,colnames(zumis_sc1)%in% allbc])==0)/ncol(zumis_sc1[,colnames(zumis_sc1)%in% allbc]),
                     preprocess="zUMIs"),
          data.frame(total_counts_per_gene =log10(rowSums(counts(kb_sc1[,colnames(kb_sc1)%in% allbc]))) ,
                     pct.zeros= rowSums(counts(kb_sc1[,colnames(kb_sc1)%in% allbc])==0)/ncol(kb_sc1[,colnames(kb_sc1)%in% allbc]),
                     preprocess="kallisto bustools"),
          data.frame(total_counts_per_gene =log10(rowSums(counts(alevin_sc[,colnames(alevin_sc)%in% allbc]))) ,
                     pct.zeros= rowSums(counts(alevin_sc[,colnames(alevin_sc)%in% allbc])==0)/ncol(alevin_sc[,colnames(alevin_sc)%in% allbc]),
                     preprocess="salmon alevin")
          
) ->zd


ggplot(zd,aes(x=total_counts_per_gene,y=pct.zeros)) +facet_grid(.~preprocess)+
  geom_point(aes(col=preprocess),size=0.1) +scale_color_manual(values = mycol) +geom_smooth(col="black") +
  theme(text = element_text(size=20))+theme_bw()

ggplot(zd) +geom_density(aes(x=total_counts_per_gene,color=preprocess),adjust = 1.2,size=1.1) +scale_color_manual(values = mycol) +
  theme(text = element_text(size=20),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())
