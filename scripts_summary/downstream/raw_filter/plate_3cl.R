library(scran)
library(scater)
library(CellBench)
library(R.utils)
library(tidyverse)

setwd("/stornext/General/data/user_managed/grpu_mritchie_1/Yue/preprocess/new")
data.path <- "/stornext/General/data/user_managed/grpu_mritchie_1/Yue/preprocess/data"
write.path <- "/stornext/General/data/user_managed/grpu_mritchie_1/Yue/preprocess/new/SCE/raw"
barcode_info <- read.csv(file.path(data.path,"barcode_annotation7.csv"))
plate_info <- read.table("/stornext/General/data/user_managed/grpu_mritchie_1/Yue/preprocess/data/demuxlet/sc_celseq2.metadata.csv",sep=",",header = TRUE)

plate_info <- plate_info[,c("cell_line" ,"cell_line_demuxlet" ,"demuxlet_cls")]
plate_info$SNG.1ST <- plate_info$cell_line
plate_info$samplename <- rownames(plate_info)
library(AnnotationHub)
ah <- AnnotationHub()
EnsDb.Hsapiens.v98 <- query(ah, c("EnsDb", "Homo Sapiens", 98))[[1]]

gene_filter <- function(sce){
  keep1 = (apply(counts(sce), 1, function(x) mean(x[x>0])) > 1)  # average count larger than 1
  keep2 = (rowSums(counts(sce)>0) > 3)  
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
  is_ercc <- grepl("ERCC",rownames(sce))
  sce <- addPerCellQC(sce, subsets = list(Mito=which(is_mito),ercc = which(is_ercc)))
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



plotpath <- "/stornext/General/data/user_managed/grpu_mritchie_1/Yue/preprocess/new/results/raw"
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
  ribo_set <- rowData(sce)$symbol[grep("^RP(S|L)", rowData(sce)$symbol)]
  
  is_mito <- rowData(sce)$symbol %in% mito_set
  is_ribo <- rowData(sce)$symbol %in% ribo_set
  png(paste0(plotpath,"/",name,".png"))
  
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
 
  write.table(text,paste0(plotpath,"/fc/",name,".txt"),quote = FALSE,row.names = FALSE,sep="\t")
  
}

##scpipe 

library(scPipe)
sc1_p <- "/stornext/General/data/user_managed/grpu_mritchie_1/Yue/preprocess/scpipe/nn84-4"

sc1 <- create_sce_by_dir(sc1_p, organism = "hsapiens_gene_ensembl", gene_id_type="ensembl_gene_id")
colData(sc1) = cbind(colData(sc1), DataFrame(plate_info[match(colnames(sc1),as.character(plate_info$samplename)),]))
#check match
#identical(as.character(plate_info[match(colnames(sc1),as.character(plate_info$wall_position)),][,3]),colnames(sc1))
sc1$barcode <- barcode_info$index[match(colnames(sc1),barcode_info$samplename)]


calculateQCMet(sc1) ->sc1
#zero
sc1 <- sc1[rowSums(counts(sc1))>0,]
data.frame(total_counts_per_gene =log10(rowSums(counts(sc1))) ,
           pct.zeros= rowSums(counts(sc1)==0)/ncol(sc1),
           preprocess="scPipe",
           biotype=rowData(sc1)$gene_biotype) -> zero.d.scpipe

#biotype
library(tidyverse)
as.data.frame(rowData(sc1)) %>% dplyr::group_by(gene_biotype) %>% summarise(count=n()) %>% mutate(preprocess="scPipe") ->gene_biotype_scpipe


#filter by remove outlier
sc1 <- calculate_QC_metrics(sc1)
scp_sc1_qc <- detect_outlier(sc1, type= "low", comp=2)
table(QC_metrics(scp_sc1_qc)$outliers)


scp_sc1_qc$keep <- FALSE
scp_sc1_qc$keep[scp_sc1_qc$outliers=="FALSE"] <- TRUE


#plotfc(scp_sc1_qc,"plate_3cl/scpipe_rmout_fc")
data.frame(total_counts_per_cell = log10(scp_sc1_qc$total),
           kept = scp_sc1_qc$keep,
           Mito_percent_per_cell= scp_sc1_qc$subsets_Mito_percent,
           detected = scp_sc1_qc$detected,
           preprocess="scPipe_rm") -> filter.scpipe.rm



#use scater
scater_filter(sc1) -> sc1_f
table(sc1_f)

sc1$keep <-sc1_f

saveRDS(gene_filter(sc1[,sc1$keep]),file.path(write.path,"plate3cl_scpipe.rds"))

plotfc(sc1,"plate_3cl/scpipe_scater_fc")

data.frame(total_counts_per_cell = log10(sc1$total),
           kept = sc1$keep,
           Mito_percent_per_cell= sc1$subsets_Mito_percent,
           detected = sc1$detected,
           preprocess="scPipe") -> filter.scpipe



#zumis

library(SingleCellExperiment)
zumis_path <- "/stornext/General/data/user_managed/grpu_mritchie_1/Yue/preprocess/zumis2/nn84-4/zUMIs_output/expression"
zumis_sc1_data <- readRDS(file.path(zumis_path,"SC1_S1.dgecounts.rds"))
count_sc1 <- as.matrix(zumis_sc1_data$umicount$inex$all)
colnames(count_sc1) <-barcode_info$samplename[match(colnames(zumis_sc1_data$readcount$exon$all),barcode_info$index)]
zumis_sc1 <- SingleCellExperiment(
  assays = list(counts = as.matrix(count_sc1)), 
  colData = plate_info[match(colnames(count_sc1),as.character(plate_info$samplename)),]
)

colnames(zumis_sc1) <- colnames(count_sc1)

zumis_sc1$barcode <- barcode_info$index[match(colnames(zumis_sc1),barcode_info$samplename)]

zumis_sc1 <- zumis_sc1[rowSums(counts(zumis_sc1))>0,]
calculateQCMet(zumis_sc1) ->zumis_sc1
data.frame(total_counts_per_gene =log10(rowSums(counts(zumis_sc1))) ,
           pct.zeros= rowSums(counts(zumis_sc1)==0)/ncol(zumis_sc1),
           preprocess="zUMIs",
           biotype=rowData(zumis_sc1)$gene_biotype) -> zero.d.zumis



zumis_f<- scater_filter(zumis_sc1)

as.data.frame(rowData(zumis_sc1)) %>% dplyr::group_by(gene_biotype) %>% summarise(count=n()) %>% mutate(preprocess="zUMIs") ->gene_biotype_zumis


zumis_sc1$keep <-zumis_f
#plotfc(zumis_sc1,"plate_3cl/zumis_scater_fc")


#only exon 
count_sc2 <- as.matrix(zumis_sc1_data$umicount$exon$all)
colnames(count_sc2) <-barcode_info$samplename[match(colnames(zumis_sc1_data$readcount$exon$all),barcode_info$index)]
zumis_sc2 <- SingleCellExperiment(
  assays = list(counts = as.matrix(count_sc2)),
  colData = plate_info[match(colnames(count_sc2),as.character(plate_info$samplename)),]
)
colnames(zumis_sc2) <- colnames(count_sc2)

zumis_sc2 <- zumis_sc2[rowSums(counts(zumis_sc2))>0,]
data.frame(total_counts_per_gene =log10(rowSums(counts(zumis_sc2))) ,
           pct.zeros= rowSums(counts(zumis_sc2)==0)/ncol(zumis_sc2),
           preprocess="zUMIs_exon") -> zero.d.zumis2
calculateQCMet(zumis_sc2) ->zumis_sc2

zumis_sc2$barcode <- barcode_info$index[match(colnames(zumis_sc2),barcode_info$samplename)]

colData = plate_info[match(colnames(count_sc2),as.character(plate_info$samplename)),]

as.data.frame(rowData(zumis_sc2)) %>% dplyr::group_by(gene_biotype) %>% summarise(count=n()) %>% mutate(preprocess="zUMIs_exon") ->gene_biotype_zumis2

zumis_f<- scater_filter(zumis_sc1)
zumis_sc1$keep <-zumis_f
#plotfc(zumis_sc1,"plate_3cl/zumis_scater_fc")

data.frame(total_counts_per_cell = log10(zumis_sc1$total),
           kept = zumis_sc1$keep,
           Mito_percent_per_cell= zumis_sc1$subsets_Mito_percent,
           detected = zumis_sc1$detected,
           preprocess="zUMIs") -> filter.zumis

saveRDS(gene_filter(zumis_sc1[,zumis_sc1$keep]),file.path(write.path,"plate3cl_zumis.rds"))

#celseq2
cs_path <-"/stornext/General/data/user_managed/grpu_mritchie_1/Yue/preprocess/celseq2/nn84-4/result/expr/nn84-4"
cs_mat <- read.csv(file.path(cs_path,"expr.csv"))
rownames(cs_mat)<-cs_mat$X
cs_mat <- cs_mat[,-1]
colnames(cs_mat)<- unlist(strsplit(colnames(cs_mat),".",fix=TRUE))[seq(3,ncol(cs_mat)*3,3)]
colnames(cs_mat) <-barcode_info$samplename[match(colnames(cs_mat),barcode_info$index)]
celseq_sc1 <- SingleCellExperiment(
  assays = list(counts = as.matrix(cs_mat)), 
  colData = plate_info[match(colnames(cs_mat),as.character(plate_info$samplename)),]
)

colnames(celseq_sc1) <- colnames(cs_mat)

celseq_sc1 <- celseq_sc1[rowSums(counts(celseq_sc1))>0,]
calculateQCMet(celseq_sc1) ->celseq_sc1
data.frame(total_counts_per_gene =log10(rowSums(counts(celseq_sc1))) ,
           pct.zeros= rowSums(counts(celseq_sc1)==0)/ncol(celseq_sc1),
           preprocess="celseq2",
           biotype=rowData(celseq_sc1)$gene_biotype) -> zero.d.celseq2


celseq_sc1$barcode <- barcode_info$index[match(colnames(celseq_sc1),barcode_info$samplename)]
celseq_f<- scater_filter(celseq_sc1)

as.data.frame(rowData(celseq_sc1)) %>% dplyr::group_by(gene_biotype) %>% summarise(count=n())  %>% mutate(preprocess="celseq2") ->gene_biotype_celseq2


celseq_sc1$keep <-celseq_f
plotfc(celseq_sc1,"plate_3cl/celseq_scater_fc")

data.frame(total_counts_per_cell = log10(celseq_sc1$total),
           kept = celseq_sc1$keep,
           Mito_percent_per_cell= celseq_sc1$subsets_Mito_percent,
           detected = celseq_sc1$detected,
           preprocess="celseq2") -> filter.celseq

saveRDS(gene_filter(celseq_sc1[,celseq_sc1$keep]),file.path(write.path,"plate3cl_celseq.rds"))



###scruff
library(scruff)
scruff.path <- "/stornext/General/data/user_managed/grpu_mritchie_1/Yue/preprocess/scruff/nn84-4/Count"
load(file.path(scruff.path,"20200818_072310_countUMI_sce.rda"))
scruff_sc1 <- SingleCellExperiment(
  assays = list(counts = counts(scruffsce)))
scruff_sc1$barcode <- scruffsce$barcode

scruff_sc1$samplename <- barcode_info$samplename[match(scruff_sc1$barcode,barcode_info$index)]
colData(scruff_sc1) = cbind(colData(scruff_sc1),plate_info[match(scruff_sc1$samplename,as.character(plate_info$samplename)),])

colnames(scruff_sc1) <- scruff_sc1$samplename


calculateQCMet(scruff_sc1) ->scruff_sc1

scruff_sc1 <- scruff_sc1[rowSums(counts(scruff_sc1))>0,]
data.frame(total_counts_per_gene =log10(rowSums(counts(scruff_sc1))) ,
           pct.zeros= rowSums(counts(scruff_sc1)==0)/ncol(scruff_sc1),
           preprocess="scruff",
           biotype=rowData(scruff_sc1)$gene_biotype) -> zero.d.scruff
scruff_f<- scater_filter(scruff_sc1)

as.data.frame(rowData(scruff_sc1)) %>% dplyr::group_by(gene_biotype) %>% summarise(count=n())  %>% mutate(preprocess="scruff") ->gene_biotype_scruff

scruff_sc1$keep <-scruff_f
plotfc(scruff_sc1,"plate_3cl/scruff_scater_fc")

data.frame(total_counts_per_cell = log10(scruff_sc1$total),
           kept = scruff_sc1$keep,
           Mito_percent_per_cell= scruff_sc1$subsets_Mito_percent,
           detected = scruff_sc1$detected,
           preprocess="scruff") -> filter.scruff

saveRDS(gene_filter(scruff_sc1[,scruff_sc1$keep]),file.path(write.path,"plate3cl_scruff.rds"))

#kallisto 

library(Matrix)
kb_path <-"/stornext/General/data/user_managed/grpu_mritchie_1/Yue/preprocess/kallisto_bus/nn84-4/genecount"
kb_mat <- as.data.frame(t(as.matrix(readMM(file=file.path(kb_path,"genes.mtx")))))
kb_well <- read.table(file.path(kb_path,"genes.barcodes.txt"))
kb_genes <- read.table(file.path(kb_path,"genes.genes.txt"))
rownames(kb_mat) <- kb_genes[,1]
colnames(kb_mat) <- kb_well[,1]

colnames(kb_mat) <-barcode_info$samplename[match(colnames(kb_mat),barcode_info$index)]
kb_sc1 <- SingleCellExperiment(
  assays = list(counts = as.matrix(kb_mat)), 
  colData = plate_info[match(colnames(kb_mat),as.character(plate_info$samplename)),]
)
colnames(kb_sc1) <- colnames(kb_mat)
rownames(kb_sc1) <- unlist(strsplit(as.character(rownames(kb_sc1[])),".",fix=TRUE))[seq(1,nrow(kb_sc1)*2,2)]

calculateQCMet(kb_sc1) ->kb_sc1

kb_sc1$barcode <- barcode_info$index[match(colnames(kb_sc1),barcode_info$samplename)]

kb_sc1 <- kb_sc1[rowSums(counts(kb_sc1))>0,]
data.frame(total_counts_per_gene =log10(rowSums(counts(kb_sc1))) ,
           pct.zeros= rowSums(counts(kb_sc1)==0)/ncol(kb_sc1),
           preprocess="kallisto bustools",
           biotype=rowData(kb_sc1)$gene_biotype) -> zero.d.kb


as.data.frame(rowData(kb_sc1)) %>% dplyr::group_by(gene_biotype) %>% summarise(count=n()) %>% mutate(preprocess="kallisto bustools") ->gene_biotype_kb

kb_sc1 <- kb_sc1[,colSums(counts(kb_sc1))>0,]
kb_f<- scater_filter(kb_sc1)
kb_sc1$keep <-kb_f
plotfc(kb_sc1,"plate_3cl/kb_scater_fc")


data.frame(total_counts_per_cell = log10(kb_sc1$total),
           kept = kb_sc1$keep,
           Mito_percent_per_cell= kb_sc1$subsets_Mito_percent,
           detected = kb_sc1$detected,
           preprocess="kallisto bustools") -> filter.kb

saveRDS(gene_filter(kb_sc1[,kb_sc1$keep]),file.path(write.path,"plate3cl_kb.rds"))


###all together

#zero
bind_rows(zero.d.celseq2,zero.d.scpipe,zero.d.kb,zero.d.scruff,zero.d.zumis) ->zero.d
saveRDS(zero.d,"SCE/raw/zeros/plate_3cl_zero.rds")

pdf("results/raw/plate_3cl/zero_all.pdf",width = 5,height = 4)
ggplot(zero.d[zero.d$total_counts_per_gene>0,],aes(x=10^(total_counts_per_gene),y=pct.zeros)) +
  scale_color_manual(values = mycol) +
  scale_x_log10()+
  geom_smooth(aes(col=preprocess),span = 0.8) +
  theme(text = element_text(size=20))+theme_bw()
dev.off()


pdf("results/raw/plate_3cl/zero.pdf",width = 8,height = 4)
ggplot(zero.d[zero.d$total_counts_per_gene>0,],aes(x=total_counts_per_gene,y=pct.zeros)) +
  facet_grid(.~preprocess)+
  geom_point(aes(col=preprocess),size=0.5) +
  scale_color_manual(values = mycol) +
  geom_smooth(span = 0.8,color="black") +
  theme(text = element_text(size=20))+theme_bw()
dev.off()


#total counts density
pdf("results/raw/plate_3cl/total_counts_density.pdf",width = 5,height = 4)
ggplot(zero.d) +geom_density(aes(x=total_counts_per_gene,color=preprocess),adjust = 1.2,size=1.1) +scale_color_manual(values = mycol) +
  theme(text = element_text(size=20)) +theme_bw()
dev.off()


#biotype
bind_rows(gene_biotype_kb,gene_biotype_scpipe,gene_biotype_zumis,gene_biotype_scruff,gene_biotype_celseq2) ->gene_biotype
saveRDS(gene_biotype,"SCE/raw/biotypes/plate_3cl_gene_biotype.rds")


full_join(gene_biotype,gene_biotype %>% group_by(preprocess) %>% 
            summarise(n=sum(count)),by="preprocess") %>% mutate(proportion=count/n) %>% 
  group_by(preprocess) %>% top_n(7,proportion) -> gene_biotype

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
pdf("results/raw/plate_3cl/biotype_count.pdf",width =6, height = 4)
ggplot(gene_biotype[!is.na(gene_biotype$gene_biotype),]) + 
  geom_col(aes(x=reorder(preprocess, count),y=count,fill=reorder(gene_biotype,count))) + 
  scale_fill_manual(values=col1) +
  theme(text = element_text(size=12),axis.text.x = element_text(angle = 30,hjust = 1),
        panel.background = element_rect(fill = "white"),axis.line = element_line(colour = "black")) 
dev.off()



#for bcs kept in all
intersect(intersect(colnames(sc1)[sc1$keep],colnames(scruff_sc1)[scruff_sc1$keep]),colnames(celseq_sc1)[celseq_sc1$keep]) -> allbc
intersect(colnames(zumis_sc1)[zumis_sc1$keep],allbc) -> allbc
#glmpca
library(glmpca)

plotglmpca <- function(sce,name){
  sce$demuxlet_cls <- as.character(sce$demuxlet_cls)
  sce$demuxlet_cls[is.na(sce$demuxlet_cls)] <- "no"
  sce <- sce[,!sce$demuxlet_cls=="DBL"]
  sce <- sce[,sce$keep=="TRUE"]
  sce <- sce[rowSums(counts(sce))>0,]
  res<-glmpca(counts(sce), 2,minibatch='stochastic')
  factors<-res$factors
  
  color <- as.character(sce$SNG.1ST)
  color[which(color=="H1975")] <-"red"
  color[which(color=="H2228")] <- "blue"
  color[which(color=="HCC827")] <- "green"
  color[is.na(color)] <-"grey"
  
  #color[!colnames(sce) %in% allbc & color=="green"] <- "darkgreen"
  #color[!colnames(sce) %in% allbc & color=="blue"] <- "darkblue"
  #color[!colnames(sce) %in% allbc & color=="red"] <- "darkred"
  plotpath="/stornext/General/data/user_managed/grpu_mritchie_1/Yue/preprocess_update/results/plate_based/"
  #png(paste0(plotpath,"/",name,"all.png"),res=250,width = 5,height = 5, units = "in")
  png(paste0(plotpath,"/","kallisto_plate3cl.png"),res=250,width = 5,height = 5, units = "in")
  plot(factors[,1],factors[,2],col= color,pch=19, cex = .8)
  dev.off()
}


plotglmpca(zumis_sc1,"plate_3cl/glm/zumis_glm")
plotglmpca(sc1,"plate_3cl/glm/scpipe_glm")
#plotglmpca(scp_sc1_qc,"plate_3cl/glm/scpipe_rm_glm")
plotglmpca(celseq_sc1,"plate_3cl/glm/celseq_glm")
plotglmpca(scruff_sc1,"plate_3cl/glm/scruff_glm")



###upsetr plot

datasets <- list(`scPipe` = sc1,
                 `scruff`=scruff_sc1,
                 `zUMIs`= zumis_sc1,`celseq2` = celseq_sc1,`kallisto bustools` = kb_sc1)


get_barcodes <- function(sce){
  colnames(sce)[sce$keep]
}
list(bc=get_barcodes) -> bcs
apply_methods(datasets,bcs) -> bc_all
library(UpSetR)
bc_all$result -> rowname
names(rowname) <- bc_all$data
pdf("results/raw/plate_3cl/upsetr.pdf",width = 5,height = 5)
upset(fromList(rowname),nsets = 7,order.by="freq",nintersects = 10, sets.bar.color = "#56B4E9",matrix.color = "blue")
dev.off()

###all cells violin plots
get_cells_lib <- function(sce){
  return(sce$sum)
}
get_cells_gene <- function(sce){
  return(sce$detected)
}

get_cells_n <- function(sce){
  return(ncol(sce))
}
list(lib= get_cells_lib,
     genes= get_cells_gene,
     n= get_cells_n) -> libm


library(CellBench)
apply_methods(datasets,libm) -> lib
pdf("results/raw/plate_3cl/total_counts_boxplot.pdf",width = 5,height = 3.5)
ggplot(unnest(lib %>% spread(libm,result))) +geom_boxplot(aes(y=lib,x=data,color=data)) +scale_color_manual(values = mycol)+
  geom_text(data=lib %>% spread(libm,result),aes(y=3.5e+05,x=data,label=n)) +  
  theme(text = element_text(size=12),axis.text.x = element_text(angle = 30,hjust = 1),
        panel.background = element_rect(fill = "white"),axis.line = element_line(colour = "black")) 
dev.off()

apply_methods(datasets,libm) -> lib
pdf("results/raw/plate_3cl/genes_boxplot.pdf",width = 5,height = 3.5)
ggplot(unnest(lib %>% spread(libm,result))) +geom_boxplot(aes(y=genes,x=data,color=data)) +scale_color_manual(values = mycol)+ 
  theme(text = element_text(size=12),axis.text.x = element_text(angle = 30,hjust = 1),
        panel.background = element_rect(fill = "white"),axis.line = element_line(colour = "black")) 
dev.off()


#for bcs kept in all
intersect(intersect(intersect(colnames(sc1),colnames(zumis_sc1)),colnames(celseq_sc1)),colnames(scruff_sc1)) -> allbc

intersect(intersect(intersect(rownames(sc1),rownames(zumis_sc1)),rownames(celseq_sc1)),rownames(scruff_sc1)) -> allgenes


get_overlap_counts <- function(sce){
  sce[rownames(sce) %in% allgenes,colnames(sce) %in% allbc]  -> sce
  sce[match(allgenes,rownames(sce)),match(allbc,colnames(sce))] ->sce
  return(as.matrix(counts(sce)))
}

get_overlap_cells <- function(sce){
  sce[match(allgenes,rownames(sce)),match(allbc,colnames(sce))] ->sce
  return(as.matrix(counts(sce)))
}

datasets$`kallisto bustools` <- NULL
lapply(datasets, get_overlap_counts) -> counts_all

lapply(datasets, get_overlap_cells) -> cells_all


lapply(cells_all, colSums) -> umi_counts
bind_cols(umi_counts) -> umi_counts


adjustcolor(mycol,alpha.f = 0.8) ->mycola
names(mycola) <- names(mycol)
mycola[1] <- mycola[7]

#colnames(umi_counts) <- c("scPipe","Cell Ranger","kallisto bustools","salmon alevin","Optimus")



pdf("results/raw/plate_3cl/total_counts_overlap.pdf",width = 5,height = 3.5)
ggplot(umi_counts) + geom_point(aes(x=scPipe,y=scruff,color="scruff"),size=0.8) + 
  geom_point(aes(x=scPipe,y=celseq2,color="celseq2"),size=0.8) +
  geom_point(aes(x=scPipe,y=zUMIs,color="zUMIs"),size=0.8) +
  scale_color_manual(values=mycola) + 
  geom_abline(slope=1, intercept=0) +theme(text = element_text(size=20))+theme_bw()
dev.off()


lapply(cells_all, function(x){colSums(x>0)}) -> dgenes
bind_cols(dgenes) -> dgenes

pdf("results/raw/plate_3cl/genes_overlap.pdf",width = 5,height = 3.5)
ggplot(dgenes) + geom_point(aes(x=scPipe,y=scruff,color="scruff"),size=0.8) + 
  geom_point(aes(x=scPipe,y=celseq2,color="celseq2"),size=0.8) +
  geom_point(aes(x=scPipe,y=zUMIs,color="zUMIs"),size=0.8) +
  scale_color_manual(values=mycola) + 
  geom_abline(slope=1, intercept=0) +theme(text = element_text(size=20))+theme_bw()
dev.off()


saveRDS(dgenes,"SCE/raw/overlap/plate_3cl_dgenes.rds")
saveRDS(umi_counts,"SCE/raw/overlap/plate_3cl_umi_counts.rds")

calculate_cor <- function(pp){
  cor = lapply(1:length(counts_all), 
               function (x){diag(cor(counts_all[[pp]],counts_all[[x]],method = "pearson"))})
  return(unlist(cor))}


cors <- data.frame(preprocess = rep(names(counts_all),each=length(allbc)),
                   correlation = calculate_cor(1))
saveRDS(cors,"SCE/raw/cor/plate_3cl_cor_pearson.rds")



#add kb
intersect(colnames(sc1),colnames(kb_sc1)) -> bc_kb
intersect(rownames(sc1),rownames(kb_sc1)) -> gene_kb

get_overlap_counts <- function(sce){
  sce[rownames(sce) %in% gene_kb,colnames(sce) %in% bc_kb]  -> sce
  sce[match(gene_kb,rownames(sce)),match(bc_kb,colnames(sce))] ->sce
  return(as.matrix(counts(sce)))
}

dd <- list(`scPipe` = sc1,`kallisto bustools` = kb_sc1)

lapply(dd, get_overlap_counts) -> counts_all

cors <- data.frame(preprocess = rep(names(counts_all),each=length(bc_kb)),
                   correlation = calculate_cor(1))

bind_rows(cors[cors$preprocess=="kallisto bustools",],plate_3cl_cor_pearson) -> cors
## 1 above was to specify scpipe
pdf("results/raw/plate_3cl/correlation.pdf",width = 5,height = 4)
ggplot(cors[!cors$preprocess=="scPipe",]) + geom_boxplot(aes(x=reorder(preprocess, -correlation),y=correlation,color=preprocess)) +
  scale_color_manual(values = mycol) + theme(text = element_text(size=12),axis.text.x = element_text(angle = 15),
                                             panel.background = element_rect(fill = "white"),axis.line = element_line(colour = "black"))  +
  labs(x="",y="Pearson Correlation coefficient")
dev.off()
saveRDS(cors,"SCE/raw/cor/plate_3cl_cor_pearson.rds")

saveRDS(umi_counts,"SCE/raw/cor/plate_3cl_cors.rds")

#filter 

bind_rows(filter.celseq,filter.kb,filter.scpipe,filter.scpipe.rm,filter.scruff,filter.zumis) ->filter.d
saveRDS(filter.d,"SCE/raw/filter/plate_3cl_filter.rds")

library(ggsci)
#mito

pdf("results/raw/plate_3cl/mito_filter.pdf")
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
pdf("results/raw/plate_3cl/total_filter.pdf",width = 6,height = 4)
ggplot(data=filter.d,aes(x=preprocess,y=total_counts_per_cell)) +
  geom_violin()  +facet_grid(.~kept)+
  ggbeeswarm::geom_quasirandom(size = 0.5, aes(colour = kept)) +
  coord_flip()+
  geom_text(data=filter.d %>% group_by(preprocess,kept) %>% summarise(n=n()) ,aes(x=preprocess,y=0.5,label=n),size=3)+
  scale_y_log10() +scale_color_jco()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
        text=element_text(size=15)) +
  #scale_fill_manual(values = mycol)+
  theme_bw() + labs(x="",y="Total counts per cell",fill="preprocess")
dev.off()

#detected
pdf("results/raw/plate_3cl/detected_filter.pdf",width = 6,height = 4)
ggplot(data=filter.d,aes(x=preprocess,y=detected)) +
  geom_violin()  +
  ggbeeswarm::geom_quasirandom(size = 0.5, aes(colour = kept)) +
  coord_flip()+facet_grid(.~kept) +
  scale_y_log10() +scale_color_jco()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
        text=element_text(size=15)) +
  geom_text(data=filter.d %>% group_by(preprocess,kept) %>% summarise(n=n()) ,aes(x=preprocess,y=5,label=n),size=3)+
  #scale_fill_manual(values = mycol)+
  theme_bw() + labs(x="",y="Number of detected genes per cell",fill="preprocess")
dev.off()











#kb compared to scpipe
identical(colnames(sc1),colnames(kb_sc1))
data.frame(total_counts=sc1$total,
           detetced_g = sc1$detected,
           bc=colnames(sc1),
           keep=sc1$keep) -> scpipe_info
data.frame(total_countsk=kb_sc1$total,
           detetced_gk = kb_sc1$detected,
           bc=colnames(kb_sc1),
           keepk=kb_sc1$keep) -> kb_info
full_join(kb_info,scpipe_info) -> p
library(ggpubr)

pdf("../plots_added/plate_3cl_kb_scpipe_total.pdf",width=8,height = 3)

p1 <- ggplot(p[!is.na(p$detetced_gk),]) + geom_point(aes(x=log10(total_counts),y=log10(total_countsk),col=keep)) +theme_bw() +
  scale_color_brewer(palette = "Set1") +
  labs(x="log10(scPipe)",y="log10(kallisto butsools)",title="Total counts per cell",col="Kept in scPipe")
p2 <- ggplot(p[!is.na(p$detetced_gk),]) + geom_point(aes(x=log10(total_counts),y=log10(total_countsk),col=keepk)) +theme_bw() +
  scale_color_brewer(palette = "Set1") +
  labs(x="log10(scPipe)",y="log10(kallisto butsools)",title=" ",col="Kept in kb")
ggarrange(p1,p2)

dev.off()


pdf("../plots_added/plate_3cl_kb_scpipe_gene.pdf",width=8,height = 3)

p1 <- ggplot(p[!is.na(p$detetced_gk),]) + geom_point(aes(x=log10(detetced_g),y=log10(detetced_gk),col=keep)) +theme_bw() +
  scale_color_brewer(palette = "Set1") +
  labs(x="log10(scPipe)",y="log10(kallisto butsools)",title="Number of detected genes per cell",col="Kept in scPipe")
p2 <- ggplot(p[!is.na(p$detetced_gk),]) + geom_point(aes(x=log10(detetced_g),y=log10(detetced_gk),col=keepk)) +theme_bw() +
  scale_color_brewer(palette = "Set1") +
  labs(x="log10(scPipe)",y="log10(kallisto butsools)",title=" ",col="Kept in kb")
ggarrange(p1,p2)

dev.off()





ggplot(unnest(lib %>% spread(libm,result))) +geom_violin(aes(y=log10(lib),x=data,color=data)) +scale_color_manual(values = mycol)+
  geom_hline(yintercept= log10(4096)) +
  theme(text = element_text(size=12),axis.text.x = element_text(angle = 30,hjust = 1),
        panel.background = element_rect(fill = "white"),axis.line = element_line(colour = "black")) 

