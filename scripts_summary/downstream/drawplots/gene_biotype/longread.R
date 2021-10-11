#only common cells should be used
read.csv("/stornext/General/data/user_managed/grpu_mritchie_1/Yue/flames/PromethION_5cl_rep1/FLTSA_output/transcript_count.csv")-> long_matrix
long_matrix$gene_id

long_matrix <- long_matrix[,-1]
aggregate(. ~ gene_id, long_matrix, sum) -> long_matrix

rownames(long_matrix) <- long_matrix$gene_id
long_matrix <- long_matrix[,-1]
unlist(strsplit(as.character(rownames(long_matrix[])),".",fix=TRUE))[seq(1,nrow(long_matrix)*2,2)] -> long_matrix$rowname
aggregate(. ~ rowname, long_matrix, sum) -> long_matrix

rownames(long_matrix) <- long_matrix$rowname

long_matrix$rowname<- NULL

library(AnnotationHub)
ah <- AnnotationHub()
EnsDb.Hsapiens.v98 <- query(ah, c("EnsDb", "Homo Sapiens", 98))[[1]]


keep1 = (apply(long_matrix, 1, function(x) mean(x[x>0])) > 1)  # average count larger than 1
keep2 = (rowSums(long_matrix>0) > 10)  
table(keep1&keep2)
long_matrix = long_matrix[(keep1 & keep2), ]

biotype<- mapIds(
  x = EnsDb.Hsapiens.v98, 
  # NOTE: Need to remove gene version number prior to lookup.
  keys = rownames(long_matrix),
  keytype = "GENEID",
  column = "GENEBIOTYPE")
table(grepl("pseudo",biotype))[2]/nrow(long_matrix) -> pct_1
sum(long_matrix[grepl("pseudo",biotype),])/sum(long_matrix) -> pct_2
tibble(data="longread",read_prop=pct_1,detect_prop=pct_2) -> long_d


SCE.path <- "/stornext/HPCScratch/home/you.y/preprocess_update/SCEs/raw/droplet_based/"
design=c("sc5cl")
dataname=c("kb","scpipe","zumis","drop","cellranger","sa","splici")
paste(apply(expand.grid(design, dataname), 1, paste, collapse="_"), ".rds",sep="") ->files


#check
files %in% dir(SCE.path)

result0 <- lapply(files, function(file){readRDS(file.path(SCE.path,file))})
rep(dataname,each=length(design))->data

tibble(design=rep(design,length(dataname)), 
       data=data,
       result=result0) -> datasets

Reduce(intersect,lapply(datasets$result, function(sce){colnames(sce)})) -> coln
intersect(coln,colnames(long_matrix)) -> use_cell


lapply(datasets$result, function(sce){ sce[,use_cell]}) -> datasets$result
reads_propor <- function(sce){
  pct <- sum(counts(sce)[grepl("pseudo",rowData(sce)$gene_biotype),])/sum(counts(sce))
  return(pct)
}

detect_propor <- function(sce) {
  pct <- table(grepl("pseudo",rowData(sce)$gene_biotype))[2] /nrow(sce)
  return(pct)
}


lapply(datasets$result, reads_propor) -> datasets$read_prop
lapply(datasets$result, detect_propor) -> datasets$detect_prop


long_matrix2 <- long_matrix[,use_cell]
table(grepl("pseudo",biotype))[2]/nrow(long_matrix2) -> pct_1
sum(long_matrix[grepl("pseudo",biotype),])/sum(long_matrix2) -> pct_2
tibble(data="longread",read_prop=pct_1,detect_prop=pct_2) -> long_d


rbind(datasets %>% dplyr::select(data,read_prop,detect_prop),long_d)  -> tmp
tmp %>% mutate(long= (data=="longread")) -> tmp

tmp$data <- recode(tmp$data,  "scpipe"="scPipe","zumis"="zUMIs","kb"="kallisto bustools",
                   "optimus" ="Optimus","drop"="dropSeqPipe","cellranger"="Cell Ranger",
                   "sa"="salmon_SA","splici"="salmon_splici")

droplet_col
c(droplet_col,"darkslategrey") -> mycol2
names(mycol2)[9] <- "longread"

saveRDS(tmp,"/stornext/HPCScratch/home/you.y/preprocess_update/SCEs/bcv/longread_sc5clv2.rds")
tmp %>%  mutate(data = fct_reorder(data, as.numeric(read_prop))) %>% 
  dplyr::filter(!data=="Cell Ranger") %>% 
  ggplot() + geom_col(aes(x=data,y=read_prop,fill=data))  +scale_fill_manual(values=mycol2) +theme_bw() +
  labs(x="",y="Proportions (reads mapped to pseudogenes)", fill="Preprocessing workflows")
ggsave("/stornext/General/data/user_managed/grpu_mritchie_1/Yue/flames/preprocess/new/reads_v2.pdf",width=4,height = 2)
tmp %>%  mutate(data = fct_reorder(data, as.numeric(detect_prop))) %>% 
  ggplot() + geom_col(aes(x=data,y=detect_prop,fill=long)) +scale_fill_brewer(palette = "Set2") +theme_bw() +
  labs(x="",y="Proportions (detectd pseudogene features)", fill="Long_read")
ggsave("/stornext/General/data/user_managed/grpu_mritchie_1/Yue/flames/preprocess/new/detect_v2.pdf",width=4,height = 2)







##10xv3

read.csv("/stornext/General/data/user_managed/grpu_mritchie_1/Yue/flames/PromethION_5cl_rep2/FLTSA_output/transcript_count.csv")-> long_matrix
long_matrix$gene_id

long_matrix <- long_matrix[,-1]
aggregate(. ~ gene_id, long_matrix, sum) -> long_matrix

rownames(long_matrix) <- long_matrix$gene_id
long_matrix <- long_matrix[,-1]
unlist(strsplit(as.character(rownames(long_matrix[])),".",fix=TRUE))[seq(1,nrow(long_matrix)*2,2)] -> long_matrix$rowname
aggregate(. ~ rowname, long_matrix, sum) -> long_matrix

rownames(long_matrix) <- long_matrix$rowname

long_matrix$rowname<- NULL

library(AnnotationHub)
ah <- AnnotationHub()
EnsDb.Hsapiens.v98 <- query(ah, c("EnsDb", "Homo Sapiens", 98))[[1]]


keep1 = (apply(long_matrix, 1, function(x) mean(x[x>0])) > 1)  # average count larger than 1
keep2 = (rowSums(long_matrix>0) > 10)  
table(keep1&keep2)
long_matrix = long_matrix[(keep1 & keep2), ]



SCE.path <- "/stornext/HPCScratch/home/you.y/preprocess_update/SCEs/raw/droplet_based/"
design=c("sc5clv3")
dataname=c("kb","scpipe","zumis","drop","cellranger","sa","splici")
paste(apply(expand.grid(design, dataname), 1, paste, collapse="_"), ".rds",sep="") ->files


#check
files %in% dir(SCE.path)

result0 <- lapply(files, function(file){readRDS(file.path(SCE.path,file))})
rep(dataname,each=length(design))->data

tibble(design=rep(design,length(dataname)), 
       data=data,
       result=result0) -> datasets

Reduce(intersect,lapply(datasets$result, function(sce){colnames(sce)})) -> coln
intersect(coln,colnames(long_matrix)) -> use_cell


lapply(datasets$result, function(sce){ sce[,use_cell]}) -> datasets$result
reads_propor <- function(sce){
  pct <- sum(counts(sce)[grepl("pseudo",rowData(sce)$gene_biotype),])/sum(counts(sce))
  return(pct)
}

detect_propor <- function(sce) {
  pct <- table(grepl("pseudo",rowData(sce)$gene_biotype))[2] /nrow(sce)
  return(pct)
}

lapply(datasets$result, reads_propor) -> datasets$read_prop
lapply(datasets$result, detect_propor) -> datasets$detect_prop


long_matrix2 <- long_matrix[,use_cell]
table(grepl("pseudo",biotype))[2]/nrow(long_matrix2) -> pct_1
sum(long_matrix[grepl("pseudo",biotype),])/sum(long_matrix2) -> pct_2
tibble(data="longread",read_prop=pct_1,detect_prop=pct_2) -> long_d


rbind(datasets %>% dplyr::select(data,read_prop,detect_prop),long_d)  -> tmp
tmp %>% mutate(long= (data=="longread")) -> tmp


tmp$data <- recode(tmp$data,  "scpipe"="scPipe","zumis"="zUMIs","kb"="kallisto bustools",
                   "optimus" ="Optimus","drop"="dropSeqPipe","cellranger"="Cell Ranger",
                   "sa"="salmon_SA","splici"="salmon_splici")



saveRDS(tmp,"/stornext/HPCScratch/home/you.y/preprocess_update/SCEs/bcv/longread_sc5clv3.rds")
tmp %>%  mutate(data = fct_reorder(data, as.numeric(read_prop))) %>% 
  dplyr::filter(!data=="Cell Ranger") %>% 
  ggplot() + geom_col(aes(x=data,y=read_prop,fill=data))  +scale_fill_manual(values=mycol2) +theme_bw() +
  labs(x="",y="Proportions (reads mapped to pseudogenes)", fill="Preprocessing workflows")
ggsave("/stornext/General/data/user_managed/grpu_mritchie_1/Yue/flames/preprocess/new/reads_v3.pdf",width=4,height = 2)
tmp %>%  mutate(data = fct_reorder(data, as.numeric(detect_prop))) %>% 
  ggplot() + geom_col(aes(x=data,y=detect_prop,fill=long)) +scale_fill_brewer(palette = "Set2") +theme_bw() +
  labs(x="",y="Proportions (detectd pseudogene features)", fill="Long_read")
ggsave("/stornext/General/data/user_managed/grpu_mritchie_1/Yue/flames/preprocess/detect_v3.pdf",width=4,height = 2)

