library(scran)
library(tidyverse)
library(scater)
library(CellBench)
SCE.path <- "/stornext/General/data/user_managed/grpu_mritchie_1/Yue/preprocess/new/SCE/raw"


#droplet tissue
design=c("mus2")
dataname=c("kb","scpipe","zumis","alevin","optimus","drop","cellranger")
paste(apply(expand.grid(design, dataname), 1, paste, collapse="_"), ".rds",sep="") ->files
#check
files %in% dir(SCE.path)

result0 <- lapply(files, function(file){readRDS(file.path(SCE.path,file))})

names(result0) <- dataname
Reduce(intersect, lapply(result0, colnames)) -> allbc
Reduce(intersect, lapply(result0, rownames)) -> allgenes

get_overlap_counts <- function(sce){
  sce[rownames(sce) %in% allgenes,colnames(sce) %in% allbc]  -> sce
  sce[match(allgenes,rownames(sce)),match(allbc,colnames(sce))] ->sce
  return(as.matrix(counts(sce)))
}

lapply(result0, get_overlap_counts) -> counts_all


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

cors_use$preprocess1 <- recode(cors_use$preprocess1,"cellranger"="Cell Ranger","drop"="dropSeqPipe",
                               "optimus" ="Optimus","alevin"="salmon alevin","zumis"="zUMIs","scpipe"="scPipe","kb"="kallisto bustools")


cors_use$preprocess2 <- recode(cors_use$preprocess2,"cellranger"="Cell Ranger","drop"="dropSeqPipe",
                               "optimus" ="Optimus","alevin"="salmon alevin","zumis"="zUMIs","scpipe"="scPipe","kb"="kallisto bustools")

setwd("/stornext/General/data/user_managed/grpu_mritchie_1/Yue/preprocess/new/SCE")
saveRDS(cors_use,"raw/cor/all/mus2.rds")


cors_use$preprocess1_f <- factor(cors_use$preprocess1,levels=c("Cell Ranger","dropSeqPipe","Optimus",
                                                               "salmon alevin","zUMIs","scPipe"))

pdf("../results/raw/cors/mus2.pdf",width=6,heigh=5)
full_join(cors_use  %>% group_by(preprocess1_f,preprocess2) %>% summarise(mean=mean(correlation)), cors_use) %>% 
  ggplot() + geom_boxplot(aes(y=correlation)) + facet_grid(preprocess1_f ~ preprocess2) +
  geom_text(aes(0,0.2,label=round(mean,2)),size=3.5) + labs(y="Pearson correlation coefficient")
theme_bw() 
dev.off()


#sc3cl



design=c("sc3cl")
dataname=c("kb","scpipe","zumis","alevin","optimus","drop","cellranger")
paste(apply(expand.grid(design, dataname), 1, paste, collapse="_"), ".rds",sep="") ->files
#check
files %in% dir(SCE.path)

result0 <- lapply(files, function(file){readRDS(file.path(SCE.path,file))})

names(result0) <- dataname
Reduce(intersect, lapply(result0, colnames)) -> allbc
Reduce(intersect, lapply(result0, rownames)) -> allgenes

get_overlap_counts <- function(sce){
  sce[rownames(sce) %in% allgenes,colnames(sce) %in% allbc]  -> sce
  sce[match(allgenes,rownames(sce)),match(allbc,colnames(sce))] ->sce
  return(as.matrix(counts(sce)))
}

lapply(result0, get_overlap_counts) -> counts_all


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

cors_use$preprocess1 <- recode(cors_use$preprocess1,"cellranger"="Cell Ranger","drop"="dropSeqPipe",
                               "optimus" ="Optimus","alevin"="salmon alevin","zumis"="zUMIs","scpipe"="scPipe","kb"="kallisto bustools")


cors_use$preprocess2 <- recode(cors_use$preprocess2,"cellranger"="Cell Ranger","drop"="dropSeqPipe",
                               "optimus" ="Optimus","alevin"="salmon alevin","zumis"="zUMIs","scpipe"="scPipe","kb"="kallisto bustools")

setwd("/stornext/General/data/user_managed/grpu_mritchie_1/Yue/preprocess/new/SCE")
saveRDS(cors_use,"raw/cor/all/sc3cl.rds")


cors_use$preprocess1_f <- factor(cors_use$preprocess1,levels=c("Cell Ranger","dropSeqPipe","Optimus",
                                                               "salmon alevin","zUMIs","scPipe"))

pdf("../results/raw/cors/sc3cl.pdf",width=6,heigh=5)
full_join(cors_use  %>% group_by(preprocess1_f,preprocess2) %>% summarise(mean=mean(correlation)), cors_use) %>% 
  ggplot() + geom_boxplot(aes(y=correlation)) + facet_grid(preprocess1_f ~ preprocess2) +
  geom_text(aes(0,0.8,label=round(mean,2)),size=3.5) + labs(y="Pearson correlation coefficient") +
  theme_bw() 
dev.off()




#sc5clv3


design=c("sc5clv3")
dataname=c("kb","scpipe","zumis","alevin","optimus","drop","cellranger")
paste(apply(expand.grid(design, dataname), 1, paste, collapse="_"), ".rds",sep="") ->files
#check
files %in% dir(SCE.path)

result0 <- lapply(files, function(file){readRDS(file.path(SCE.path,file))})

names(result0) <- dataname
Reduce(intersect, lapply(result0, colnames)) -> allbc
Reduce(intersect, lapply(result0, rownames)) -> allgenes

get_overlap_counts <- function(sce){
  sce[rownames(sce) %in% allgenes,colnames(sce) %in% allbc]  -> sce
  sce[match(allgenes,rownames(sce)),match(allbc,colnames(sce))] ->sce
  return(as.matrix(counts(sce)))
}

lapply(result0, get_overlap_counts) -> counts_all


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

cors_use$preprocess1 <- recode(cors_use$preprocess1,"cellranger"="Cell Ranger","drop"="dropSeqPipe",
                               "optimus" ="Optimus","alevin"="salmon alevin","zumis"="zUMIs","scpipe"="scPipe","kb"="kallisto bustools")


cors_use$preprocess2 <- recode(cors_use$preprocess2,"cellranger"="Cell Ranger","drop"="dropSeqPipe",
                               "optimus" ="Optimus","alevin"="salmon alevin","zumis"="zUMIs","scpipe"="scPipe","kb"="kallisto bustools")

setwd("/stornext/General/data/user_managed/grpu_mritchie_1/Yue/preprocess/new/SCE")
saveRDS(cors_use,"raw/cor/all/sc5clv3.rds")


cors_use$preprocess1_f <- factor(cors_use$preprocess1,levels=c("Cell Ranger","dropSeqPipe","Optimus",
                                                               "salmon alevin","zUMIs","scPipe"))

pdf("../results/raw/cors/sc5clv3.pdf",width=6,heigh=5)
full_join(cors_use  %>% group_by(preprocess1_f,preprocess2) %>% summarise(mean=mean(correlation)), cors_use) %>% 
  ggplot() + geom_boxplot(aes(y=correlation)) + facet_grid(preprocess1_f ~ preprocess2) +
  geom_text(aes(0,0.7,label=round(mean,2)),size=3.5) + labs(y="Pearson correlation coefficient") +
  theme_bw() 
dev.off()



