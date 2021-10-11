write.path <- "/stornext/HPCScratch/home/you.y/preprocess_update/results/droplet-based/"

#biotype_count
cal_biotype <- function(sce){
  return(table(rowData(sce)$gene_biotype))
}
designs=c("pbmc5k","pbmc10k","sc5cl","sc3cl","sc5clv3","mus1","mus2")

biotype_all <- tibble()

for (d in designs) {
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
lapply(datasets$result,cal_biotype) -> biotype_n

data.frame(count=unlist(biotype_n),
           data=rep(datasets$data,unlist(lapply(biotype_n, length))),
           biotype=unlist(names(unlist(biotype_n))),
           design=unique(datasets$design)) -> df


df %>% dplyr::filter(count>50) -> df_filter

left_join(df_filter,df_filter %>% group_by(data) %>% summarise(n=sum(count))) -> df_filter

df_filter$data <- recode(df_filter$data, "scpipe"="scPipe","zumis"="zUMIs","kb"="kallisto bustools",
                              "optimus" ="Optimus","drop"="dropSeqPipe","cellranger"="Cell Ranger",
                              "alevin"="salmon alevin")
pdf(paste0(write.path,"biotype/",d,"_bar.pdf"),width=7,height=4)
p <- ggplot(df_filter) + geom_col(aes(x=reorder(data,n),y=count,fill=biotype)) + 
  scale_fill_manual(values=mycol_biotype[unique(df_filter$biotype)]) +
  theme(text = element_text(size=14),
        axis.ticks.x=element_blank(),panel.background = element_rect(fill = "white", colour = NA), 
        panel.border = element_rect(fill = NA,colour = "grey20"), panel.grid = element_line(colour = "grey92"), 
        panel.grid.minor = element_line(size = rel(0.5)), 
        strip.background = element_rect(fill = "grey85",colour = "grey20"))+
  scale_x_discrete(guide = guide_axis(angle = 15)) +
  labs(x="",y="Number of detected features",title=unique(df$design),fill="Gene biotypes")
print(p)
dev.off()
#biotype_umi_count_pct
calumi_count <- function(sce){
  sce[!is.na(rowData(sce)$gene_biotype),] ->sce
  bt <- unique(rowData(sce)$gene_biotype)
  umi<- c()
  for (i in bt){
    sum(as.matrix(counts(sce))[rowData(sce)$gene_biotype==i,]) -> tmp_umi
    umi<-c(umi,tmp_umi)}
  data.frame(biotype=bt,n=umi) -> df
  return(df)
}

lapply(datasets$result,calumi_count) -> biotype_umi

bind_cols(bind_rows(biotype_umi),
          data.frame(data=rep(datasets$data,unlist(lapply(biotype_umi, nrow))))) -> df


df$biotype[grepl("IG",df$biotype) | grepl("TR",df$biotype) |
             grepl("TEC",df$biotype)] <- "others"
df$biotype[grepl("pseudogene",df$biotype)] <- "pseudogene"
df$biotype[!df$biotype %in% c("pseudogene","protein_coding","snoRNA","snRNA","lncRNA","miRNA","misc_RNA")] <- "others"

df %>% group_by(biotype,data) %>% summarise(n=sum(n)) -> df

left_join(df, df %>% group_by(data) %>% summarise(all=sum(n))) %>% mutate(pct=n/all) -> biotype_tmp
biotype_tmp$design <- unique(datasets$design)
bind_rows(biotype_tmp,biotype_all) -> biotype_all
}



biotype_all <- biotype_all[!is.na(biotype_all$pct),]
biotype_all$design <- factor(biotype_all$design,levels = c("10xv2_3cellline","10xv2_5cellline","10xv3_5cellline","10xv2_lungtissue1","10xv2_lungtissue2",
                                                           "10xv3_pbmc5k","10xv3_pbmc10k"))

bind_rows(biotype_all1 %>% filter(!design=="10xv2_lungtissue1"), biotype_all) -> biotype_all2
table(biotype_all2$data,biotype_all2$design)

saveRDS(biotype_all2,"/stornext/HPCScratch/home/you.y/preprocess_update/SCEs/raw/summary/biotype_all.rds")


pdf(file.path(write.path,"biotype/summary_pct_umi.pdf"),width = 13,height = 6.5)
ggplot(biotype_all2) +geom_col(aes(x=data,y=pct*100,fill=data),position = "dodge") +
  facet_grid(reorder(biotype,-pct)~design,scales = "free") +
  scale_fill_manual(values = droplet_col) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),text = element_text(size=14),
        axis.ticks.x=element_blank(),panel.background = element_rect(fill = "white", colour = NA), 
        panel.border = element_rect(fill = NA,colour = "grey20"), panel.grid = element_line(colour = "grey92"), 
        panel.grid.minor = element_line(size = rel(0.5)), 
        strip.background = element_rect(fill = "grey85",colour = "grey20"))+
  labs(x="",y="Percentage (UMI counts)",fill="Preprocessing workflows")
dev.off()          






#plate_3cl
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
       data=droplet_recode(data),
       result=result0) -> datasets
lapply(datasets$result,cal_biotype) -> biotype_n

data.frame(count=unlist(biotype_n),
           data=rep(datasets$data,unlist(lapply(biotype_n, length))),
           biotype=unlist(names(unlist(biotype_n))),
           design=unique(datasets$design)) -> df


df %>% dplyr::filter(count>50) -> df_filter

left_join(df_filter,df_filter %>% group_by(data) %>% summarise(n=sum(count))) -> df_filter

df_filter$data <- recode(df_filter$data, "scpipe"="scPipe","zumis"="zUMIs","kb"="kallisto bustools","celseq"="celseq2")
pdf(paste0(write.path,"biotype/",d,"_bar.pdf"),width=7,height=4)
p <- ggplot(df_filter) + geom_col(aes(x=reorder(data,n),y=count,fill=biotype)) + 
  scale_fill_manual(values=mycol_biotype[unique(df_filter$biotype)]) +
  theme(text = element_text(size=14),
        axis.ticks.x=element_blank(),panel.background = element_rect(fill = "white", colour = NA), 
        panel.border = element_rect(fill = NA,colour = "grey20"), panel.grid = element_line(colour = "grey92"), 
        panel.grid.minor = element_line(size = rel(0.5)), 
        strip.background = element_rect(fill = "grey85",colour = "grey20"))+
  scale_x_discrete(guide = guide_axis(angle = 15)) +
  labs(x="",y="Number of detected features",title="plate_3cellline",fill="Gene biotypes")
print(p)
dev.off()
