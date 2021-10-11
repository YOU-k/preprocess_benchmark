designs=c("sc3cl","sc5cl","sc5clv3","pbmc5k","pbmc10k","mus1","mus2")
data=c("scPipe","zUMIs","kallisto bustools","Optimus","dropSeqPipe","Cell Ranger","salmon_SA","salmon_splici")
biotype=c("pro","lnc","pseudo")
setwd("/stornext/HPCScratch/home/you.y/preprocess_update/SCEs/disp/")
dir() -> files
files[-1] -> files
result0 <- lapply(files, function(file){readRDS(file)})

data.frame(cellnum=unlist(lapply(result0, function(x) length(x$trended.dispersion))),
           file=files) -> df


df$design="mus1"
df$data="a"
df$biotype="b"


df$design[grepl("mus1",files)] <- "mus1"
df$design[grepl("mus2",files)] <- "mus2"
df$design[grepl("sc3cl",files)] <- "sc3cl"
df$design[grepl("sc5cl",files)] <- "sc5cl"
df$design[grepl("sc5clv3",files)] <- "sc5clv3"
df$design[grepl("pbmc5k",files)] <- "pbmc5k"
df$design[grepl("pbmc10k",files)] <- "pbmc10k"


df$data[grepl("scPipe",files)] <- "scPipe"
df$data[grepl("zUMIs",files)] <- "zUMIs"
df$data[grepl("kallisto",files)] <- "kallisto bustools"
df$data[grepl("Optimus",files)] <- "Optimus"
df$data[grepl("dropSeqPipe",files)] <- "dropSeqPipe"
df$data[grepl("Cell Ranger",files)] <- "Cell Ranger"
df$data[grepl("salmon_SA",files)] <- "salmon_SA"
df$data[grepl("salmon_splici",files)] <- "salmon_splici"



df$biotype[grepl("pro",files)] <- "protein_coding"
df$biotype[grepl("lnc",files)] <- "lncRNA"
df$biotype[grepl("pseudo",files)] <- "pseudogene"

df
library(tidyverse)
df$design <- recode(df$design,"mus1"="10xv2_lungtissue1","mus2"="10xv2_lungtissue2","sc3cl"="10xv2_3cellline",
                          "sc5cl"="10xv2_5cellline","sc5clv3"="10xv3_5cellline","pbmc5k"="10xv3_pbmc5k",
                          "pbmc10k"="10xv3_pbmc10k")

library(ggplot2)
pdf("/stornext/HPCScratch/home/you.y/preprocess_update/results/droplet-based/biotype/cellnum_summary.pdf",height = 3,width = 12)
ggplot(df) + geom_col(aes(x=reorder(data,cellnum),y=cellnum,fill=data)) +
  facet_grid(biotype~design,scales = "free") +
  scale_fill_manual(values=droplet_col[unique(df$data)])+
  labs(x="",y="Number of cells",fill="Preprocessing methods") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),text = element_text(size=14),
        axis.ticks.x=element_blank(),panel.background = element_rect(fill = "white", colour = NA), 
        panel.border = element_rect(fill = NA,colour = "grey20"), panel.grid = element_line(colour = "grey92"), 
        panel.grid.minor = element_line(size = rel(0.5)), 
        strip.background = element_rect(fill = "grey85",colour = "grey20"))
dev.off()
