setwd("/stornext/HPCScratch/home/you.y/preprocess_update/SCEs/glmpc_sil/")
result0 <- lapply(dir(), function(file){readRDS(file)})

result0
library(tidyverse)
bind_rows(result0) -> result1
result1$dt=result1$design
result1$dt[grepl("lung",result1$dt)] <- "Lung_tissue"
result1$dt[grepl("cell",result1$dt)] <- "Cell_line"
result1$dt[grepl("pbmc",result1$dt)] <- "PBMC"

as.character(result1$glmpca_sil) -> result1$glmpca_sil
result1$glmpca_sil[grepl("pseu",result1$glmpca_sil)] <- "Pseudogene"
result1$glmpca_sil[grepl("lnc",result1$glmpca_sil)] <- "lncRNA"
result1$glmpca_sil[grepl("pro",result1$glmpca_sil)] <- "protein_coding"
table(result1$dt)

library(ggplot2)
pdf("/stornext/HPCScratch/home/you.y/preprocess_update/results/droplet-based/biotype/glmpc_sil.pdf",width = 7,height = 4)
ggplot(result1,aes(x=reorder(data,result),y=result,col=data)) + geom_boxplot() + 
  geom_jitter() +
  facet_grid(glmpca_sil~dt,scales = "free_x") +coord_flip() +
  scale_color_manual(values=droplet_col)+
  geom_hline(yintercept=0,linetype="dotted") +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),text = element_text(size=14),
        axis.ticks.y=element_blank(),panel.background = element_rect(fill = "white", colour = NA), 
        panel.border = element_rect(fill = NA,colour = "grey20"), panel.grid = element_line(colour = "grey92"), 
        panel.grid.minor = element_line(size = rel(0.5)), 
        strip.background = element_rect(fill = "grey85",colour = "grey20")) +
  labs(y="Silhouette widths",col="Preprocessing ")
dev.off()
