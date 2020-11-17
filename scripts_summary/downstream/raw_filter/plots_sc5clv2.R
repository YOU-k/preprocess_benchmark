#plots
###all together

#zero
bind_rows(zero.d.scpipe,zero.d.kb,zero.d.zumis,zero.d.alevin,zero.d.cellranger,zero.d.optimus,zero.d.drop) ->zero.d
saveRDS(zero.d,"SCE/raw/zeros/sc5clv2_zero.rds")

pdf("results/raw/sc_5clv2/zero_all.pdf",width = 5,height = 4)
ggplot(zero.d[zero.d$total_counts_per_gene>0,],aes(x=10^(total_counts_per_gene),y=pct.zeros)) +
  scale_color_manual(values = mycol) +
  scale_x_log10()+
  geom_smooth(aes(col=preprocess),span = 0.8) +
  theme(text = element_text(size=20))+theme_bw()
dev.off()

pdf("results/raw/sc_5clv2/zero.pdf",width = 10,height = 4)
ggplot(zero.d[zero.d$total_counts_per_gene>0,],aes(x=total_counts_per_gene,y=pct.zeros)) +
  facet_grid(.~preprocess)+
  geom_point(aes(col=preprocess),size=0.5) +
  scale_color_manual(values = mycol) +
  geom_smooth(span = 0.8,color="black") +
  theme(text = element_text(size=20))+theme_bw()
dev.off()


#total counts density
pdf("results/raw/sc_5clv2/total_counts_density.pdf",width = 5,height = 4)
ggplot(zero.d) +geom_density(aes(x=total_counts_per_gene,color=preprocess),adjust = 1.2,size=1.1) +scale_color_manual(values = mycol) +
  theme(text = element_text(size=20)) +theme_bw()
dev.off()


#biotype
pdf("results/raw/sc_5clv2/biotype.pdf",width = 10)
bind_rows(gene_biotype_kb,gene_biotype_scpipe,gene_biotype_zumis,gene_biotype_alevin,gene_biotype_cellranger,
          gene_biotype_optimus,gene_biotype_drop) ->gene_biotype
saveRDS(gene_biotype,"SCE/raw/biotypes/sc5clv2_gene_biotype.rds")


full_join(gene_biotype,gene_biotype %>% group_by(preprocess) %>% 
            summarise(n=sum(count)),by="preprocess") %>% mutate(proportion=count/n) %>% 
  group_by(preprocess) %>% top_n(7,proportion) -> gene_biotype

col1 <- c(col1,"violetred","red")
names(col1)[13:14] <- c("TR_C_gene","IG_C_gene")
pdf("results/raw/sc_5clv2/biotype_count.pdf",width =6, height = 4)
ggplot(gene_biotype[!is.na(gene_biotype$gene_biotype),]) + 
  geom_col(aes(x=reorder(preprocess, count),y=count,fill=reorder(gene_biotype,count))) + 
  scale_fill_manual(values=col1) +
  theme(text = element_text(size=12),axis.text.x = element_text(angle = 30,hjust = 1),
        panel.background = element_rect(fill = "white"),axis.line = element_line(colour = "black")) 
dev.off()



#glmpca
library(glmpca)
plotpath <- "/stornext/General/data/user_managed/grpu_mritchie_1/Yue/preprocess/new/results/raw/"

intersect(intersect(colnames(sc1)[sc1$keep],colnames(zumis_sc1)[zumis_sc1$keep]),colnames(drop_sc1)[drop_sc1$keep]) -> allbcf
intersect(intersect(colnames(kb_sc1)[kb_sc1$keep],intersect(colnames(cellrangr_sc)[cellrangr_sc$keep],allbcf)),
          colnames(optimus_sc)[optimus_sc$keep]) -> allbcf


intersect(allbcf,colnames(alevin_sc)[alevin_sc$keep]) -> allbcf

plotglmpca <- function(sce,name){
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
  color[which(color=="A549")] <- "orange"
  color[which(color=="H8383")] <- "purple"
  
  color[is.na(color)] <-"grey"
  
  color[!colnames(sce) %in% allbcf & color=="green"] <- "darkgreen"
  color[!colnames(sce) %in% allbcf & color=="blue"] <- "darkblue"
  color[!colnames(sce) %in% allbcf & color=="red"] <- "darkred"
  color[!colnames(sce) %in% allbcf & color=="orange"] <- "orange4"
  color[!colnames(sce) %in% allbcf & color=="purple"] <- "purple4"
  png(paste0(plotpath,"/",name,"all.png"),res=250,width = 5,height = 5, units = "in")
  plot(factors[,1],factors[,2],col= color,pch=19, cex = .8)
  dev.off()
}


plotglmpca(zumis_sc1,"sc_5clv2/glm/zumis_glm")
plotglmpca(sc1,"sc_5clv2/glm/scpipe_glm")
plotglmpca(cellrangr_sc,"sc_5clv2/glm/cellranger_glm")
plotglmpca(kb_sc1,"sc_5clv2/glm/kb_glm")
plotglmpca(alevin_sc,"sc_5clv2/glm/alevin_glm")
plotglmpca(optimus_sc,"sc_5clv2/glm/optimus_glm")
plotglmpca(drop_sc1,"sc_5clv2/glm/drop_glm")


###upsetr plot
datasets <- list(`scPipe` = sc1,
                 `Cell Ranger` = cellrangr_sc,
                 `kallisto bustools` = kb_sc1,
                 `salmon alevin` = alevin_sc,
                 `Optimus`= optimus_sc,
                 `dropSeqPipe` = drop_sc1,
                 `zUMIs`= zumis_sc1)
get_barcodes <- function(sce){
  colnames(sce)[sce$keep]
}
list(bc=get_barcodes) -> bcs

apply_methods(datasets,bcs) -> bc_all
library(UpSetR)
bc_all$result -> rowname
names(rowname) <- bc_all$data
pdf("results/raw/sc_5clv2/upsetr.pdf",width = 5,height = 5)
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
pdf("results/raw/sc_5clv2/total_counts_violin.pdf",width = 5,height = 3.5)
ggplot(unnest(lib %>% spread(libm,result))) +geom_violin(aes(y=lib,x=data,color=data)) +scale_color_manual(values = mycol)+
  geom_text(data=lib %>% spread(libm,result),aes(y=1.5e+05,x=data,label=n)) +  
  theme(text = element_text(size=12),axis.text.x = element_text(angle = 30,hjust = 1),
        panel.background = element_rect(fill = "white"),axis.line = element_line(colour = "black")) 
dev.off()

apply_methods(datasets,libm) -> lib
pdf("results/raw/sc_5clv2/genes_violin.pdf",width = 5,height = 3.5)
ggplot(unnest(lib %>% spread(libm,result))) +geom_violin(aes(y=genes,x=data,color=data)) +scale_color_manual(values = mycol)+ 
  theme(text = element_text(size=12),axis.text.x = element_text(angle = 30,hjust = 1),
        panel.background = element_rect(fill = "white"),axis.line = element_line(colour = "black")) 
dev.off()


#for bcs kept in all
intersect(intersect(intersect(colnames(sc1),colnames(zumis_sc1)),colnames(kb_sc1)),colnames(alevin_sc)) -> allbc

intersect(intersect(intersect(allbc,colnames(optimus_sc)),colnames(cellrangr_sc)),colnames(drop_sc1)) -> allbc

intersect(intersect(intersect(rownames(sc1),rownames(zumis_sc1)),rownames(kb_sc1)),rownames(alevin_sc)) -> allgenes

intersect(intersect(intersect(allgenes,rownames(optimus_sc)),rownames(cellrangr_sc)),rownames(drop_sc1)) -> allgenes


get_overlap_counts <- function(sce){
  sce[rownames(sce) %in% allgenes,colnames(sce) %in% allbc]  -> sce
  sce[match(allgenes,rownames(sce)),match(allbc,colnames(sce))] ->sce
  return(as.matrix(counts(sce)))
}

get_overlap_cells <- function(sce){
  sce[match(allgenes,rownames(sce)),match(allbc,colnames(sce))] ->sce
  return(as.matrix(counts(sce)))
}

lapply(datasets, get_overlap_counts) -> counts_all

lapply(datasets, get_overlap_cells) -> cells_all


lapply(cells_all, colSums) -> umi_counts
bind_cols(umi_counts) -> umi_counts


adjustcolor(mycol,alpha.f = 0.8) ->mycola
names(mycola) <- names(mycol)
mycola[1] <- mycola[7]

#colnames(umi_counts) <- c("scPipe","Cell Ranger","kallisto bustools","salmon alevin","Optimus")



pdf("results/raw/sc_5clv2/total_counts_overlap.pdf",width = 5,height = 3.5)
ggplot(umi_counts) + geom_point(aes(x=`Cell Ranger`,y=scPipe,color="scPipe"),size=0.8) + 
  geom_point(aes(x=`Cell Ranger`,y=`kallisto bustools`,color="kallisto bustools"),size=0.8) +
  geom_point(aes(x=`Cell Ranger`,y=`salmon alevin`,color="salmon alevin"),size=0.8) +
  geom_point(aes(x=`Cell Ranger`,y=Optimus,color="Optimus"),size=0.8) +
  geom_point(aes(x=`Cell Ranger`,y=dropSeqPipe,color="dropSeqPipe"),size=0.8) +
  geom_point(aes(x=`Cell Ranger`,y=zUMIs,color="zUMIs"),size=0.8) +
  scale_color_manual(values=mycola) + 
  geom_abline(slope=1, intercept=0) +theme(text = element_text(size=20))+theme_bw()
dev.off()


lapply(cells_all, function(x){colSums(x>0)}) -> dgenes
bind_cols(dgenes) -> dgenes

pdf("results/raw/sc_5clv2/genes_overlap.pdf",width = 5,height = 3.5)
ggplot(dgenes) + geom_point(aes(x=`Cell Ranger`,y=scPipe,color="scPipe"),size=0.8) + 
  geom_point(aes(x=`Cell Ranger`,y=`kallisto bustools`,color="kallisto bustools"),size=0.8) +
  geom_point(aes(x=`Cell Ranger`,y=`salmon alevin`,color="salmon alevin"),size=0.8) +
  geom_point(aes(x=`Cell Ranger`,y=Optimus,color="Optimus"),size=0.8) +
  geom_point(aes(x=`Cell Ranger`,y=dropSeqPipe,color="dropSeqPipe"),size=0.8) +
  geom_point(aes(x=`Cell Ranger`,y=zUMIs,color="zUMIs"),size=0.8) +
  scale_color_manual(values=mycola) + theme(text = element_text(size=20))+theme_bw() +
  geom_abline(slope=1, intercept=0)
dev.off()


saveRDS(dgenes,"SCE/raw/overlap/sc5clv2_dgenes.rds")
saveRDS(umi_counts,"SCE/raw/overlap/sc5clv2_umi_counts.rds")

calculate_cor <- function(pp){
  cor = lapply(1:length(counts_all), 
               function (x){diag(cor(counts_all[[pp]],counts_all[[x]],method = "pearson"))})
  return(unlist(cor))}


cors <- data.frame(preprocess = rep(names(counts_all),each=length(allbc)),
                   correlation = calculate_cor(2))
## 2 above was to specify cellranger 
pdf("results/raw/sc_5clv2/correlation.pdf",width = 5,height = 3.5)
ggplot(cors[!cors$preprocess=="Cell Ranger",]) + geom_boxplot(aes(x=reorder(preprocess, -correlation),y=correlation,color=preprocess)) +
  scale_color_manual(values = mycol) + theme(text = element_text(size=12),axis.text.x = element_text(angle = 30,hjust = 1),
                                             panel.background = element_rect(fill = "white"),axis.line = element_line(colour = "black")) 
dev.off()

saveRDS(umi_counts,"SCE/raw/cor/sc5clv2_cors.rds")



#filter 

bind_rows(filter.kb,filter.scpipe,filter.scpipe.rm,filter.zumis,filter.zumis.auto,filter.alevin,
          filter.cellranger,filter.optimus,filter.drop) ->filter.d
saveRDS(filter.d,"SCE/raw/filter/sc5clv2_filter.rds")

library(ggsci)

#mito

pdf("results/raw/sc_5clv2/mito_filter.pdf")
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
pdf("results/raw/sc_5clv2/total_filter.pdf",width = 6,height = 4)

ggplot(data=filter.d,aes(x=preprocess,y=total_counts_per_cell)) +
  geom_violin()  +
  facet_grid(.~kept) +
  ggbeeswarm::geom_quasirandom(size = 0.05, aes(colour = kept)) +
  coord_flip()+scale_color_jco()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
        text=element_text(size=15)) +
  geom_text(data=filter.d %>% group_by(preprocess,kept) %>% summarise(n=n()) ,aes(x=preprocess,y=1,label=n),size=3) +
  #scale_fill_manual(values = mycol)+
  theme_bw() + labs(x="",y="Total counts per cell",fill="preprocess")
dev.off()


pdf("results/raw/sc_5clv2/total_filter_density.pdf",width = 6,height = 4)
ggplot(data=filter.d,aes(x=total_counts_per_cell,color=preprocess)) +
  geom_density(adjust = 0.5,size=1)  +
  #facet_grid(.~kept) +
  #ggbeeswarm::geom_quasirandom(size = 0.05, aes(colour = kept)) +
  #coord_flip()+scale_color_jco()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
        text=element_text(size=15)) +
  scale_color_manual(values = mycol)+
  theme_bw() + labs(x="",y="Total counts per cell",fill="preprocess")
dev.off()

#detected
pdf("results/raw/sc_5clv2/detected_filter.pdf",width = 6,height = 4)

ggplot(data=filter.d,aes(x=preprocess,y=detected)) +
  geom_violin()  +
  facet_grid(.~kept) +
  ggbeeswarm::geom_quasirandom(size = 0.05, aes(colour = kept)) +
  coord_flip()+scale_color_jco()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
        text=element_text(size=15)) +
  geom_text(data=filter.d %>% group_by(preprocess,kept) %>% summarise(n=n()) ,aes(x=preprocess,y=300,label=n),size=3) +
  #scale_fill_manual(values = mycol)+
  theme_bw() + labs(x="",y="Number of detected genes per cell",fill="preprocess")
dev.off()



panel.background = element_rect(fill = "white"),
panel.grid.major = element_line(size = 0.1, colour = "gainsboro"), 
panel.grid.minor = element_line(size = 0.1,colour = "gainsboro"),axis.line = element_line(colour = "black")
