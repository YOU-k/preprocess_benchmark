library(scran)
library(scater)
library(CellBench)
library(tidyverse)

read.path<- "/stornext/HPCScratch/home/you.y/preprocess_update/SCEs/norm/"
save.path <- "/stornext/HPCScratch/home/you.y/preprocess_update/SCEs/norm/hvg/"
scran_high_var = function(sce,topn=1500){
  if (nrow(sce)>1500){
    dec <- modelGeneVar(sce)
    # Ordering by most interesting genes for inspection.
    chosen <- getTopHVGs(dec, n=topn)
    sce <- sce[chosen,]
  }
  return(sce)
}

get_hvg_curve = function(sce){
  if (nrow(sce)>1500){
    dec <- modelGeneVar(sce) 
  }else{dec<- NA}
  return(dec)
}






hivar_method = list(scran_hi = scran_high_var,
                    curve=get_hvg_curve)
droplet_recode <- function(vec){
  return(recode(vec, "scpipe"="scPipe","zumis"="zUMIs","kb"="kallisto bustools",
                "optimus" ="Optimus","drop"="dropSeqPipe","cellranger"="Cell Ranger",
                "sa"="salmon_SA","splici"="salmon_splici"))
}

dir(read.path)

designs=c("sc3cl","sc5cl","sc5clv3","mus1","mus2","pbmc5k","pbmc10k")
designs=c("sc5clv3")
for (d in designs) {
  tmp <- readRDS(paste0(read.path,d,"_norm.rds"))
  tmp$design <- recode(tmp$design,"mus1"="10xv2_lungtissue1","mus2"="10xv2_lungtissue2","sc3cl"="10xv2_3cellline",
                       "sc5cl"="10xv2_5cellline","sc5clv3"="10xv3_5cellline","pbmc5k"="10xv3_pbmc5k",
                       "pbmc10k"="10xv3_pbmc10k")
  tmp$data <- droplet_recode(tmp$data)
  tmp %>% arrange(design,data,norm_method) %>% apply_methods(hivar_method) -> tmp_hvged
  
  tmp_hvged %>% spread(hivar_method,result) -> tmp_all
  saveRDS(tmp_all %>% dplyr::select(-curve),paste0(save.path,d,"_hvg.rds"))
  
  tmp_all %>% 
    #dplyr::select(-scran_hi) %>% 
    dplyr::filter(!is.na(curve)) -> tmp_curve

  for (n in unique(tmp_curve$norm_method)) {
    tmp_curve_sub <- tmp_curve[tmp_curve$norm_method==n,]
    pdf(paste0("/stornext/HPCScratch/home/you.y/preprocess_update/results/droplet-based/hvg/",n,"_",d,".pdf"),width = 8,height =8)
    par(mfrow=c(3, 3))
    for (i in 1:nrow(tmp_curve_sub)){
      decX <- tmp_curve_sub$curve[[i]]
      plot(decX$mean, decX$total, xlab="Mean log-expression", 
      ylab="Variance", main=tmp_curve_sub$data[i])
      curve(metadata(decX)$trend(x), col="red", add=TRUE) }
    dev.off()
  }

  saveRDS(tmp_curve,paste0(save.path,d,"_curve.rds"))
  
}


