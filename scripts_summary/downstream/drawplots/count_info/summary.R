#get qc metric
get_qc <- function(sce){
  return(df <- data.frame(total_max= log10(max(sce$total)),
                          total_min=log10(min(sce$total)),
                          detect_max= log10(max(sce$detected)),
                          detect_min=log10(min(sce$detected)),
                          mt_max=sce$subsets_Mito_percent[order(sce$subsets_Mito_percent,decreasing = TRUE)[1]] ,
                          mt_min=sce$subsets_Mito_percent[order(sce$subsets_Mito_percent)[1]] ))
}
plate_path <- "/stornext/HPCScratch/home/you.y/preprocess_update/SCEs/raw/plate_based/"
dir(plate_path)
design=c("rnamix","plate3cl","plate5cl1","plate5cl2","plate5cl3")
dataname=c("scpipe","zumis","celseq","scruff","kb")
paste(apply(expand.grid(design, dataname), 1, paste, collapse="_"), ".rds",sep="") ->files
files %in% dir(plate_path)

result0 <- lapply(files, function(file){readRDS(file.path(plate_path,file))})
remove_doublet <- function(sce){
  if("logcounts" %in% assayNames(sce)){
  as.character(sce$demuxlet_cls) ->sce$demuxlet_cls
  sce$demuxlet_cls[is.na(sce$demuxlet_cls)] <- "no"
  sce <- sce[,!sce$demuxlet_cls=="DBL"]}
  return(sce)
}
lapply(result0, remove_doublet) -> result0
lapply(result0, get_qc) -> qc_plate
tibble(design=rep(design,5),method=rep(dataname,each=5),qc=qc_plate) %>% unnest_wider(qc) -> tmp1

write.csv(tmp1,
          file = "/home/users/allstaff/you.y/plate_qc.csv",quote = FALSE)







drop_path <- "/stornext/HPCScratch/home/you.y/preprocess_update/SCEs/raw/droplet_based/"
dir(drop_path)
design=c("mus1","mus2","pbmc5k","pbmc10k","sc3cl","sc5cl","sc5clv3")
dataname=c("kb","zumis","scpipe","sa","splici","drop","cellranger")
paste(apply(expand.grid(design, dataname), 1, paste, collapse="_"), ".rds",sep="") ->files

files %in% dir(drop_path)

result0 <- lapply(files, function(file){readRDS(file.path(drop_path,file))})

lapply(result0, get_qc) -> qc_drop
tibble(design=rep(design,7),method=rep(dataname,each=7),qc=qc_drop) %>% unnest_wider(qc) -> tmp2

bind_rows(tmp1,tmp2) -> tmp




#optimus
design=c("mus1","mus2")
dataname=c("optimus")
paste(apply(expand.grid(design, dataname), 1, paste, collapse="_"), ".rds",sep="") ->files

files %in% dir(drop_path)

result0 <- lapply(files, function(file){readRDS(file.path(drop_path,file))})

lapply(result0, get_qc) -> qc_drop2
tibble(design=rep(design,1),method=rep(dataname,each=1),qc=qc_drop2) %>% unnest_wider(qc) -> tmp3

bind_rows(tmp1,tmp2,tmp3) -> tmp

tmp$design <- recode(tmp$design,"mus1"="10xv2_lungtissue1","mus2"="10xv2_lungtissue2","sc3cl"="10xv2_3cellline",
                     "sc5cl"="10xv2_5cellline","sc5clv3"="10xv3_5cellline","pbmc5k"="10xv3_pbmc5k",
                     "pbmc10k"="10xv3_pbmc10k","plate3cl"="plate_3cellline","plate5cl1"="plate_5cellline1",
                     "plate5cl2"="plate_5cellline2","plate5cl3"="plate_5cellline3","rnamix"="plate_RNAmix")
tmp$method <- recode(tmp$method, "scpipe"="scPipe","zumis"="zUMIs","kb"="kallisto bustools",
                                "optimus" ="Optimus","drop"="dropSeqPipe","cellranger"="Cell Ranger",
                                "sa"="alevin","splici"="alevin-fry","celseq"="celseq2")
write.csv(tmp %>% arrange(design,method),
          file = "/home/users/allstaff/you.y/table_s3.csv",quote = FALSE)
