library(scPipe)
library(SingleCellExperiment)
library(Rsubread)
project="pbmc5k"
data_dir = paste0("/stornext/HPCScratch/home/you.y/preprocess_update/raw_results/scpipe/",project)
fa_fn  = "/stornext/General/data/user_managed/grpu_mritchie_1/Yue/preprocess/data/genome.fa"
gff3_fn="/stornext/General/data/user_managed/grpu_mritchie_1/Yue/preprocess/data/Homo_sapiens.GRCh38.98.gff3"
index_dir="/stornext/HPCScratch/home/you.y/preprocess_update/raw_results/scpipe/"
##file path 
fq_R2="/stornext/General/data/user_managed/grpu_mritchie_1/Yue/preprocess_update/pbmc_data/raw_data/5k_pbmc/com5kpbmc_S1_L001_R1.fastq.gz"
fq_R1="/stornext/General/data/user_managed/grpu_mritchie_1/Yue/preprocess_update/pbmc_data/raw_data/5k_pbmc/com5kpbmc_S1_L001_R2.fastq.gz"
print("###   trim barcode")
read_structure <- list(
  bs1 = -1,   # barcode start position in fq_R1, -1 indicates no barcode
  bl1 = 0,    # barcode length in fq_R1, 0 since no barcode present
  bs2 = 0,    # barcode start position in fq_R2
  bl2 = 16,   # barcode length in fq_R2
  us = 16,    # UMI start position in fq_R2
  ul = 12     # UMI length
)
sc_trim_barcode(file.path(data_dir,"combined.fastq.gz"),
                fq_R1,
                fq_R2,
                read_structure = read_structure)
print("###   alignment")
align(
  index = file.path(index_dir, "human_index"),
  readfile1 = file.path(data_dir,"combined.fastq.gz"),
  output_file = file.path(data_dir,"out_aln.bam"),
  nthreads = 16
)


print("###   detect barcode")
sc_detect_bc(
  infq = file.path(data_dir,"combined.fastq.gz"),
  outcsv = file.path(data_dir, "barcode_anno.csv"),       # bacode annotation output file name
  bc_len = read_structure$bl2, # 16bp barcode
  max_reads = 50000000, 
  white_list_file="/stornext/General/data/user_managed/grpu_mritchie_1/Yue/preprocess/data/10xv3_whitelist.txt"
)

print("###   map to exon and count")
sc_count_aligned_bam(
  inbam = file.path(data_dir,"out_aln.bam"),
  outbam = file.path(data_dir,"out_map.bam"),
  annofn = gff3_fn,
  bc_len = read_structure$bl2,
  UMI_len = read_structure$ul,
  outdir = data_dir,
  bc_anno = file.path(data_dir, "barcode_anno.csv")
)

