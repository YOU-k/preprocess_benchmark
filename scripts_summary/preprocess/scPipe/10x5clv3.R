library(scPipe)
library(SingleCellExperiment)
library(Rsubread)
data_dir = "/stornext/General/data/user_managed/grpu_mritchie_1/Yue/preprocess/scpipe/sc_5cl_v3"
fa_fn  = "/stornext/General/data/user_managed/grpu_mritchie_1/Yue/preprocess/data/homo_ercc.fa"
gff3_fn="/stornext/General/data/user_managed/grpu_mritchie_1/Yue/preprocess/data/Homo_sapiens.GRCh38_ercc.gff3"
index_dir="/stornext/General/data/user_managed/grpu_mritchie_1/Yue/preprocess/computation/scruff/sub8"
##file path 
fq_R2="/stornext/General/data/user_managed/grpu_mritchie_1/SCmixologyV3/luyiT_10X_260319/H2YJLBGXB/outs/fastq_path/H2YJLBGXB/Lib90_combined_S1_L001_R1.fastq.gz"
fq_R1="/stornext/General/data/user_managed/grpu_mritchie_1/SCmixologyV3/luyiT_10X_260319/H2YJLBGXB/outs/fastq_path/H2YJLBGXB/Lib90_combined_S1_L001_R2.fastq.gz"
print("###   trim barcode")
read_structure <- list(
  bs1 = -1,   # barcode start position in fq_R1, -1 indicates no barcode
  bl1 = 0,    # barcode length in fq_R1, 0 since no barcode present
  bs2 = 0,    # barcode start position in fq_R2
  bl2 = 16,   # barcode length in fq_R2
  us = 16,    # UMI start position in fq_R2
  ul = 10     # UMI length
)
sc_trim_barcode(file.path(data_dir,"combined.fastq.gz"),
                fq_R1,
                fq_R2,
                read_structure = read_structure)
print("###   alignment")
align(
  index = file.path(index_dir, "homo_ercc_index"),
  readfile1 = file.path(data_dir,"combined.fastq.gz"),
  output_file = file.path(data_dir,"out_aln.bam"),
  nthreads = 12
)


print("###   detect barcode")
sc_detect_bc(
  infq = file.path(data_dir,"combined.fastq.gz"),
  outcsv = file.path(data_dir, "barcode_anno.csv"),       # bacode annotation output file name
  bc_len = read_structure$bl2, # 16bp barcode
  max_reads = 5000000, 
  number_of_cells = 5000,
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

