library(scPipe)
library(SingleCellExperiment)
library(Rsubread)
data_dir = "/stornext/General/data/user_managed/grpu_mritchie_1/Yue/preprocess/scpipe/cellmix1"
fa_fn  = "/stornext/General/data/user_managed/grpu_mritchie_1/Yue/preprocess/data/homo_ercc.fa"
gff3_fn="/stornext/General/data/user_managed/grpu_mritchie_1/Yue/preprocess/data/Homo_sapiens.GRCh38_ercc.gff3"
index_dir="/stornext/General/data/user_managed/grpu_mritchie_1/Yue/preprocess/computation/scruff/sub8"
##file path 
fq_R2="/stornext/General/data/user_managed/grpu_mritchie_1/SCmixology/NN79/daniela_nn79_NEXS2512_171117/SC1_S1_R1.fastq.gz"
fq_R1="/stornext/General/data/user_managed/grpu_mritchie_1/SCmixology/NN79/daniela_nn79_NEXS2512_171117/SC1_S1_R2.fastq.gz"
print("###   trim barcode")
read_structure <- list(
  bs1 = -1,   # barcode start position in fq_R1, -1 indicates no barcode
  bl1 = 0,    # barcode length in fq_R1, 0 since no barcode present
  bs2 = 6,    # barcode start position in fq_R2
  bl2 = 8,   # barcode length in fq_R2
  us = 0,    # UMI start position in fq_R2
  ul = 6     # UMI length
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


print("###   map to exon and count")
sc_count_aligned_bam(
  inbam = file.path(data_dir,"out_aln.bam"),
  outbam = file.path(data_dir,"out_map.bam"),
  annofn = gff3_fn,
  bc_len = read_structure$bl2,
  UMI_len = read_structure$ul,
  outdir = data_dir,
  bc_anno = "/stornext/General/data/user_managed/grpu_mritchie_1/Yue/preprocess/data/barcode_annotation.csv"
)

