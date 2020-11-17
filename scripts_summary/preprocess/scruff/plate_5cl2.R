library(scruff)
# NOTE: Requires Rsubread index and TxDb objects for the reference genome.
n110_r1<- "/stornext/General/data/user_managed/grpu_mritchie_1/SCmixology/NN110/daniela_NN110_060818/5_S3_R1.fastq.gz"
n110_r2<- "/stornext/General/data/user_managed/grpu_mritchie_1/SCmixology/NN110/daniela_NN110_060818/5_S3_R2.fastq.gz"
fasta <- "/stornext/General/data/user_managed/grpu_mritchie_1/Yue/preprocess/data/homo_ercc.fa"
gtf <- "/stornext/General/data/user_managed/grpu_mritchie_1/Yue/preprocess/data/homo_ercc.gtf"

barcode <- read.table("/stornext/General/data/user_managed/grpu_mritchie_1/Yue/preprocess/data/barcode7base.txt")
barcode <- as.character(barcode[,1])
de <- demultiplex(project = "n110-5",
experiment = "n110-5",
lane = "l1",
read1Path = n110_r1,
read2Path = n110_r2,
barcode,
bcStart = 7,
bcStop = 13,
bcEdit = 1,
umiStart = 1,
umiStop = 6,
keep = 75,
minQual = 10,
yieldReads = 1e+06,
verbose = TRUE,
overwrite = TRUE,
cores = 8)

if (!requireNamespace("Rsubread", quietly = TRUE)) {
  message("Package \"Rsubread\" needed.",
          " Please install it if you are using Linux or macOS systems.",
          " The function is not available in Windows environment.\n")
} else {
  # For details, please refer to Rsubread user manual.
  # Specify the basename for Rsubread index
  index_dir="/stornext/General/data/user_managed/grpu_mritchie_1/Yue/preprocess/computation/scruff/sub8"
  indexBase <- file.path(index_dir, "homo_ercc_index")
}


# Align the reads using Rsubread
if (requireNamespace("Rsubread", quietly = TRUE)) {
al <- alignRsubread(de,
indexBase,
unique = FALSE,
nBestLocations = 1,
format = "BAM",
overwrite = TRUE,
verbose = TRUE,
cores = 8)
}


if (requireNamespace("Rsubread", quietly = TRUE)) {
sce = countUMI(al,
gtf,
umiEdit = 1,
format = "BAM",
cellPerWell = rep(1,383),
verbose = TRUE,
cores = 8)
}

