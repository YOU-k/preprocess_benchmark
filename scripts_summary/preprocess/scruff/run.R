library(scruff)
# NOTE: Requires Rsubread index and TxDb objects for the reference genome.
sc1_s1_r1<- "/stornext/General/data/user_managed/grpu_mritchie_1/SCmixology/NN79/daniela_nn79_NEXS2512_171117/SC1_S1_R1.fastq.gz"
sc1_s1_r2<- "/stornext/General/data/user_managed/grpu_mritchie_1/SCmixology/NN79/daniela_nn79_NEXS2512_171117/SC1_S1_R2.fastq.gz"
fasta <- "/stornext/General/data/user_managed/grpu_mritchie_1/Yue/preprocess/data/homo_ercc.fa"
gtf <- "/stornext/General/data/user_managed/grpu_mritchie_1/Yue/preprocess/data/homo_ercc.gtf"

barcode <- read.table("/stornext/General/data/user_managed/grpu_mritchie_1/Yue/preprocess/data/barcode_annotation")
barcode <- as.character(barcode[,1])
de <- demultiplex(project = "sc1_s1",
experiment = "sc1",
lane = "l1",
read1Path = sc1_s1_r1,
read2Path = sc1_s1_r2,
barcode,
bcStart = 7,
bcStop = 14,
bcEdit = 1,
umiStart = 1,
umiStop = 6,
keep = 75,
minQual = 10,
yieldReads = 1e+06,
verbose = TRUE,
overwrite = TRUE,
cores = 2)

if (!requireNamespace("Rsubread", quietly = TRUE)) {
message("Package \"Rsubread\" needed.",
" Please install it if you are using Linux or macOS systems.",
" The function is not available in Windows environment.\n")
} else {
# For details, please refer to Rsubread user manual.
# Specify the basename for Rsubread index
indexBase <- "homo_ercc_index"
Rsubread::buildindex(basename = indexBase,
reference = fasta,
indexSplit = FALSE)
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
cores = 2)
}


if (requireNamespace("Rsubread", quietly = TRUE)) {
sce = countUMI(al,
gtf,
umiEdit = 1,
format = "BAM",
cellPerWell = c(0,0,rep(1,20),0,0,rep(1,335),0,0,rep(1,20),0,0),
verbose = TRUE,
cores = 2)
}

