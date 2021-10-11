#!/bin/bash
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=100G
#SBATCH --mail-user=you.y@wehi.edu.au
#SBATCH --mail-type=BEGIN,END,FAIL


module load cellranger/6.0.0
project="pbmc10k"
work_p="/stornext/HPCScratch/home/you.y/preprocess_update/raw_results/cellranger/cellranger"
fastq_p="/stornext/General/data/user_managed/grpu_mritchie_1/Yue/preprocess_update/pbmc_data/raw_data/10k_pbmc/pbmc_10k_v3_fastqs"
data_p="/stornext/General/data/user_managed/grpu_mritchie_1/Yue/preprocess/data"

cd $work_p
cd $project

cellranger count --id=pbmc_10k_v3 --sample pbmc_10k_v3 \
--transcriptome=$data_p/hg_update --fastqs=$fastq_p \
--localcores=8 --localmem=100 --expect-cells=10000
