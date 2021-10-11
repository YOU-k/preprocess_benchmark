#!/bin/bash
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=300G
#SBATCH --mail-user=you.y@wehi.edu.au
#SBATCH --mail-type=BEGIN,END,FAIL


module load cellranger/6.0.0
project="pbmc5k"
work_p="/stornext/HPCScratch/home/you.y/preprocess_update/raw_results/cellranger/cellranger"
fastq_p="/stornext/General/data/user_managed/grpu_mritchie_1/Yue/preprocess_update/pbmc_data/raw_data/5k_pbmc/5k_pbmc_v3_fastqs"
data_p="/stornext/General/data/user_managed/grpu_mritchie_1/Yue/preprocess/data"

cd $work_p
cd $project

cellranger count --id=5k_pbmc_v3 --sample 5k_pbmc_v3 \
--transcriptome=$data_p/hg_update --fastqs=$fastq_p \
--localcores=16 --localmem=300 --expect-cells=5000
