#!/bin/bash
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=100G
#SBATCH --mail-user=you.y@wehi.edu.au
#SBATCH --mail-type=BEGIN,END,FAIL


module load cellranger/6.0.0
project="sc_5clv3"
work_p="/stornext/HPCScratch/home/you.y/preprocess_update/raw_results/cellranger/cellranger"
fastq_p="/stornext/General/data/user_managed/grpu_mritchie_1/SCmixologyV3/luyiT_10X_260319/H2YJLBGXB/outs/fastq_path/H2YJLBGXB/Lib90"
data_p="/stornext/General/data/user_managed/grpu_mritchie_1/Yue/preprocess/data"

cd $work_p
cd $project

cellranger count --id=Lib90 --sample Lib90 \
--transcriptome=$data_p/hg_update --fastqs=$fastq_p \
--localcores=8 --localmem=100 --expect-cells=5000
