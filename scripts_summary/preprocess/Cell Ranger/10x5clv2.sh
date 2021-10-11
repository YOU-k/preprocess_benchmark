#!/bin/bash
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=100G
#SBATCH --mail-user=you.y@wehi.edu.au
#SBATCH --mail-type=BEGIN,END,FAIL


module load cellranger/6.0.0
project="sc_5cl"
work_p="/stornext/HPCScratch/home/you.y/preprocess_update/raw_results/cellranger/cellranger"
fastq_p="/stornext/General/data/user_managed/grpu_mritchie_1/SCmixology/scbench_5cellline_10x/AGRF_CAGRF17555_CCHA5ANXX/fastqs"
data_p="/stornext/General/data/user_managed/grpu_mritchie_1/Yue/preprocess/data"

cd $work_p
cd $project

cellranger count --id=Lib90_combined --sample Lib90_combined \
--transcriptome=$data_p/hg_update --fastqs=$fastq_p \
--localcores=8 --localmem=100 --expect-cells=5000
