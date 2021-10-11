#!/bin/bash
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=12
#SBATCH --mem=80G
#SBATCH --mail-user=you.y@wehi.edu.au
#SBATCH --mail-type=BEGIN,END,FAIL


module load anaconda3
source activate salmon

project="sc_5cl"
seq_p="/stornext/General/data/user_managed/grpu_mritchie_1/SCmixology/scbench_5cellline_10x/AGRF_CAGRF17555_CCHA5ANXX/fastqs"
fastq1=$seq_p/Lib90_combined_S1_L001_R1_001.fastq.gz
fastq2=$seq_p/Lib90_combined_S1_L001_R2_001.fastq.gz

work_p="/stornext/General/data/user_managed/grpu_mritchie_1/Yue/preprocess_update/raw_results/SA"

data_p="/stornext/General/data/user_managed/grpu_mritchie_1/Yue/preprocess/data"

cd $work_p/$project

salmon alevin -lISR -1 $fastq1 -2 $fastq2 --chromium -i $work_p/salmon_index_human -p 12 -o output --tgMap $data_p/test12.tsv --keepCBFraction 1