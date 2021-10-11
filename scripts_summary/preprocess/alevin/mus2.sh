#!/bin/bash
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=12
#SBATCH --mem=80G
#SBATCH --mail-user=you.y@wehi.edu.au
#SBATCH --mail-type=BEGIN,END,FAIL


module load anaconda3
source activate salmon

project="mus2"
ref_p="/stornext/General/data/user_managed/grpu_mritchie_1/Yue/preprocess/mouse/reference"
seq_p="/stornext/General/data/user_managed/grpu_mritchie_1/Yue/preprocess/mouse/data/result2/combined"
fastq1=$seq_p/bamtofastq_S1_L002_R1.fastq.gz
fastq2=$seq_p/bamtofastq_S1_L002_R2.fastq.gz

work_p="/stornext/General/data/user_managed/grpu_mritchie_1/Yue/preprocess_update/raw_results/SA"

data_p="/stornext/General/data/user_managed/grpu_mritchie_1/Yue/preprocess/data"

cd $work_p/$project

salmon alevin -lISR -1 $fastq1 -2 $fastq2 --chromium -i $work_p/salmon_index_mus -p 12 -o output --tgMap $ref_p/txg12.tsv --keepCBFraction 1