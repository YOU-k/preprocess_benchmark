#!/bin/bash
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=50G
#SBATCH --mail-user=you.y@wehi.edu.au
#SBATCH --mail-type=BEGIN,END,FAIL

project="pbmc5k"
work_p="/stornext/General/data/user_managed/grpu_mritchie_1/Yue/preprocess/kallisto_bus"
data_p="/stornext/General/data/user_managed/grpu_mritchie_1/Yue/preprocess/data"
seq_p="/stornext/General/data/user_managed/grpu_mritchie_1/Yue/preprocess_update/pbmc_data/raw_data/5k_pbmc"
dest_p="/stornext/General/data/user_managed/grpu_mritchie_1/Yue/preprocess_update/raw_results/kb"

fastq1=$seq_p/com5kpbmc_S1_L001_R1.fastq.gz
fastq2=$seq_p/com5kpbmc_S1_L001_R2.fastq.gz

module load anaconda3
source activate kb

#index has been built before
#pseudoalign the reads
kallisto bus -i $work_p/homo_ercc.idx -o $dest_p/$project -x 10xv3 -t 4 $fastq1 $fastq2

#run bustools
cd $dest_p/
cd $project
mkdir genecount/ tmp/
bustools correct -w $data_p/10xv3_whitelist.txt -p output.bus | bustools sort -t 4 -T ./tmp -p - | bustools count -o genecount/genes -g $data_p/test.tsv -e matrix.ec -t transcripts.txt --genecounts - 


