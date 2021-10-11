#!/bin/bash
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=12
#SBATCH --mem=150G
#SBATCH --mail-user=you.y@wehi.edu.au
#SBATCH --mail-type=BEGIN,END,FAIL


#human index
module load anaconda3
source activate salmon

data_p="/stornext/General/data/user_managed/grpu_mritchie_1/Yue/preprocess/data"
work_p="/stornext/General/data/user_managed/grpu_mritchie_1/Yue/preprocess_update/raw_results/SA"

cd $work_p
grep "^>" <$data_p/genome.fa | cut -d " " -f 1 > human_decoys.txt
sed -i.bak -e 's/>//g' human_decoys.txt

cat $data_p/cdna/Homo_sapiens.GRCh38.cdna.all.fa $data_p/cdna/Homo_sapiens.GRCh38.ncrna.fa $data_p/genome.fa > $data_p/gentrome.fa

salmon index -t $data_p/gentrome.fa -d human_decoys.txt -p 12 -i salmon_index_human 

