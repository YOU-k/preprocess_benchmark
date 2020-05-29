#!/bin/bash
#PBS -l nodes=1:ppn=8,mem=50gb,walltime=100:00:00
#PBS -m abe
#PBS -M you.y@wehi.edu.au

bustools="/stornext/General/data/user_managed/grpu_mritchie_1/Yue/preprocess/tools/bustools/bustools"
work_p="/stornext/General/data/user_managed/grpu_mritchie_1/Yue/preprocess/kallisto_bus"
data_p="/stornext/General/data/user_managed/grpu_mritchie_1/Yue/preprocess/data"

cd $work_p
module load anaconda3
module load kallisto
START=$(date +%s)
#build an index
kallisto index -i homo_ercc.idx -k 31 $data_p/cdna/Homo_sapiens.GRCh38.cdna.all.fa $data_p/Homo_sapiens.GRCh38.ncrna.fa
END=$(date +%s)
DIFF=$(( $END - $START ))
echo "It took $DIFF seconds for kallisto_bus"
