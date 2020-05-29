#!/bin/bash
#PBS -l nodes=1:ppn=8,mem=250gb,walltime=20:00:00
#PBS -m abe

START=$(date +%s)
work_p="/stornext/General/data/user_managed/grpu_mritchie_1/Yue/preprocess/dropseqpipe"
cd $work_p
module load anaconda3
module load java/1.8.0_211
source activate snakemake
snakemake --cores 8 --use-conda --directory sc_3cl

END=$(date +%s)
DIFF=$(( $END - $START ))
echo "It took $DIFF seconds for dropseqpipe"
