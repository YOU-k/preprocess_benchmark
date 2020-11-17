#!/bin/bash
#PBS -l nodes=1:ppn=8,mem=60gb,walltime=100:00:00
#PBS -m abe
#PBS -M you.y@wehi.edu.au


work_p="/stornext/General/data/user_managed/grpu_mritchie_1/Yue/preprocess/celseq2/nn84-4"
cd $work_p
module load bowtie2
module load anaconda3
source activate conda-celseq2

START=$(date +%s)
celseq2 --config-file /stornext/General/data/user_managed/grpu_mritchie_1/Yue/preprocess/celseq2/nn84-4/config.yaml \
    --experiment-table /stornext/General/data/user_managed/grpu_mritchie_1/Yue/preprocess/celseq2/nn84-4/experiment_design.txt \
    --output-dir /stornext/General/data/user_managed/grpu_mritchie_1/Yue/preprocess/celseq2/nn84-4/result \
    --j 10

END=$(date +%s)
DIFF=$(( $END - $START ))
echo "It took $DIFF seconds for celseq2"


