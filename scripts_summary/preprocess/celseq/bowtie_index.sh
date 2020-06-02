#!/bin/bash
#PBS -l nodes=1:ppn=8,mem=30gb,walltime=100:00:00
#PBS -m abe
#PBS -M you.y@wehi.edu.au

work_p="/stornext/General/data/user_managed/grpu_mritchie_1/Yue/preprocess/celseq2"
cd $work_p

module load bowtie2
START=$(date +%s)

bowtie2-build ../data/homo_toplevel_ercc.fa ./index/homo_bowtie_index
END=$(date +%s)
DIFF=$(( $END - $START ))
echo "It took $DIFF seconds for bowtie2 alignment"
