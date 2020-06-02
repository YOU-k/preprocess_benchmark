#!/bin/bash
#PBS -l nodes=1:ppn=8,mem=60gb,walltime=100:00:00
#PBS -m abe
#PBS -M you.y@wehi.edu.au

module load cellranger
work_p="/stornext/General/data/user_managed/grpu_mritchie_1/Yue/preprocess/cellranger/mus1"

cd $work_p
START=$(date +%s)
cellranger count --id=bamtofastq --sample bamtofastq --transcriptome=/stornext/General/data/user_managed/grpu_mritchie_1/Yue/preprocess/mouse/reference/cellranger/mm10 --fastqs=/stornext/General/data/user_managed/grpu_mritchie_1/Yue/preprocess/mouse/data/result1 --localcores=16 --localmem=16 --expect-cells=5000 --chemistry=SC3Pv2
END=$(date +%s)
DIFF=$(( $END - $START ))
echo "It took $DIFF seconds for cellranger"
