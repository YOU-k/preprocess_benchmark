#!/bin/bash
#PBS -l nodes=1:ppn=8,mem=60gb,walltime=100:00:00
#PBS -m abe
#PBS -M you.y@wehi.edu.au

module load cellranger
work_p="/stornext/General/data/user_managed/grpu_mritchie_1/Yue/preprocess/cellranger/sc_3cl"

cd $work_p
START=$(date +%s)
cellranger count --id=Undetermined --sample Undetermined --transcriptome=/stornext/General/data/user_managed/grpu_mritchie_1/SCmixologyV3/refdata-cellranger-GRCh38-3.0.0 --fastqs=/stornext/General/data/user_managed/grpu_mritchie_1/SCmixology/NN76/fastq --localcores=16 --localmem=16 --expect-cells=5000
END=$(date +%s)
DIFF=$(( $END - $START ))
echo "It took $DIFF seconds for cellranger"
