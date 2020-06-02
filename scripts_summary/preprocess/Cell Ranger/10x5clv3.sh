#!/bin/bash
#PBS -l nodes=1:ppn=8,mem=60gb,walltime=100:00:00
#PBS -m abe
#PBS -M you.y@wehi.edu.au

module load cellranger
work_p="/stornext/General/data/user_managed/grpu_mritchie_1/Yue/preprocess/cellranger/sc_5cl_v3"

cd $work_p
START=$(date +%s)
cellranger count --id=Lib90_comb --transcriptome=/stornext/General/data/user_managed/grpu_mritchie_1/SCmixologyV3/refdata-cellranger-GRCh38-3.0.0 --fastqs=/stornext/General/data/user_managed/grpu_mritchie_1/SCmixologyV3/luyiT_10X_260319/H2YJLBGXB/outs/fastq_path/H2YJLBGXB/Lib90 --localcores=16 --localmem=16 --expect-cells=5000
END=$(date +%s)
DIFF=$(( $END - $START ))
echo "It took $DIFF seconds for cellranger"
