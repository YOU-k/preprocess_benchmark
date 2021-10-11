#!/bin/bash
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=60G
#SBATCH --mail-user=you.y@wehi.edu.au
#SBATCH --mail-type=BEGIN,END,FAIL

project="mus1"

module load singularity
singular_p="/stornext/General/data/user_managed/grpu_mritchie_1/Yue/preprocess_update/tools/salmon"

seq_p="/stornext/General/data/user_managed/grpu_mritchie_1/Yue/preprocess/mouse/data/result1/combined"
fastq1="bamtofastq_S1_L002_R1.fastq.gz"
fastq2="bamtofastq_S1_L002_R2.fastq.gz"

cp $seq_p/$fastq1 $singular_p/$fastq1

cp $seq_p/$fastq2 $singular_p/$fastq2

cd $singular_p 
singularity exec --cleanenv \
--bind $singular_p:/workdir \
--pwd /usefulaf/bash usefulaf_latest.sif \
./simpleaf quant \
-1 /workdir/$fastq1 \
-2 /workdir/$fastq2 \
-i /workdir/mus_splici/index \
-o /workdir/quants/$project \
-f u -c v2 -r cr-like \
-m /workdir/mus_splici/ref/transcriptome_splici_fl86_t2g_3col.tsv \
-t 16


rm $singular_p/$fastq1
rm $singular_p/$fastq2
