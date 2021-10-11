#!/bin/bash
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=100G
#SBATCH --mail-user=you.y@wehi.edu.au
#SBATCH --mail-type=BEGIN,END,FAIL


module load cellranger/6.0.0
work_p="/stornext/HPCScratch/home/you.y/preprocess_update/raw_results/cellranger/sc_3cl"
data_p="/stornext/General/data/user_managed/grpu_mritchie_1/Yue/preprocess/data"
fastq_p="/stornext/General/data/user_managed/grpu_mritchie_1/SCmixology/NN76/fastq"
cd $data_p

cellranger mkgtf Homo_sapiens.GRCh38.98.gtf Homo_sapiens.GRCh38.98.filtered.gtf \
                 --attribute=gene_biotype:protein_coding \
                 --attribute=gene_biotype:lncRNA \
                 --attribute=gene_biotype:antisense
                 
cellranger mkref --genome=hg_update \
                 --fasta=genome.fa \
                 --genes=Homo_sapiens.GRCh38.98.filtered.gtf \


cd $work_p

cellranger count --id=combined --sample combined \
--transcriptome=$data_p/hg_update --fastqs=$fastq_p \
--localcores=8 --localmem=100 --expect-cells=5000

