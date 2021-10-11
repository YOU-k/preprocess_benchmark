#!/bin/bash
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=100G
#SBATCH --mail-user=you.y@wehi.edu.au
#SBATCH --mail-type=BEGIN,END,FAIL

module load cellranger/6.0.0
project="mus2"
work_p="/stornext/General/data/user_managed/grpu_mritchie_1/Yue/preprocess_update/raw_results/cellranger"
fastq_p="/stornext/General/data/user_managed/grpu_mritchie_1/Yue/preprocess/mouse/data/result2"
mouse_ref="/stornext/General/data/user_managed/grpu_mritchie_1/Yue/preprocess/mouse/reference/cellranger/mm10"
cd $work_p
cd $project

cellranger count --id=bamtofastq --sample bamtofastq \
--transcriptome=$mouse_ref --fastqs=$fastq_p \
--localcores=8 --localmem=100 --expect-cells=5000



