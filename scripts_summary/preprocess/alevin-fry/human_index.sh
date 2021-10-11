#!/bin/bash
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=60G
#SBATCH --mail-user=you.y@wehi.edu.au
#SBATCH --mail-type=BEGIN,END,FAIL

##human index

module load singularity
singular_p="/stornext/General/data/user_managed/grpu_mritchie_1/Yue/preprocess_update/tools/salmon"

data_p="/stornext/General/data/user_managed/grpu_mritchie_1/Yue/preprocess/data"
work_p="/stornext/General/data/user_managed/grpu_mritchie_1/Yue/preprocess_update/raw_results/splici"


cp $data_p/genome.fa $singular_p/genome.fa

cp $data_p/Homo_sapiens.GRCh38.98.gtf $singular_p/Homo_sapiens.GRCh38.98.gtf

cd $singular_p
singularity exec --cleanenv \
--bind $singular_p:/workdir \
--pwd /usefulaf/bash usefulaf_latest.sif \
./simpleaf index \
-f /workdir/genome.fa \
-g /workdir/Homo_sapiens.GRCh38.98.gtf \
-l 91 -t 16 -o /workdir/human_splici


rm $singular_p/genome.fa
rm $singular_p/Homo_sapiens.GRCh38.98.gtf