#!/bin/bash
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=150G
#SBATCH --mail-user=you.y@wehi.edu.au
#SBATCH --mail-type=BEGIN,END,FAIL

##mouse index

module load singularity
singular_p="/stornext/General/data/user_managed/grpu_mritchie_1/Yue/preprocess_update/tools/salmon"

data_p="/stornext/General/data/user_managed/grpu_mritchie_1/Yue/preprocess/mouse/reference"
work_p="/stornext/General/data/user_managed/grpu_mritchie_1/Yue/preprocess_update/raw_results/splici"


cp $data_p/Mus_musculus.GRCm38.dna.primary_assembly.fa $singular_p/Mus_musculus.GRCm38.dna.primary_assembly.fa

cp $data_p/cellranger/Mus_musculus.GRCm38.99.gtf $singular_p/Mus_musculus.GRCm38.99.gtf
cd $singular_p
singularity exec --cleanenv \
--bind $singular_p:/workdir \
--pwd /usefulaf/bash usefulaf_latest.sif \
./simpleaf index \
-f /workdir/Mus_musculus.GRCm38.dna.primary_assembly.fa \
-g /workdir/Mus_musculus.GRCm38.99.gtf \
-l 91 -t 16 -o /workdir/mus_splici


rm $singular_p/Mus_musculus.GRCm38.dna.primary_assembly.fa
rm $singular_p/Mus_musculus.GRCm38.99.gtf