#!/bin/bash
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=32
#SBATCH --mem=400G
#SBATCH --mail-user=you.y@wehi.edu.au
#SBATCH --mail-type=BEGIN,END,FAIL


project="mus1"
tool_p="/stornext/General/data/user_managed/grpu_mritchie_1/Yue/preprocess_update/tools/zUMIs"
work_p="/stornext/HPCScratch/home/you.y/preprocess_update/raw_results/zumis"
cd $work_p
cd $project
module load samtools
module load R/4.0.3
module load STAR/2.6.1c
module load pigz
module load anaconda3
source activate myenv


$tool_p/zUMIs.sh -y zUMIs.yaml

