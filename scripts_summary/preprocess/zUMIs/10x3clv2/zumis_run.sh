#!/bin/bash
#PBS -l nodes=1:ppn=8,mem=100gb,walltime=100:00:00
#PBS -m abe
#PBS -M you.y@wehi.edu.au

tool_p="/stornext/General/data/user_managed/grpu_mritchie_1/Yue/preprocess/tools/zUMIs"
work_p="/stornext/General/data/user_managed/grpu_mritchie_1/Yue/preprocess/zumis/sc_3cl"
cd $work_p
module load samtools
module load R/3.6.1
module load STAR
module load pigz
START=$(date +%s)
$tool_p/zUMIs-master.sh -y test.yaml
END=$(date +%s)
DIFF=$(( $END - $START ))
echo "It took $DIFF seconds for zumis"
