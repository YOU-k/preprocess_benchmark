#!/bin/bash
#PBS -l nodes=1:ppn=8,mem=120gb,walltime=200:00:00
#PBS -m abe
#PBS -M you.y@wehi.edu.au

tool_p="/stornext/General/data/user_managed/grpu_mritchie_1/Yue/preprocess/tools/zUMIs"
work_p="/stornext/General/data/user_managed/grpu_mritchie_1/Yue/preprocess/zumis2/sc_5cl_v3"
cd $work_p
module load samtools
module load R/3.6.2
module load STAR/2.6.1c
module load pigz
START=$(date +%s)
$tool_p/zUMIs-master.sh -y test.yaml
END=$(date +%s)
DIFF=$(( $END - $START ))
echo "It took $DIFF seconds for zumis"
