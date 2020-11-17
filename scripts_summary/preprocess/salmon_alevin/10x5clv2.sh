#!/bin/bash
#PBS -l nodes=1:ppn=8,mem=50gb,walltime=100:00:00
#PBS -m abe
#PBS -M you.y@wehi.edu.au

salmon="/stornext/General/data/user_managed/grpu_mritchie_1/Yue/preprocess/tools/salmon-latest_linux_x86_64/bin/salmon"
work_p="/stornext/General/data/user_managed/grpu_mritchie_1/Yue/preprocess/sal_ale"
data_p="/stornext/General/data/user_managed/grpu_mritchie_1/Yue/preprocess/data"
seq_p="/stornext/General/data/user_managed/grpu_mritchie_1/Yue/preprocess/dsp/sc_5cl_v2/RAW_DATA"

cd $work_p
START=$(date +%s)
#build an index
$salmon alevin --dumpMtx -l ISR -1 $seq_p/Lib90_combined_S1_L001_R1.fastq.gz -2 $seq_p/Lib90_combined_S1_L001_R2.fastq.gz --chromium -i $work_p/homo_ercc2.idx -p 10 -o $work_p/sc_5cl_v2/alevin_output --tgMap $data_p/test12.tsv --keepCBFraction 1
END=$(date +%s)
DIFF=$(( $END - $START ))
echo "It took $DIFF seconds for salmon_alevin"