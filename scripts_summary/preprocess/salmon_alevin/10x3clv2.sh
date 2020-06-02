#!/bin/bash
#PBS -l nodes=1:ppn=8,mem=50gb,walltime=100:00:00
#PBS -m abe
#PBS -M you.y@wehi.edu.au

salmon="/stornext/General/data/user_managed/grpu_mritchie_1/Yue/preprocess/tools/salmon-latest_linux_x86_64/bin/salmon"
work_p="/stornext/General/data/user_managed/grpu_mritchie_1/Yue/preprocess/sal_ale"
data_p="/stornext/General/data/user_managed/grpu_mritchie_1/Yue/preprocess/data"
seq_p="/stornext/General/data/user_managed/grpu_mritchie_1/SCmixology/NN76/fastq"

cd $work_p
START=$(date +%s)
#build an index
$salmon alevin -l ISR -1 $seq_p/combined_S0_L001_R1_001.fastq.gz -2 $seq_p/combined_S0_L001_R2_001.fastq.gz --chromium -i $work_p/homo_ercc2.idx -p 10 -o $work_p/sc_3cl/alevin_output --tgMap $data_p/test12.tsv --keepCBFraction 1
END=$(date +%s)
DIFF=$(( $END - $START ))
echo "It took $DIFF seconds for salmon_alevin"