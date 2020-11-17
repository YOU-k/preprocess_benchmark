#!/bin/bash
#PBS -l nodes=1:ppn=8,mem=50gb,walltime=100:00:00
#PBS -m abe
#PBS -M you.y@wehi.edu.au

salmon="/stornext/General/data/user_managed/grpu_mritchie_1/Yue/preprocess/tools/salmon-latest_linux_x86_64/bin/salmon"
work_p="/stornext/General/data/user_managed/grpu_mritchie_1/Yue/preprocess/sal_ale/mouse"
data_p="/stornext/General/data/user_managed/grpu_mritchie_1/Yue/preprocess/data"
seq_p="/stornext/General/data/user_managed/grpu_mritchie_1/Yue/preprocess/mouse/data/result2/combined"
ref_p="/stornext/General/data/user_managed/grpu_mritchie_1/Yue/preprocess/mouse/reference"
cd $work_p
START=$(date +%s)
#build an index
$salmon alevin --dumpMtx -l ISR -1 $seq_p/bamtofastq_S1_L002_R1.fastq.gz -2 $seq_p/bamtofastq_S1_L002_R2.fastq.gz --chromium -i $work_p/mus2.idx -p 10 -o $work_p/mus2/alevin_output --tgMap $ref_p/txg12.tsv --keepCBFraction 1
END=$(date +%s)
DIFF=$(( $END - $START ))
echo "It took $DIFF seconds for salmon_alevin"