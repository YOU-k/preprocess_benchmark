#!/bin/bash
#PBS -l nodes=1:ppn=8,mem=50gb,walltime=100:00:00
#PBS -m abe
#PBS -M you.y@wehi.edu.au

salmon="/stornext/General/data/user_managed/grpu_mritchie_1/Yue/preprocess/tools/salmon-latest_linux_x86_64/bin/salmon"
work_p="/stornext/General/data/user_managed/grpu_mritchie_1/Yue/preprocess/sal_ale/mouse"
data_p="/stornext/General/data/user_managed/grpu_mritchie_1/Yue/preprocess/mouse/reference"

cd $work_p
START=$(date +%s)
#build an index
$salmon index --keepDuplicates -t $data_p/mus_all.fa -i mus2.idx --decoys decoys.txt  -k 31
END=$(date +%s)
DIFF=$(( $END - $START ))
echo "It took $DIFF seconds for salmon"
