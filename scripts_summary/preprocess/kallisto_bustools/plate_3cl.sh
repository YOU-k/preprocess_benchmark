#!/bin/bash
#PBS -l nodes=1:ppn=8,mem=20gb,walltime=100:00:00
#PBS -m abe
#PBS -M you.y@wehi.edu.au
bustools="/stornext/General/data/user_managed/grpu_mritchie_1/Yue/preprocess/tools/bustools/build/src/bustools"
work_p="/stornext/General/data/user_managed/grpu_mritchie_1/Yue/preprocess/kallisto_bus"
data_p="/stornext/General/data/user_managed/grpu_mritchie_1/Yue/preprocess/data"
seq_p="/stornext/General/data/user_managed/grpu_mritchie_1/SCmixology/NN84/luyiT030118_nn84_NEXS2595"
cd $work_p
module load anaconda3
source activate conda-env
START=$(date +%s)
#index has been built before
#pseudoalign the reads
kallisto bus -i $work_p/homo_ercc.idx -o nn84-4/ -x 0,6,13:0,0,6:1,0,0 -t 4 $seq_p/RPI4-TGACCA_S2_R1.fastq.gz $seq_p/RPI4-TGACCA_S2_R2.fastq.gz

#run bustools
cd nn84-4/
mkdir genecount/ tmp/
$bustools correct -w $data_p/barcode7base.txt -p output.bus | $bustools sort -t 4 -T ./tmp -p - | $bustools count -o genecount/genes -g $data_p/test.tsv -e matrix.ec -t transcripts.txt --genecounts - 

END=$(date +%s)
DIFF=$(( $END - $START ))
echo "It took $DIFF seconds for kallisto_bustools"
