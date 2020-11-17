#!/bin/bash
#PBS -l nodes=1:ppn=8,mem=20gb,walltime=100:00:00
#PBS -m abe
#PBS -M you.y@wehi.edu.au
bustools="/stornext/General/data/user_managed/grpu_mritchie_1/Yue/preprocess/tools/bustools/build/src/bustools"
work_p="/stornext/General/data/user_managed/grpu_mritchie_1/Yue/preprocess/kallisto_bus"
data_p="/stornext/General/data/user_managed/grpu_mritchie_1/Yue/preprocess/data"
seq_p="/stornext/General/data/user_managed/grpu_mritchie_1/Yue/preprocess/dropseqpipe/sc_5cl_v2/RAW_DATA"
cd $work_p
module load anaconda3
source activate conda-env
START=$(date +%s)
#index has been built before
#pseudoalign the reads
kallisto bus -i $work_p/homo_ercc.idx -o sc_5cl_v2/ -x 10xv2 -t 4 $seq_p/Lib90_combined_S1_L001_R1.fastq.gz $seq_p/Lib90_combined_S1_L001_R2.fastq.gz

#run bustools
cd sc_5cl_v2/
mkdir genecount/ tmp/
$bustools correct -w $data_p/10xv2_whitelist.txt -p output.bus | $bustools sort -t 4 -T ./tmp -p - | $bustools count -o genecount/genes -g $data_p/test.tsv -e matrix.ec -t transcripts.txt --genecounts - 

END=$(date +%s)
DIFF=$(( $END - $START ))
echo "It took $DIFF seconds for kallisto_bustools"
