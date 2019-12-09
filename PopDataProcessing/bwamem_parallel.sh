#!/bin/bash   

ID=$1
line1=$(( 2 * (${ID} - 1) + 1 ))
line2=$(( 2 * (${ID} - 1) + 2 ))
POP=$2

cd /scratch/peter114/PopDataProcessing/${POP}

URL1=`sed -n ''"$line1"'p' FASTQurls_${POP}list.txt`
FILENAMEreturn1=${URL1:(-21)}
FILENAME1=$(echo "$FILENAMEreturn1" | tr -d '\r')

URL2=`sed -n ''"$line2"'p' FASTQurls_${POP}list.txt`
FILENAMEreturn2=${URL2:(-21)}
FILENAME2=$(echo "$FILENAMEreturn2" | tr -d '\r')

bwa mem GRCh38_full_analysis_set_plus_decoy_hla.fa ${FILENAME1} ${FILENAME2} | gzip -3 > aln-se.sam.gz
