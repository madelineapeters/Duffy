#!/bin/bash   
cd ~/PopDataProcessing

line=$1
POP=$2

URL=sed -n ''"$line"'p' < FASTQurls_${POP}list.txt
FILENAME=${URL:(-20)}

                       
wget $URL
scp ~/PopDataProcessing/${FILENAME} peter114@cedar.computecanada.ca:/scratch/peter114/PopDataProcessing/${POP}                                            
rm $FILENAME
