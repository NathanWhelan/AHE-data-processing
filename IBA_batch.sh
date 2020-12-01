#!/bin/bash




####Most of the renaming is very specific to your directory tree and file names. Please do not mindlessly follow this script.

for FILE in ./trimmed_F_Breinholt_et_al_2/*R1*.fq
do
NAME=`echo $FILE |awk -F "R1" '{print $1}'` 
NAME2=`echo $FILE |awk -F "R1" '{print $1}'|awk -F "_" '{print $10 "_" $11}'`

##I used ELithasia to get more loci for -label
python IBA.py -raw1 $FILE -raw2 $NAME\R2\_combo_val_2.fq -d ./FWS_123701_REF -n 6 -t 1 -p 30 -g 200 -c 10 -taxa $NAME2 -label Biomphalaria -K 25
done
