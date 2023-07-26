#!/bin/bash

dir1="/nfs/picower001/lhtsailab/djuna/ABCA7_LOF_BMC_data/10x-4819F/10x-fastqs-L1"
curr_arr=("$(ls -I "Undetermined_*" $dir1)")
curr_arr2="${curr_arr::${#curr_arr}-14}"

for FILE in $curr_arr2
do 
    echo $FILE
    md5sum /nfs/picower001/lhtsailab/djuna/ABCA7_LOF_BMC_data/10x-4819F/10x-fastqs-L1/$FILE >> /nfs/picower001/lhtsailab/djuna/ABCA7_LOF_BMC_data/md5sum_picower001/md5sum_ABCA7_lof_fastqs.txt
done
