#!/bin/bash


## Set file paths
filePath=/exports/molepi/users/ljsinke/GoDMC/17/*
logdir=/exports/molepi/users/ljsinke/LLS/IL6_Data/mQTLs/logs/


## Loop through all files in directory
for file in $filePath
do
    base="$(basename $file .txt.gz)"

    ## Run CpG script for each file
    sbatch \
    --job-name=$base \
    --time=12:0:0 \
    --cpus-per-task=1 \
    --mem=20G \
    --partition=all,highmem \
    --output="$logdir/$base.log" \
    --wrap "
        bash /exports/molepi/users/ljsinke/LLS/IL6_Data/mQTLs/scripts/17_mQTL-Extract.sh \
            -f $file \
            -o /exports/molepi/users/ljsinke/LLS/IL6_Data/mQTLs/out/ \
            -b $base \
            -c /exports/molepi/users/ljsinke/LLS/IL6_Data/mQTLs/IL6_cgList.txt
    "
done
