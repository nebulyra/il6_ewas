#!/bin/bash


## Get arguments
while getopts f:o:b:c: flag
do
    case "${flag}" in
        f) file=${OPTARG};;
        o) outputdir=${OPTARG};;
        b) base=${OPTARG};;
        c) cpgList=${OPTARG};;
    esac
done


## Start unzipping file (keeping original)
echo "Unzipping file..."
gunzip -c $file > $outputdir/$base
echo "- File unzipped"


## Loop through all CpGs
echo "Searching for CpGs..."
while read cpg; do
    echo - "Checking $cpg"
    LC_ALL=C fgrep $cpg $outputdir/$base >> $outputdir/$cpg.txt
done < $cpgList


## Remove zip
echo "Removing unzipped file"
rm $outputdir/$base


echo "Done!"
