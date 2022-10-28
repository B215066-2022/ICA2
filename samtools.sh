#!/bin/bash

# This should sort our sam files accordinly before they can all be used to map against our bedfile with bedtools
for k in $(ls /home/$USER/ICA1/Mapping/*.sam); 
do
	samtools sort -o ${j}Sorted.bam ${j};
done

# Get rid of the annoying .samSorted.bam file names
for k in *.bam;
do
	basename $k .bam;
done | sed ':a;N;$!ba;s/\n/ /g'
ls

# This should index all our bam files
ls *.bam | parallel samtools index '{}'

# Get rid of the annoying .samSorted.bam.bai file names
for l in *.bai;
do
        basename $l .bai;
done | sed ':a;N;$!ba;s/\n/ /g'
ls


source /home/$USER/bedtools.sh
