#!/bin/bash

for k in $(ls /home/$USER/ICA1/Mapping/*.sam); 
do
	samtools sort -o ${j}Sorted.bam ${j};
done
ls

for k in *.bam;
do
	basename $k .bam;
done | sed ':a;N;$!ba;s/\n/ /g' # Get rid of the annoying .samSorted.bam file names

ls *.bam | parallel samtools index '{}'

for l in *.bai;
do
        basename $l .bai;
done | sed ':a;N;$!ba;s/\n/ /g' # Get rid of the annoying .samSorted.bam.bai file names

source /home/$USER/bedtools.sh
