#!/bin/bash

for j in $(ls /home/$USER/ICA1/Mapping/*.sam); 
do
	samtools sort -o ${j}Sorted.bam ${j};
done
ls

ls *.bam | parallel samtools index '{}'



source bedtools.sh
