#!/bin/bash

# This should sort our sam files accordinly before they can all be used to map against our bedfile with bedtools
for k in $(ls /home/$USER/ICA1/Mapping/*.sam); 
do
 	samtools sort ${k} -o ${k}Sorted.bam;
done

# This should get rid of all the annoying .fq.gzSorted.bam names
for l in *.bam;
do
        basename $l.bam;
done | sed ':a;N;$!ba;s/\n/ /g'
ls

rm -f *.fq.gz # Remove all fq gz-compressed files
ls

# Index our bam files
ls *.bam | parallel samtools index '{}'

while true; do
    read -p "We will now proceed to mapping our read sequences with Bedtools. Please press Y to proceed or N to cancel." yn
    case $yn in
        [Yy]* ) source /home/$USER/bedtools.sh; break;;
        [Nn]* ) exit;;
        * ) echo "Please answer Y or N.";;
    esac
done

