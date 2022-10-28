#!/bin/bash

# Import the T. congolense bedfile
cp -r /localdisk/data/BPSM/ICA1/TriTrypDB-46_TcongolenseIL3000_2019.bed /home/$USER/ICA1/Mapping
mv TriTrypDB-46_TcongolenseIL3000_2019.bed Tco-genome3000.bed # As usual, we don't like long, complicated names
ls # The bedfile should have been renamed by now

# Convert bam files to bed first - this will ease looking at intersect
for o in $(ls /home/$USER/ICA1/Mapping/*.bam);
do 
	bedtools bamtobed -i ${o} > ${o}.bed
done
ls

# Finds all unique overlaps in the genomes with the read sequences
for p in $(ls /home/$USER/ICA1/Mapping/*.bed);
	bedtools intersect -wa -a Tco-genome3000.bed -b ${p} | sort | uniq > ${p}IntersectOutput.txt
done
ls

# REFERENCES
# https://bedtools.readthedocs.io/en/latest/content/tools/map.html 
# https://www.protocols.io/view/week-6-mapping-with-bowtie2-kqdg3qyev25z/v2?step=12 
# https://sciberg.com/resources/bioinformatics-scripts/manipulation-with-fasta-and-fastq-files-in-linux 
# https://onestopdataanalysis.com/fastq-to-fasta/ 
# https://wikis.utexas.edu/display/CoreNGSTools/Analysis+using+BEDTools
# https://www.biostars.org/p/213700/ 
# https://www.biostars.org/p/179234/
# https://www.biostars.org/p/467896/
# https://www.biostars.org/p/123383/
# https://stackoverflow.com/questions/59777276/bash-script-loop-through-files-and-add-to-different-variables
# https://biohpc.cornell.edu/lab/doc/linux_examples_exercises_v5.html

