#!/bin/bash

# Check where we are right now
pwd

# Import raw T. congolense genome sequence (reference data) from the BPSM directory to the ICA1 directory
chmod -R 777 ICA1 # change permission first, otherwise there would be an issue later when we try to import files into this directory
cp -r /localdisk/data/BPSM/ICA1/Tcongo_genome/ /home/$USER/ICA1/Mapping
cd Tcongo_genome
ls

# Let's unzip the gzip compressed fasta file for our reference data
# gunzip -c TriTrypDB-46_TcongolenseIL3000_2019_Genome.fasta.gz | awk '{if(NR%4==1) {printf(">%s\n",substr($0,2));} else if(NR%4==2) print;}' > /home/$USER/ICA1/Mapping/Tco-genome2.fasta # Convert it into fasta and store it to the previously created Mapping directory
cp TriTrypDB-46_TcongolenseIL3000_2019_Genome.fasta.gz /home/$USER/ICA1/Mapping # Copy this fasta.gz file into the previously created Mapping folder 

# Let's look into our Mapping directory now
cd ..
ls
# Let's rename this complicated file name
mv TriTrypDB-46_TcongolenseIL3000_2019_Genome.fasta.gz Tco-genome.fasta.gz # File name is too long. Let's shorten it, shall we?
ls

# Prepare an index of our reference data so that we can map our reads to it
bowtie2-build Tco-genome.fasta.gz Tco-genome.btindex
ls

# Well that works! Let's perform a run-dry
# for i in $(ls /home/$USER/ICA1/Mapping/*1.fq.gz); 
# do 
# 	echo "bowtie2 -p 4 -x Tco-genome.btindex -1 $i -2 ${i/_1.fq.gz/_2.fq.gz} -S ${i/_1.fq.gz/_2.fq.gz}Aligned.sam";
# done 
# This will show the paths to our pair sequence files, the corresponding SAM files
# Hooray, that works! Let's do this
for i in $(ls /home/$USER/ICA1/Mapping/*1.fq.gz); 
do 
 	bowtie2 --threads 48 -x Tco-genome.btindex -1 $i -2 ${i/_1.fq.gz/_2.fq.gz} -S ${i/_1.fq.gz}Aligned.sam;
done
ls

for j in *.sam;
do 
	basename $j.sam;
done | sed ':a;N;$!ba;s/\n/ /g'

rm -f *.fq.gz # Remove all fq gz-compressed files

#source /home/$USER/samtools.sh
