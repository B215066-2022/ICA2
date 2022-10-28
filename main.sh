#!/bin/bash

pwd
mkdir ICA1
cd ICA1

git init # Initialising git

cp -r /localdisk/data/BPSM/ICA1/fastq/ /home/$USER/ICA1
ls
cd fastq
ls

# looks like there are soxme gzip files. Let's perform FATAqc all the gzip compressed fastq files
parallel fastqc :::  *.fq.gz
ls	

# Let's make a new directory ICA1/Mapping so that all our read sequences can be grouped together to map with the reference genome
cd ..
mkdir Mapping
chmod -R 777 Mapping

# unzip -o *.zip
# gunzip -c *.fq.gz | awk '{if(NR%4==1) {printf(">%s\n",substr($0,2));} else if(NR%4==2) print;}' > /home/$USER/ICA1/Mapping/*.fasta

# mv *.zip /home/$USER/ICA1/Mapping
# mv *.html /home/$USER/ICA1/Mapping

# Let's move our read sequences into the same directory in ICA/Mapping, one pair at a time
for filefq in "ls /home/$USER/ICA1/fastq/*.fq.gz"
do
	dir="/home/$USER/ICA1/fastq"
	base1=$(basename $file "*_1.fq.gz")
	base2=$(basename $file "*_2.fq.gz")
	cp ${dir}/${base1} /home/$USER/ICA1/Mapping
	cp ${dir}/${base2} /home/$USER/ICA1/Mapping
done

ls # Check if the files are here

cd fastq # We forgot to copy one more file Tco.fqfiles
mv Tco.fqfiles /home/$USER/ICA1/Mapping
cd ..
cd Mapping
ls

source bowtie2.sh
