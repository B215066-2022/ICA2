#!/usr/bin/python3

import os
import sys
import shutil
import subprocess # We will be running a lot of background functions that are non python-based
import time # We need this module to gap each print statement

# export NCBI_API_KEY=54b38a82050346352dbaf27e8e02a65f0408" >> ~/.bashrc
# export NCBI_API_KEY=54b38a82050346352dbaf27e8e02a65f0408" >> ~/.bash_profile
# echo ${NCBI_API_KEY}

# Make a working directory first, and let's call it ICA2
os.mkdir("ICA2")

Option1 = 
print("This programme will help you analyse protein sequences of your choice.")
time.sleep(1)
print("First, I need you to specify the name of the organism where your protein may be found e.g mosquito.")
time.sleep(1)
print("If you know the organism taxonomy ID, you can write it in e.g TXID7157.")
time.sleep(1)
QueryA = input('Please input your protein family now: ').upper()
print("Now, I need you to specify the family of protein you are interested in e.g glucose-6-phosphate dehydrogenase.")
time.sleep(1)
QueryB = input('Please input your protein family name now: ').upper()

OptionA = ['yes','y']
OptionB = ['no','n']
while True:        
        print(f"You have specified the protein family, {QueryB} in {QueryA}. Is this correct?")
        QueryX = input("Answer yes or no: ")

        if QueryX.lower() in OptionA:
            print(f"We will begin searching NCBI database for {QueryB} in {QueryA}.")
            # Search against the NCBI database for all Trehalose-6-Phosphate Synthase protein sequences across all mosquitoes species under the taxonomic group Culicidae (TaxonID 7157)
            print(f"esearch -db protein -query \"{QueryA}[Organism:exp] AND {QueryB}.")
            subprocess.call(f"esearch -db protein -query \"{QueryA}[Organism:exp] AND {QueryB}\"",shell=True)
            break # Must have a way of parsing an error message if no sequence count is found
        if QueryX.lower() in OptionB:
            print("Too bad...")
            break
        else:
            print("Answer yes or no only.")
            continue

            # Search against the NCBI database for all Trehalose-6-Phosphate Synthase protein sequences across all mosquitoes species under the taxonomic group Culicidae (TaxonID 7157) and fetch their UID#
            # subprocess.call(f"esearch -db protein -query \"{QueryA}[Organism:exp] AND {QueryB}\" | efetch -format uid > ICA2/{QueryB}.gis",shell=True)
            # TempFileA = open(f"ICA2/{QueryB}.gis").read().upper()
            # len(TempFileA)
            # print(TempFileA)

            # Search against the NCBI database for all Trehalose-6-Phosphate Synthase protein sequences across all mosquito species under the taxonomic group Culicidae (TaxonID 7157) and fetch their Accession numbers
            subprocess.call(f"esearch -db protein -query \"{QueryA}[Organism:exp] AND {QueryB}\" | efetch -format acc > ICA2/{QueryB}{QueryA}.acc",shell=True)
            TempFileB = open("ICA2/{QueryB}{QueryA}.acc").read().upper()
            len(TempFileB)
            print(TempFileB)

            # Search against the NCBI database for all Trehalose-6-Phosphate Synthase protein sequences across all mosquito species under the taxonomic group Culicidae (TaxonID 7157) and fetch their Fasta sequences
            subprocess.call(f"esearch -db protein -query \"{QueryA}[Organism:exp] AND {QueryB}\" | efetch -format fasta > ICA2/{QueryB}{QueryA}.fasta",shell=True)
            TempFileC = open("ICA2/{QueryB}{QueryA}.fasta").read().upper()
            len(TempFileC)
            print(TempFileC)

            # Fetch the protein sequence in fasta based on the user-input Accession numbers
            QueryS = input('Select the accession number for the protein sequence to align against the rest of the sequences: ')
            print(f"You selected the protein sequence {QueryS}")
            subprocess.call(f"esearch -db protein -query \"{QueryS}[Accession]\" | efetch -format fasta > ICA2/{QueryS}.fasta",shell=True)
            TempFileD = open(f"ICA2/{QueryS}.fasta").read().upper()
            len(TempFileD)
            print(f"Showing you the protein sequence for the protein sequence {QueryS}...")
            print(TempFileD)

            # Create a Blast database by indexing all the protein sequences
            # Indexing our database will make our alignment work so much faster
            subprocess.call("makeblastdb -in \"ICA2/{QueryS}.fasta\" -title \"{QueryS}\" -dbtype prot -out ICA2/{QueryS}",shell=True)

            # Align the user-input protein sequence with all the other protein sequences using Blastp
            # Print the blast alignment output in text
            subprocess.call(f"blastp -query \"ICA2/{QueryS}.fasta\" -db \"ICA2/{QueryS}\" -out \"ICA2/{QueryS}Blasted.txt\" -outfmt 7",shell=True)
            TempFileE = open(f"ICA2/{QueryS}Blasted.txt").read().upper()
            print(TempFileE)

            # Align the user-input protein sequence with all the other protein sequences using Clustalo
            # Print the Clustalo alignment output in fasta
            subprocess.call(f"clustalo -i \"ICA2/{QueryS}.fasta\" -o \"ICA2/{QueryS}Alignment.fasta\" --auto -v",shell=True)
            TempFileF = open(f"ICA2/{QueryS}Alignment.fasta").read().upper()
            print(TempFileF)

