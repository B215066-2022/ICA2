#!/usr/bin/python3

import os, sys, subprocess, shutil

# export NCBI_API_KEY=54b38a82050346352dbaf27e8e02a65f0408" >> ~/.bashrc
# export NCBI_API_KEY=54b38a82050346352dbaf27e8e02a65f0408" >> ~/.bash_profile
# echo ${NCBI_API_KEY}

# Make a working directory first, and let's call it ICA2
os.mkdir("ICA2")

# Search against the NCBI database for all Trehalose-6-Phosphate Synthase protein sequences across all mosquitoes species under the taxonomic group Culicidae (TaxonID 7157)
subprocess.call("esearch -db protein -query \"txid7157[Organism:exp] AND Trehalose-6-Phosphate Synthase\"",shell=True)

# Search against the NCBI database for all Trehalose-6-Phosphate Synthase protein sequences across all mosquitoes species under the taxonomic group Culicidae (TaxonID 7157) and fetch their UID#
subprocess.call("esearch -db protein -query \"txid7157[Organism:exp] AND Trehalose-6-Phosphate Synthase\" | efetch -format uid > ICA2/T6PSynthase.gis",shell=True)
TempFileA = open("ICA2/T6PSynthase.gis").read().upper()
len(TempFileA)
print(TempFileA)

# Search against the NCBI database for all Trehalose-6-Phosphate Synthase protein sequences across all mosquito species under the taxonomic group Culicidae (TaxonID 7157) and fetch their Accession numbers
subprocess.call("esearch -db protein -query \"txid7157[Organism:exp] AND Trehalose-6-Phosphate Synthase\" | efetch -format acc > ICA2/T6PSynthase.acc",shell=True)
TempFileB = open("ICA2/T6PSynthase.acc").read().upper()
len(TempFileB)
print(TempFileB)

# Search against the NCBI database for all Trehalose-6-Phosphate Synthase protein sequences across all mosquito species under the taxonomic group Culicidae (TaxonID 7157) and fetch their Fasta sequences
subprocess.call("esearch -db protein -query \"txid7157[Organism:exp] AND Trehalose-6-Phosphate Synthase\" | efetch -format fasta > ICA2/T6PSynthase.fasta",shell=True)
TempFileC = open("ICA2/T6PSynthase.fasta").read().upper()
len(TempFileC)
print(TempFileC)

# Fetch the protein sequence in fasta based on the user-input Accession numbers
Seq = input('Select the accession number for the protein sequence to align against the rest of the sequences: ')
print(f"You selected the protein sequence {Seq}")
subprocess.call(f"esearch -db protein -query \"{Seq}[Accession]\" | efetch -format fasta > ICA2/{Seq}.fasta",shell=True)
TempFileD = open(f"ICA2/{Seq}.fasta").read().upper()
len(TempFileD)
print(TempFileD)

# Create a Blast database for the protein sequences we will align later
# Indexing our database will make our alignment work so much faster
subprocess.call("makeblastdb -in \"ICA2/T6PSynthase.fasta\" -title \"T6PSynthase\" -dbtype prot -out ICA2/T6PSynthase",shell=True)

# Align the user-input protein sequence with all the other protein sequences using Blastp
# Print the blast alignment output in text
subprocess.call(f"blastp -query \"ICA2/{Seq}.fasta\" -db \"ICA2/T6PSynthase\" -out \"ICA2/T6PSynthaseBAlignment.txt\" -outfmt 7",shell=True)
TempFileE = open("ICA2/T6PSynthaseBAlignment.txt").read().upper()
print(TempFileE)

# Align the user-input protein sequence with all the other protein sequences using Clustalo
# Print the Clustalo alignment output in fasta
subprocess.call("clustalo -i \"ICA2/T6PSynthase.fasta\" -o \"ICA2/T6PSynthaseCAlignment.fasta\" --auto -v",shell=True)
TempFileF = open("ICA2/T6PSynthaseCAlignment.fasta").read().upper()
print(TempFileF)



# OptionA = ['yes', 'y']
# options = ['no', 'n']
#
# while True:
#    user_input = input('Would you like to select a sequence to run an alignment for? Answer yes or no: ')
#    if user_input.lower() in OptionA:
#        print('user typed yes')
#        break
#    elif user_input.lower() in OptionB:
#        print('user typed no')
#        break
#    else:
#        print('Type yes or no')
#        continue
