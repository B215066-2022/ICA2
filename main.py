#!/usr/bin/python3

import os
import sys
import shutil
import re # Useful to help us locate a specific motif in our output files
import subprocess # We will be running a lot of background functions that are non python-based
import time # We will be putting some time gaps between each response or print statement

# export NCBI_API_KEY=54b38a82050346352dbaf27e8e02a65f0408" >> ~/.bashrc
# export NCBI_API_KEY=54b38a82050346352dbaf27e8e02a65f0408" >> ~/.bash_profile
# echo ${NCBI_API_KEY}

# Make a working directory first, and let's call it ICA2
if not os.path.exists("ICA2"):
    os.mkdir("ICA2")

# We need a prompt-dialogue here to ask whether the user wants to input an organism name or TaxonID.

print("This programme will help you analyse protein sequences of your choice.")
time.sleep(1)
print("-----------------------------")
print("First, I need you to specify the NCBI taxonomy ID for the organism where your protein may be found.")
time.sleep(1)
print("For example, 7157 or 4130.")
time.sleep(1)
QueryA = input('Please input the taxononmy ID of your organism now: ')
print("------------------------------")
print("Now, I need you to specify the family of protein you are interested in.")
time.sleep(1)
print("For example, hydrolase or glucose-6-phosphate dehydrogenase.")
time.sleep(1)
QueryB = input('Please input the name of your protein family now: ')
# Programming language has a weird way of parsing user-input responses containing spaces into filenames in which the output names  will either be truncated or miss the intended file extensions
# Replacing these spaces into underscores should fix the problem
QueryBr = QueryB.replace(' ','_')

OptionA = ['yes','y']
OptionB = ['no','n']
while True:
        print(f"------------------------------")
        print(f"You have specified the protein family {QueryB} in the organism under the TaxonID {QueryA}.")
        time.sleep(1)
        print(f"Is this correct?")
        time.sleep(1)
        QueryX = input("Answer yes or no: ")
        if QueryX.lower() in OptionA:
            
            print(f"------------------------------")
            print(f"We will now begin searching the NCBI database for protein entries of {QueryB} in the organism under the TaxonID {QueryA}.")
            time.sleep(3)
            # Search against the NCBI database for the user-input protein family across all species under the specified taxon group
            subprocess.call(f"esearch -db protein -query \"txid{QueryA}[Organism:exp] and {QueryB}\" > \"ICA2/{QueryBr}.txt\" ",shell=True)
            TempFile0 = open(f"ICA2/{QueryBr}.txt").read().upper()
            print(TempFile0)
            # Check if the search returns any counts for the user-input protein family
            # for row in TempFile0:
            #    column = row[0]  # The first column of this row.
            #    value = float(column)  # The text returns strings, so we should turn them into floats for numeric comparison.
            #    if value > 0:
            #        print(f"We found {value} counts for this protein.")
            #    elif value < 1:
            #        print (f"We could not find any protein based on the search results. Try a different protein family.")
            #    else:
            #        print("...")

            print(f"------------------------------")
            print(f"We will now begin fetching all the entries in the form of UID for {QueryB} in the organism under the TaxonID {QueryA}")
            time.sleep(3)
            # Search against the NCBI database for the user-input protein family across all species under the specified taxon group and fetch their UID#
            subprocess.call(f"esearch -db protein -query \"txid{QueryA}[Organism:exp] AND {QueryB}\" | efetch -format uid > \"ICA2/{QueryBr}.gis\" ",shell=True)
            TempFileA = open(f"ICA2/{QueryBr}.gis").read().upper()
            print(TempFileA)
            # with TempFileA,'r' as TempFileAB:
            #    TempFileAC = len(TempFileAB.readiness()) # Count how many UID lines are  in the .gis file, meaning how many protein entries we have
            #    print('Total number of protein counts found: ', TempFileAC)
            # Must have a way of parsing an error message if no sequence count is found using re.findall

            print(f"------------------------------")
            print(f"We will now begin fetching all the protein entries in the form of accession numbers for {QueryB} in the organism under the TaxonID {QueryA}")
            time.sleep(3)   
            # Search against the NCBI database for the user-input protein family across all species under the specified taxon group and fetch their Accession numbers
            subprocess.call(f"esearch -db protein -query \"txid{QueryA}[Organism:exp] AND {QueryB}\" | efetch -format acc > \"ICA2/{QueryBr}.acc\" ",shell=True)
            TempFileB = open(f"ICA2/{QueryBr}.acc").read().upper()
            print(f"I will now list down all the accession numbers for {QueryB} in the organism under the TaxonID {QueryA}")
            print(TempFileB)

            

            print(f"------------------------------")
            print(f"We will now begin fetching all the protein sequences in FASTA.")
            time.sleep(3)
            # Fetch all protein sequences in FASTA for all available accession entries
            subprocess.call(f"esearch -db protein -query \"txid{QueryA}[Organism:exp] AND {QueryB}\" | efetch -format fasta > \"ICA2/{QueryBr}.fasta\" ",shell=True)
            TempFileC = open(f"ICA2/{QueryBr}.fasta").read().upper()
            print(f"Would you like me to list down all the FASTA sequences? Just a warning, it can be overwhelming if you choose to: ")
            time.sleep(1)
            QueryC = input("Answer yes or no:")
            
            while True:
                if QueryC.lower() in OptionA:
                    time.sleep(1)
                    print(f"This is the list of all protein sequences for {QueryB} in the organism under the TaxonID {QueryA}.")
                    time.sleep(1)
                    print(TempFileC)
                    break
                elif QueryC.lower() in OptionB:
                    print("We can continue next time.")
                    break
                else:
                    print("Answer yes or no only.")
                    continue

            # Fetch the protein sequence in fasta based on the user-input Accession number
            print("Now, I need to know which protein sequence are you interested in analysing.") 
            time.sleep(1)
            QueryD = input('Input the accession number for the protein sequence of your choice: ')
            print(f"You selected the protein sequence {QueryD}")
            time.sleep(3)
            print(f"Fetching the protein sequence {QueryD}...")
            subprocess.call(f"esearch -db protein -query \"{QueryD}[Accession]\" | efetch -format fasta > \"ICA2/{QueryD}.fasta\" ",shell=True)
            TempFileD = open(f"ICA2/{QueryD}.fasta").read().upper()
            
            print(f"Would you like me to show you the protein sequence {QueryD} you selected?")
            time.sleep(1)
            QueryE = input("Answer yes or no:")

            while True:
                if QueryE.lower() in OptionA:
                    time.sleep(1)
                    print(f"This is the protein sequence for {QueryD}: ")
                    time.sleep(1)
                    print(TempFileD)
                    break
                elif QueryE.lower() in OptionB:
                    print("We can continue next time.")
                    break
                else:
                    print("Answer yes or no only.")
                    continue

            # Create a Blast database by indexing all the protein sequences
            # Indexing our database will make our alignment work so much faster
            # subprocess.call("makeblastdb -in \"ICA2/{QueryS}.fasta\" -title \"{QueryS}\" -dbtype prot -out \"ICA2/{QueryS}\" ",shell=True)
            # Add a print statement here stating that indexing database is needed
            # Align the user-input protein sequence with all the other protein sequences using Blastp
            # Print the blast alignment output in text
            # subprocess.call(f"blastp -query \"ICA2/{QueryS}.fasta\" -db \"ICA2/{QueryS}\" -out \"ICA2/{QueryS}Blasted.txt\" -outfmt 7",shell=True)
            # TempFileE = open(f"ICA2/{QueryS}Blasted.txt").read().upper()
            # print(TempFileE)

            # Align the user-input protein sequence with all the other protein sequences using Clustalo
            # Print the Clustalo alignment output in fasta
            # subprocess.call(f"clustalo -i \"ICA2/{QueryS}.fasta\" -o \"ICA2/{QueryS}Alignment.fasta\" --auto -v",shell=True)
            break

        elif QueryX.lower() in OptionB:
            print("You can run this programme again when you are ready.")
            break
        else:
            print("Answer yes or no only.")
            continue
