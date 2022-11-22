#!/usr/bin/python3

import os
import sys
import shutil
import re # Useful to help us locate a specific motif in our output files
import subprocess # We will be running a lot of background functions that are non python-based
import time # We will be putting some time gaps between each response or print statement

# export "NCBI_API_KEY=54b38a82050346352dbaf27e8e02a65f0408" >> ~/.bashrc
# export "NCBI_API_KEY=54b38a82050346352dbaf27e8e02a65f0408" >> ~/.bash_profile

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
        print("Is this correct?")
        time.sleep(1)
        QueryX = input("Answer yes or no: ")
        if QueryX.lower() in OptionA:
            
            # Make a directory for the given protein family name so that we can store all analysis output files into the same directory
            if not os.path.exists(f"ICA2/{QueryBr}"):
                os.mkdir(f"ICA2/{QueryBr}")
            print(f"------------------------------")
            print(f"We will now begin searching the NCBI database for protein entries of {QueryB} in the organism under the TaxonID {QueryA}.")
            time.sleep(1)
            print(f"Fetching protein entries...")
            time.sleep(3)
            # Search against the NCBI database for the user-input protein family across all species under the specified taxon group
            subprocess.call(f"esearch -db protein -query \"txid{QueryA}[Organism:exp] and {QueryB}\" > \"ICA2/{QueryBr}/CountsCheck.txt\" ",shell=True)
            TempFile0 = open(f"ICA2/{QueryBr}/CountsCheck.txt").read().upper()
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

            print("------------------------------")
            print(f"We will now begin fetching all the entries in the form of UID for {QueryB} in the organism under the TaxonID {QueryA}.")
            time.sleep(1)
            print("Fetching UID entries (this may take a moment)...")
            time.sleep(3)
            # Search against the NCBI database for the user-input protein family across all species under the specified taxon group and fetch their UID#
            subprocess.call(f"esearch -db protein -query \"txid{QueryA}[Organism:exp] AND {QueryB}\" | efetch -format uid > \"ICA2/{QueryBr}/{QueryBr}.gis\" ",shell=True)
            TempFileA = open(f"ICA2/{QueryBr}/{QueryBr}.gis").read().upper()
            print(TempFileA)
            # with TempFileA,'r' as TempFileAB:
            #    TempFileAC = len(TempFileAB.readiness()) # Count how many UID lines are  in the .gis file, meaning how many protein entries we have
            #    print('Total number of protein counts found: ', TempFileAC)
            # Must have a way of parsing an error message if no sequence count is found using re.findall

            print("------------------------------")
            print(f"We will now begin fetching all the protein entries in the form of accession numbers for {QueryB} in the organism under the TaxonID {QueryA}.")
            time.sleep(1)
            print("Fetching accession numbers (this may take a moment)...")
            time.sleep(3)   
            # Search against the NCBI database for the user-input protein family across all species under the specified taxon group and fetch their Accession numbers
            subprocess.call(f"esearch -db protein -query \"txid{QueryA}[Organism:exp] AND {QueryB}\" | efetch -format acc > \"ICA2/{QueryBr}/{QueryBr}.acc\" ",shell=True)
            TempFileB = open(f"ICA2/{QueryBr}/{QueryBr}.acc").read().upper()
            print(f"I will now list down all the accession numbers for {QueryB} in the organism under the TaxonID {QueryA}.")
            time.sleep(1)
            print("Listing down the accession numbers...")
            time.sleep(1)
            print(TempFileB)
            time.sleep(1)
            print(f"I have saved this list as <ICA2/{QueryBr}/{QueryBr}.acc>")
            time.sleep(3)

            print("------------------------------")
            print(f"We will now begin fetching all the protein sequences in for {QueryB} across all organism species under the TaxonID {QueryA}.")
            time.sleep(1)
            print("Fetching the protein sequences (this may take a moment)... ")
            time.sleep(3)
            # Fetch all protein sequences in FASTA for all available accession entries
            subprocess.call(f"esearch -db protein -query \"txid{QueryA}[Organism:exp] AND {QueryB}\" | efetch -format fasta > \"ICA2/{QueryBr}/{QueryBr}.fasta\" ",shell=True)
            subprocess.call(f"esearch -db protein -query \"txid{QueryA}[Organism:exp] AND {QueryB}\" | efetch -format gb > \"ICA2/{QueryBr}/{QueryBr}.gb\" ",shell=True)
            TempFileC = open(f"ICA2/{QueryBr}/{QueryBr}.fasta").read().upper()
            print(f"Would you like me to list down all the protein sequences in available in Fasta? Just a warning, it can be overwhelming if you choose to: ")
            time.sleep(1)
            QueryC = input("Answer yes or no: ")
            
            while True:
                if QueryC.lower() in OptionA:
                    time.sleep(1)
                    print(f"This is the list of all protein sequences in FASTA for {QueryB} across all organim species under the TaxonID {QueryA}.")
                    time.sleep(1)
                    print(TempFileC)
                    time.sleep(1)
                    print("I have saved this list of protein sequences in Fasta format as <ICA2/{QueryBr}/{QueryBr}.fasta>")
                    time.sleep(1)
                    print(f"And also a copy in Genbank format as <ICA2/{QueryBr}/{QueryBr}.gb")
                    time.sleep(3)
                    break
                elif QueryC.lower() in OptionB:
                    print("We can continue next time.")
                    break
                else:
                    print("Answer yes or no only.")
                    continue

            # Fetch the protein sequence in fasta based on the user-input Accession number
            print(f"------------------------------")
            print("Now, I need to know which protein sequence you are interested in analysing.") 
            time.sleep(1)
            QueryD = input('Input the accession number for the protein sequence of your choice: ')
            time.sleep(1)
            print(f"You selected the protein sequence {QueryD}.")
            time.sleep(3)
            print(f"Fetching the protein sequence {QueryD}...")
            subprocess.call(f"esearch -db protein -query \"{QueryD}[Accession]\" | efetch -format fasta > \"ICA2/{QueryBr}/{QueryD}.fasta\" ",shell=True)
            subprocess.call(f"esearch -db protein -query \"{QueryD}[Accession]\" | efetch -format gb > \"ICA2/{QueryBr}/{QueryD}.gb\" ",shell=True)
            TempFileD = open(f"ICA2/{QueryBr}/{QueryD}.fasta").read().upper()
            
            print(f"Would you like me to show you the protein sequence {QueryD} you selected in Fasta?")
            time.sleep(1)
            QueryE = input("Answer yes or no: ")

            while True:
                if QueryE.lower() in OptionA:
                    time.sleep(1)
                    print(f"This is the protein sequence for {QueryD} in Fasta: ")
                    time.sleep(1)
                    print(TempFileD)
                    time.sleep(1)
                    print(f"I have saved this protein sequence in Fasta format as <ICA2/{QueryBr}/{QueryD}.fasta>")
                    print(f"And also a copy in Genbak format as <ICA2/{QueryBr}/{QueryD}.gb>")
                    time.sleep(3) # Allow some time for the user to look at and process the input sequence
                    break
                elif QueryE.lower() in OptionB:
                    print("We can continue next time.")
                    break
                else:
                    print("Answer yes or no only.")
                    continue
            
            print(f"------------------------------")
            print("Great, we have all the data we need for our protein analysis!")
            time.sleep(1)
            print("But before that, we will need to create a database and index our protein sequences.")
            time.sleep(1)
            print(f"Indexing our {QueryB} protein sequences in the organisms under the TaxonID {QueryA}...")
            time.sleep(3)
            print("Generating database...")
            time.sleep(3)
            # Create a Blast database by indexing all the protein sequences in the Fasta file we initially fetched
            # Indexing our database will make our alignment work so much faster
            subprocess.call(f"makeblastdb -in \"ICA2/{QueryBr}/{QueryBr}.fasta\" -title \"{QueryBr}\" -dbtype prot -out \"ICA2/{QueryBr}/{QueryBr}\" ", shell=True) # {QueryBr}.fasta because this file contained the fasta sequence for all protein entries based on the user-input {QueryB}
            print("Database indexed and created!")
            time.sleep(3)

            # Align the user-input protein sequence with all the other protein sequences using Blastp
            # Print the blast alignment output in text
            print("------------------------------")
            print(f"Now, we will begin aligning {QueryD} against the other protein sequences with Blast.")
            time.sleep(1)
            print(f"Running Blast alignment between {QueryD} and all the other {QueryB} protein sequences...")
            time.sleep(3)
            subprocess.call(f"blastp -query \"ICA2/{QueryBr}/{QueryD}.fasta\" -db \"ICA2/{QueryBr}/{QueryBr}\" -out \"ICA2/{QueryBr}/{QueryD}Blasted.txt\" -outfmt 7",shell=True)
            print("Blast alignment complete!")
            TempFileE = open(f"ICA2/{QueryBr}/{QueryD}Blasted.txt").read().upper()
            time.sleep(1)
            print("Would you like to view the result of our Blast alignment?")
            time.sleep(1)
            QueryF = input("Answer yes or no: ")

            while True:
                if QueryF.lower() in OptionA:
                    time.sleep(1)
                    print(f"This is the Blast alignment output for {QueryD}: ")
                    time.sleep(1)
                    print(TempFileE)
                    time.sleep(3)
                    print(f"I have saved this Blast alignment output as <ICA2/{QueryBr}/{QueryD}Blasted.txt>")
                    time.sleep(1)
                    break
                elif QueryF.lower() in OptionB:
                    print("Okay, but you can still view this Blast alignment manually in the directory I showed you, next time if you want!")
                    break
                else:
                    print("Answer yes or no only.")
                    continue

            print("------------------------------")
            print(f"Would you like to run mutiple sequence alignment (MSA) for all the protein sequences for {QueryD} in the organism  under the TaxonID {QueryA} with ClustalO now?")
            time.sleep(1)
            QueryG = input("Answer yes or no: ")

            while True:
                if QueryG.lower() in OptionA:
                    time.sleep(1)
                    # Align the user-input protein sequence with all the other protein sequences using Clustalo
                    # Print the Clustalo alignment output in fasta
                    print(f"We will begin running MSA with ClustalO.")
                    time.sleep(1)
                    print(f"Running MSA with ClustalO on all {QueryB} protein sequences...")
                    time.sleep(3)
                    subprocess.call(f"clustalo -i \"ICA2/{QueryBr}/{QueryBr}.fasta\" -o \"ICA2/{QueryBr}/{QueryBr}Aligned.msf\" --auto -v",shell=True)
                    print("ClustalO alignment complete!")
                    TempFileF = open(f"ICA2/{QueryBr}/{QueryBr}Aligned.msf").read().upper()
                    time.sleep(1)
                    print(f"Would you like to view the result of our ClustalO alignment?")
                    time.sleep(1)
                    QueryH = input("Answer yes or no: ")

                    while True:
                        if QueryH.lower() in OptionA:
                            time.sleep(1)
                            print(f"This is the ClustalO alignment output for {QueryBr}: ")
                            time.sleep(1)
                            print(TempFileF)
                            time.sleep(3)
                            print(f"I have saved this ClustalO alignment output as <ICA2/{QueryBr}/{QueryBr}Aligned.msf>")
                            time.sleep(1)
                            break
                        elif QueryH.lower() in OptionB:
                            print("Okay, but you can still view this ClustalO alignment manually in the directory I showed you, next time if you want!")
                            break
                        else:
                            print("Answer yes or no only.")
                            continue    
                    break
                elif QueryG.lower() in OptionB:
                    print("We can run MSA with ClustalO some other time.")
                    break
                else:
                    print("Answer yes or no only.")
                    continue
            
            # We should make this run in a separate python shell later
            print("------------------------------")
            print(f"If you would like to take a step further, I can help you to perform further analysis using your selected protein sequence {Queryd}.")
            time.sleep(1)
            QueryH = input("Select one analysis from the following:\n1) Plot conserved regions\n2) Search motifs\n3) Extract motifs")

            while True:
                if QueryH == '1': # For "Plot conserved regions" with EMBOSS plotcon
                    time.sleep(1)
                    print("You have selected to plot conserved regions.")
                    time.sleep(1)
                    print(f"I will help you to identify all possible conserved regions between all the protein sequences for {QueryB} across all organism species under the TaxonID {QueryA} in the form of a graph.")
                    time.sleep(1)
                    print("Plotting the graph....")
                    time.sleep(3)
                    # EMBOSS plotcon
                    subprocess.call(f"plotcon -sequences \"ICA2/{QueryBr}/{QueryBr}Aligned.msf\" -winsize 100 -graph png ",shell=True)
                    # Something is wrong here:
                    # The graph output plotcon.1.png is being saved in the home directory rather than the ICA2 directory 
                    # Quick fix is to move the plotcon.1.png file to the directory we want
                    os.rename("plotcon.1.png","ICA2/PlotCon1.png") 
                    # But the problem with this is that the second time we run plotcon, it will generate another plotcon.1.png in the home directory and we won't be able to move it the the same directory we want as a separate file because it has the same name as the one we moved earlier!
                    # Unless, if we have a way of making the programme know that we run a second plotcon and make the output file saved as plotcon2.png. Hmmm but how do we do that?
                    print("Plotting complete!")
                    time.sleep(1)
                    print("Showing you the graph output...")
                    time.sleep(1)
                    # Show the output graph in a pop-up window
                    import matplotlib.pyplot as plt
                    import matplotlib.image as mpimg
                    img = mpimg.imread(f'ICA2/PlotCon1.png')
                    plt.imshow(img)
                    plt.show()
                    print("I have saved this graph output as <ICA2/PlotCon1.png>")
                    time.sleep(3)
                    break

                elif QueryH == '2': # For "Search motifs" with EMBOSS patmatmotifs
                    time.sleep(1)
                    print("You have selected to search sequence motifs.")
                    time.sleep(1)
                    print(f"I will help you to search for motifs in the protein sequence {QueryD} against PROSITE database.")
                    time.sleep(1)
                    print("Indexing PROSITE database...") # By indexing, we mean downloading the motifs file prosite.dat
                    time.sleep(3)
                    subprocess.call("wget -P \"ICA2\" \"ftp://ftp.expasy.org/databases/prosite/prosite.dat\" ",shell=True) # Download the file prosite.dat and save it in the ICA2 directory
                    print("PROSITE databased indexed!")
                    time.sleep(1)
                    # EMBOSS patmatmotifs
                    print(f"Searching for motifs in the protein sequence {QueryD}...")
                    time.sleep(3)
                    subprocess.call(f"patmatmotifs -sprotein1 -sformat gb \"ICA2/{QueryBr}/{QueryD}.gb\" -full -outfile \"ICA2/{QueryBr}/{QueryD}.patmatmotifs\" ", shell=True) # Let's standardise the input file in Genbank (.gb) format
                    print
                    print("Motif search complete!")
                    time.sleep(1)
                    break

                elif QueryH == '3': # For "Extract motifs" with EMBOSS extractalign
                    time.sleep(1)
                    print("You have selected to extract motifs")
                    break
                
                else:
                    print("Input your selection number 1-3")
                    continue

            break
        elif QueryX.lower() in OptionB:
            print("You can run this programme again when you are ready.")
            break
        else:
            print("Answer yes or no only.")
            continue
