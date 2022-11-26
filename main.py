#!/usr/bin/python3

import os, sys, shutil
import re # Useful to help us locate a specific motif in our output files
import subprocess # We will be running a lot of background functions that are non python-based
import time # We will be putting some time gaps between each response/print statement

RunLoopA = True
RunLoopB = True
RunLoopC = True
RunLoopD = True
RunLoopE = True
RunLoopF = True
RunLoopG = True

# Check if the directory ICA2 does not exist, we create one
if not os.path.exists("ICA2"):
	os.mkdir("ICA2")

print("-----------------------------")
print("This programme will help you analyse protein sequences in the organism of your choice.")
time.sleep(1)
print("If you know the NCBI Taxonomy ID for the organism you are interested in, you can specify it or search through the NCBI database for the ID.")
time.sleep(1)

OptionA = ['yes','y']
OptionB = ['no','n']

while RunLoopA:
    print("-----------------------------")
    print("I need you to specify the TaxonID for the organism where your protein may be found.")
    time.sleep(1)
    print("For example: For the mosquito family Cullinicea, input the TaxonID 7157.")
    time.sleep(1)    
    QueryA = input('Please input the taxononmy ID of your organism now: ')
    print("------------------------------")
    print("Now, I need you to specify the family of protein you are interested in.")
    time.sleep(1)
    print("For example: Hydrolase or trehalose-6-phosphate synthase.")
    time.sleep(1)
    QueryB = input('Please input the name of your protein family now: ')
    # Python has a weird way of parsing user-input responses containing spaces into filenames in which the output names  will either be truncated or miss the intended file extensions
    # Replacing these spaces into underscores should fix the problem
    QueryBr = QueryB.replace(' ','_')
    
    print("------------------------------")
    print(f"You have specified the protein family {QueryB} across all organism species under the TaxonID {QueryA}.")
    time.sleep(1)
    print("Is this correct?")
    time.sleep(1)
    QueryX = input("Answer yes or no: ")
    
    if QueryX.lower() in OptionA:
        # Make a directory for the given protein family name so that we can store all analysis output files into the same directory
        if not os.path.exists(f"ICA2/{QueryBr}"):
            os.mkdir(f"ICA2/{QueryBr}")
        print(f"------------------------------")
        print(f"We will now begin searching the NCBI database for protein entries of {QueryB} across all organism species under the TaxonID {QueryA}.")
        time.sleep(1)
        print(f"Fetching protein entries...")
        time.sleep(3)
        # Search against the NCBI database for the user-input protein family across all species under the specified taxon group
        subprocess.call(f"esearch -db protein -query \"txid{QueryA}[Organism:exp] and {QueryB}\" > \"ICA2/{QueryBr}/CountsCheck.txt\" ",shell=True)
        TempFile0 = open(f"ICA2/{QueryBr}/CountsCheck.txt").read().upper()
        print(TempFile0)
        print("------------------------------")
        print(f"We will now begin fetching all the entries in the form of UID for {QueryB} across all organism species under the TaxonID {QueryA}.")
        time.sleep(1)
        print("Fetching UID entries (this may take a moment)...")
        time.sleep(3)
        # Search against the NCBI database for the user-input protein family across all species under the specified taxon group and fetch their UID#
        subprocess.call(f"esearch -db protein -query \"txid{QueryA}[Organism:exp] AND {QueryB}\" | efetch -format uid > \"ICA2/{QueryBr}/{QueryBr}.gis\" ",shell=True)
        TempFileA = open(f"ICA2/{QueryBr}/{QueryBr}.gis").read().upper()
        print(TempFileA)
        print("------------------------------")
        print(f"We will now begin fetching all the protein entries in the form of accession numbers for {QueryB} across all orgaism species under the TaxonID {QueryA}.")
        time.sleep(1)
        print("Fetching accession numbers (this may take a moment)...")
        time.sleep(3)   
        # Search against the NCBI database for the user-input protein family across all species under the specified taxon group and fetch their Accession numbers
        subprocess.call(f"esearch -db protein -query \"txid{QueryA}[Organism:exp] AND {QueryB}\" | efetch -format acc > \"ICA2/{QueryBr}/{QueryBr}.acc\" ",shell=True)
        TempFileB = open(f"ICA2/{QueryBr}/{QueryBr}.acc").read().upper()
        print(f"I will now list down all the accession numbers for {QueryB} across all organism species under the TaxonID {QueryA}.")
        time.sleep(1)
        print("Listing down the accession numbers...")
        time.sleep(1)
        print(TempFileB)
        time.sleep(1)
        print(f"I have saved this list as <ICA2/{QueryBr}/{QueryBr}.acc>")
        time.sleep(1)

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

        while RunLoopB:
                if QueryC.lower() in OptionA:
                    time.sleep(1)
                    print(f"This is the list of all protein sequences in FASTA for {QueryB} across all organim species under the TaxonID {QueryA}.")
                    time.sleep(1)
                    print(TempFileC)
                    time.sleep(1)
                    print(f"I have saved this protein sequence in Fasta as <<ICA2/{QueryBr}/{QueryBr}.fasta>")
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

        # Fetch the protein sequence in Fasta based on the user-input Accession number
        print("------------------------------")
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

        while RunLoopC:
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

        while RunLoopD:
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

        while RunLoopE:
            if QueryG.lower() in OptionA:
                time.sleep(1)
                # Align the user-input protein sequence with all the other protein sequences using Clustalo
                # Print the Clustalo alignment output in fasta
                print(f"We will begin running MSA with ClustalO.")
                time.sleep(1)
                print(f"Running MSA with ClustalO on all {QueryB} protein sequences...")
                time.sleep(3)
                subprocess.call(f"clustalo -i \"ICA2/{QueryBr}/{QueryBr}.fasta\" -o \"ICA2/{QueryBr}/{QueryBr}_Aligned.msf\" --auto -v",shell=True)
                print("ClustalO alignment complete!")
                TempFileF = open(f"ICA2/{QueryBr}/{QueryBr}_Aligned.msf").read().upper()
                time.sleep(1)
                print(f"Would you like to view the result of our ClustalO alignment?")
                time.sleep(1)
                QueryH = input("Answer yes or no: ")

                while RunLoopF:
                    if QueryH.lower() in OptionA:
                        time.sleep(1)
                        print(f"This is the ClustalO alignment output for {QueryBr}: ")
                        time.sleep(1)
                        print(TempFileF)
                        time.sleep(3)
                        print(f"I have saved this ClustalO alignment output as <ICA2/{QueryBr}/{QueryBr}_Aligned.msf>")
                        time.sleep(1)
                        break
                    elif QueryH.lower() in OptionB:
                        print("Okay, but you can still view this ClustalO alignment manually in the directory I showed you, next time if you want!")
                        break
                    else:
                        print("Answer yes or no only.")
                        continue

            elif QueryG.lower() in OptionB:
                print("We can run MSA with ClustalO some other time.")
                break
            else:
                print("Answer yes or no only.")
                continue

            # We should make this run in a separate python script
            print("------------------------------")
            print(f"If you would like to take a step further, I can help you to perform further analysis using your selected protein sequence {QueryD}.")
            time.sleep(1)

            while RunLoopG:
                QueryI = input("Select one protein analysis from the following:\n1) Plot interactive multiple sequence alignment (for publication purposes)\n2) Analyse consensus in multiple sequence alignment\n3) Plot conserved regions in amino acid residues\n4) Search motifs in amino acid residues\n5) Determine charges in amino residues\n6) Predict protein secondary structure\n7) Predict coiled-coil likelihood\n8) Predict helical wheel\n0) Start over\n99) Do nothing\nInput your selection: ")
                if QueryI == '1': # For "Plot interactive multiple sequence alignment” with EMBOSS showalign
                    time.sleep(1)
                    print("You have selected to plot interactive multiple sequence alignment.")
                    time.sleep(1)
                    print(f"I will help you to plot MSA for the protein sequence {QueryBr} across all organism species under the TaxonID {QueryA} and report the output in html.")
                    time.sleep(1)
                    print("Plotting the alignment...")
                    time.sleep(3)
                    # EMBOSS showalign
                    subprocess.call(f"showalign -html -high \"4-13 green 43-43 red 51-56 blue\" \"ICA2/{QueryBr}/{QueryBr}_Aligned.msf\" -outfile \"ICA2/{QueryBr}/{QueryBr}_ShowAlign.html\" ", shell=True) # Take msf alignment file for showalign and generate output in html
                    print("Plotting complete!")
                    time.sleep(1)
                    TempFileG = open(f"ICA2/{QueryBr}/{QueryBr}_ShowAlign.html").read()
                    time.sleep(1)
                    print("Showing you the alignment output...")
                    time.sleep(1)
                    print(TempFileG)
                    time.sleep(3)
                    print(f"I have saved this output as <ICA2/{QueryBr}/{QueryD}_ShowAlign.html>")
                    time.sleep(1)
                    print("------------------------------")
                    continue

                elif QueryI == '2': # For “Analyse consensus in multiple sequence alignment” with EMBOSS showalign
                    time.sleep(1)
                    print("You have selected to analyse consensus in multiple sequence alignment.")
                    time.sleep(1)
                    print(f"I will help you to analyse the consensus sequences in the MSA out for {QueryB} across all organism species under the TaxonID {QueryA} and report the output in html.")
                    time.sleep(1)
                    print("Analysing the consensus sequences...")
                    time.sleep(3)
                    # EMBOSS infoalign
                    subprocess.call(f"infoalign \"ICA2/{QueryBr}/{QueryBr}_Aligned.msf\" -html -outfile \"ICA2/{QueryBr}/{QueryBr}_InfoAlign.html\" ",shell=True) 
# Take msf alignment file for showalign and generate output in html
                    print("Analysis complete!")
                    time.sleep(1)
                    TempFileH = open(f"ICA2/{QueryBr}/{QueryBr}_InfoAlign.html").read()
                    time.sleep(1)
                    print("Showing you the analysis output...")
                    time.sleep(1)
                    print(TempFileH)
                    time.sleep(3)
                    print(f"I have saved this output as <ICA2/{QueryBr}/{QueryBr}_InfoAlign.html>")
                    time.sleep(1)
                    print("------------------------------")
                    continue

                elif QueryI == '3': # For “Plot conserved regions in amino acid residues” with EMBOSS plotcon
                    time.sleep(1)
                    print("You have selected to plot conserved regions in amino acid residues.")
                    time.sleep(1)
                    print(f"I will help you to identify the conserved amino acid residues between all the protein sequences for {QueryB} across all organism species under the TaxonID {QueryA} report the output in graph.")
                    time.sleep(1)
                    print("Plotting the graph...")
                    time.sleep(3)
                    # EMBOSS plotcon
                    subprocess.call(f"plotcon -sequences \"ICA2/{QueryBr}/{QueryBr}_Aligned.msf\" -winsize 40 -graph png -goutfile \"{QueryBr}_Plotcon\" -gdirectory \"ICA2/{QueryBr}\" -scorefile EBLOSUM62 ",shell=True) # Take msf alignment file for plotcon and generate graph output in png format (default)
                    print("Plotting complete!")
                    time.sleep(1)
                    print("Showing you the graph output...")
                    time.sleep(1)
                    # Show the output graph in a pop-up window
                    import matplotlib.pyplot as plt
                    import matplotlib.image as mpimg
                    TempGraph1 = mpimg.imread(f'ICA2/{QueryBr}/{QueryBr}_Plotcon.1.png') # By default, plotcon will add number 1 in the output file name so we have to specify it for path finding
                    plt.imshow(TempGraph1)
                    plt.show()
                    print(f"I have saved this graph output as <ICA2/{QueryBr}/{QueryBr}_Plotcon.1.png>")
                    time.sleep(3)
                    print("------------------------------")
                    continue

                elif QueryI == '4': # For "Search motifs in amino acid residues" with EMBOSS patmatmotifs
                    time.sleep(1)
                    print("You have selected to search motifs in amino acid residues.")
                    time.sleep(1)
                    print(f"I will help you to search for interesting motifs in the protein sequence {QueryD} against the motifs available in PROSITE database.")
                    time.sleep(1)
                    print(f"Searching for motifs in the protein sequence {QueryD}...")
                    time.sleep(3)
                    # EMBOSS patmatmotifs
                    subprocess.call(f"patmatmotifs -full -sprotein1 -sformat gb \"ICA2/{QueryBr}/{QueryD}.gb\" -outfile \"ICA2/{QueryBr}/{QueryD}_Motif.dbmotif\" ",shell=True) # Use input {QueryD} protein in gb format and generate output file in default dbmotif format
                    print("Motifs search complete!")
                    TempFileI = open(f"ICA2/{QueryBr}/{QueryD}_Motif.dbmotif").read()
                    time.sleep(1)
                    print("Showing you the motifs output...")
                    time.sleep(1)
                    print(TempFileI)
                    time.sleep(3)
                    print(f"I have saved this output as <ICA2/{QueryBr}/{QueryD}_Motif.dbmotif>")
                    time.sleep(1)
                    print("------------------------------")
                    continue

                elif QueryI == '5': # For "Determine charges in amino acid residues" with EMBOSS charge
                    time.sleep(1)
                    print("You have selected to determine the charges in amino acid residues.")
                    time.sleep(1)
                    print(f"I will help you to determine the charge for every 15 amino acid residues in the protein sequence {QueryD} and report the output in graph.")
                    time.sleep(1)
                    print("Determining amino acid charges...")
                    time.sleep(1)
                    print("Plotting the graph...")
                    time.sleep(3)
                    # EMBOSS charge
                    subprocess.call(f"charge -sprotein -sformat gb \"ICA2/{QueryBr}/{QueryD}.gb\" -outfile \"ICA2/{QueryBr}/{QueryD}_Charge.charge\" -plot -graph png -goutfile \"{QueryD}_Charge\" -gdirectory \"ICA2/{QueryBr}\" ",shell=True) # Take input {QueryD} protein in gb format and generate output file in default charge format
                    print("Plotting complete!")
                    time.sleep(1)
                    print("Showing you the graph output...")
                    time.sleep(1)
                    # Show the output graph in a pop-up window
                    import matplotlib.pyplot as plt
                    import matplotlib.image as mpimg
                    TempGraph2 = mpimg.imread(f'ICA2/{QueryBr}/{QueryD}_Charge.1.png') # By default, charge will add number 1 in the output file name so we have to specify it for path finding
                    plt.imshow(TempGraph2)
                    plt.show()
                    time.sleep(3)
                    print(f"I have saved this graph output as <ICA2/{QueryBr}/{QueryD}_Charge.1.png>")
                    print(f"And a copy of the raw data as <ICA2/{QueryBr}/{QueryD}_Charge.charge>")
                    time.sleep(1)
                    print("------------------------------")
                    continue

                elif QueryI == '6': # For "Predict secondary structure" with EMBOSS garnier
                    time.sleep(1)
                    print("You have selected to predict secondary structure.")
                    time.sleep(1)
                    print(f"I will help you to predict the secondary structure for the protein sequence {QueryD} and write the report in Genbank format.")
                    time.sleep(1)
                    print("Predicting the secondary structure...")
                    time.sleep(3)
                    # EMBOSS garnier
                    subprocess.call(f"garnier -sprotein1 -sformat gb \"ICA2/{QueryBr}/{QueryD}.gb\" -outfile \"{QueryD}_PredictedStructure.tagseq\" -rdirectory2 \"ICA2/{QueryBr}\" ", shell=True) # Take input {QueryD} protein in gb format and generate output in default tagseq format
                    print
                    print("Prediction complete!")
                    TempFileJ = open(f"ICA2/{QueryBr}/{QueryD}_PredictedStructure.tagseq").read()
                    time.sleep(1)
                    print("Showing you the secondary structure prediction output:")
                    time.sleep(1)
                    print(TempFileJ)
                    time.sleep(3)
                    print(f"I have saved a copy of this report as <ICA2/{QueryBr}/{QueryD}_PredictedStructure.tagseq>")
                    time.sleep(1)
                    print("------------------------------")
                    continue

                elif QueryI == '7': # For "Predict coiled-coil likelihood" with EMBOSS pepcoil
                    time.sleep(1)
                    print("You have selected to predict coiled-coil likelihood.")
                    time.sleep(1)
                    print(f"I will help you to predict the coiled-coil likelihood by calculating the probability of a coiled-coil structure for the amino acid residues in the protein sequence {QueryD}.")
                    time.sleep(1)
                    print(f"Calculating the coiled-coil probability for every 15 residues...")
                    time.sleep(3)
                    # EMBOSS pepcoil
                    subprocess.call(f"pepcoil -sprotein1 -sformat gb \"ICA2/{QueryBr}/{QueryD}.gb\" -window 15 -outfile \"ICA2/{QueryBr}/{QueryD}_Coil.pepcoil\" ",shell=True) # Take input {QueryD} protein in gb format and generate output in default pepcoil format
                    print
                    print("Calculation complete!")
                    TempFileK = open(f"ICA2/{QueryBr}/{QueryD}_Coil.pepcoil").read()
                    time.sleep(1)
                    print(f"Showing you the coiled-coil prediction output:")
                    time.sleep(1)
                    print(TempFileK)
                    time.sleep(3)
                    print(f"I have saved a copy of this output as <ICA2/{QueryBr}/{QueryD}_Coil.pepcoil>")
                    time.sleep(1)
                    print("------------------------------")
                    continue

                elif QueryI == '8': # For "Predict helical wheel” with EMBOSS pepwheel
                    time.sleep(1)
                    print("You have selected to predict helical wheel.")
                    time.sleep(1)
                    print(f"I will help you to predict the helical wheel representation of the protein sequence {QueryD}.")
                    time.sleep(1)
                    print(f"Predicting helical wheel...")
                    time.sleep(3)
                    # EMBOSS pepwheel
                    subprocess.call(f"pepwheel -sprotein1 -sformat gb \"ICA2/{QueryBr}/{QueryD}.gb\" -graph png -goutfile \"{QueryD}_Wheel\" -gdirectory \"ICA2/{QueryBr}\" ", shell=True) # Take input {QueryD} protein in gb format and generate output in default png format
                    print
                    print("Prediction complete!")
                    TempFileK = open(f"ICA2/{QueryBr}/{QueryD}_Coil.pepcoil").read()
                    time.sleep(1)
                    print(f"Showing you the coiled-coil prediction output:")
                    time.sleep(1)
                    print(TempFileK)
                    time.sleep(3)
                    print(f"I have saved a copy of this output as <ICA2/{QueryBr}/{QueryD}_Coil.pepcoil>")
                    time.sleep(1)
                    print("------------------------------")
                    continue

                elif QueryI == '0': # For "Start all over again"
                    break

                elif QueryI == '99': # For "Do nothing"
                    print("Okay, we will not do any further analysis for now...")
                    break
                 
                else:
                    print("Input your selection number from 0 to 8 or 99!")
                    continue

                break

            break

    elif QueryX.lower() in OptionB:
        print("You can run this programme again when you are ready.")
        break

    else:
        print("Answer yes or no only.")
        continue
