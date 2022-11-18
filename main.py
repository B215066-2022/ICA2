#!/usr/bin/python3

import os, sys, subprocess, shutil

# export NCBI_API_KEY=54b38a82050346352dbaf27e8e02a65f0408" >> ~/.bashrc
# export NCBI_API_KEY=54b38a82050346352dbaf27e8e02a65f0408" >> ~/.bash_profile
# echo ${NCBI_API_KEY}

os.mkdir("ICA2")

# Search against the NCBI database for all Trehalose-6-Phosphate Synthase protein entries across all mosquitoes species under the taxonomic group Culicidae (TaxonID 7157)
subprocess.call("esearch -db protein -query \"txid7157[Organism:exp] AND Trehalose-6-Phosphate Synthase\" ")


# Search against the NCBI database for all Trehalose-6-Phosphate Synthase UID# across all mosquitoes species under the taxonomic group Culicidae (TaxonID 7157)
subprocess.call("esearch -db protein -query \"txid7157[Organism:exp] AND Trehalose-6-Phosphate Synthase\" | efetch -format uid > ICA2/T6PS.gis",shell=True)
TempFileA = open("ICA2/T6PS.gis").read().upper()
len(TempFileA)
# head -5 ICA2/T6PS.gis
# wc -l ICA2/T6PS.gis


# Search against the NCBI database for all Trehalose-6-Phosphate Synthase Accession# across all mosquitoes species under the taxonomic group Culicidae (TaxonID 7157)
subprocess.call("esearch -db protein -query \"txid7157[Organism:exp] AND Trehalose-6-Phosphate Synthase\" | efetch -format acc > ICA2/T6PS.acc",shell=True)
TempFileB = open("ICA2/T6PS.acc").read().upper()
len(TempFileB)
# head -5 ICA2/T6PS.acc
# wc -l ICA2/T6PS.acc


# Search against the NCBI database for all Trehalose-6-Phosphate Synthase protein entries across all mosquito species under the taxonomic group Culicidae (TaxonID 7157) and fetch their fasta sequences
subprocess.call("esearch -db protein -query \"txid7157[Organism:exp] AND Trehalose-6-Phosphate Synthase\" | efetch -format fasta > ICA2/T6PS.fasta3",shell=True)
TempFileC = open("ICA2/T6PS.fasta3").read().upper()
len(TempFileC)
# head -5 ICA2/T6PS.fasta3
# grep -c ">" ICA2/T6PS.fasta3
