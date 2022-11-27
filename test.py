            # If the input accession number is not available in the acc file, the script will loop back to the top
            # And ask the user to re-input the correect accession number
            time.sleep(1)
            print(f"You selected the protein sequence {QueryD}.")

            def CheckAcc(f"ICA2/QueryBr}/QueryBr}.fasta",seachStr)
            isFound = False # Initialise
            with open(f"ICA2/{QueryBr}/{queryBr}.fasta") as AccFile:
                AccVal = AccFile.readlines()
                for StrLine in AccVal:
                    if {QueryD} in StrLine:
                        time.sleep(1)
                        print(f"Fetching the protein sequence {QueryD}...")
                        subprocess.call(f"esearch -db protein -query \"{QueryD}[Accession]\" | efetch -format fasta > \"ICA2/{QueryBr}/{QueryD}.fasta\" ",shell=True) # Fetch protein sequence in fasta format to be shown on screen
                        TempFileD = open(f"ICA2/{QueryBr}/{QueryD}.fasta").read().upper() # Fetch protein sequence in gb format as this will be use dfor further analysis
                        print(f"Would you like me to show you the protein sequence {QueryD} you selected in Fasta?")
                        time.sleep(1)
                        isFound = True

                    else:
                        print("Accession number not found. Input a valid accession number.")
                        return isFound

        # ASK USER TO VIEW CHOSEN PROTEIN SEQUENCE
        while RunLoopXe: 
            QueryE = input("Answer yes or no: ") 
            if QueryE.lower() in OptionA:
                time.sleep(1)
                print(f"This is the protein sequence for {QueryD} in Fasta: ")
                time.sleep(1)
                print(TempFileD)
                time.sleep(1)
                print(f"I have saved this protein sequence in Fasta format as <ICA2/{QueryBr}/{QueryD}.fasta>")
                print(f"And also a copy in Genbak format as <ICA2/{QueryBr}/{QueryD}.gb>")
                time.sleep(3)
                break 
            elif QueryE.lower() in OptionB:
                print("We can view the protein sequence later.")
                break
            else:
                print("Answer yes or no only.")
                continue
                                                                                                                                219,65        41%
