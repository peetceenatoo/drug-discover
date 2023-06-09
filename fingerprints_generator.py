# ------------------------------- Imports ------------------------------- #

from FPSim2.io import create_db_file

# ------------------------------- Code ------------------------------- #

# Specify the paths for the database
dataset_paths = ["\\dataset\\Commercial_MW\\Commercial_MWlower330.csv",
                 "\\dataset\\Commercial_MW\\Commercial_MW330-500-1.csv",
                 "\\dataset\\Commercial_MW\\Commercial_MW330-500-2.csv",
                 "\\dataset\\Commercial_MW\\Commercial_MWhigher500.csv"]
out1 = "fingerprints.h5"

# List of input smiles strings
input_smiles = []

# For each dataset in the list
i=1
print("Start reading all the smiles...")
for path in dataset_paths:

    # Open the files
    f = open(path,"r")

    # Read all the smiles from f
    for x in f:
        # Append current SMILES at index i
        input_smiles.append([x.replace("\n",""), 0+i])
        i +=1
        
    f.close()

# Create the .h5 file
print("Start creating the database... Good luck!")
create_db_file(input_smiles, out1, 'RDKit')
