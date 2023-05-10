# --------------------- Solve errors and warnings ---------------------- #

import array
import os
import sys
import logging
from contextlib import redirect_stdout
from tensorflow.python.util import deprecation

# Disable logging output of tensorflow content [May be useless] 
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'

# Do not print deprecation warnings of tensorflow content
deprecation._PRINT_DEPRECATION_WARNINGS = False

# Set ERR.txt as default error stream
f = open('ERR.txt', 'w')
sys.stderr = f

# Set tensorflow logging level to only print fatal errors
logging.getLogger('tensorflow').setLevel(logging.FATAL) 

# ------------------------------- Imports ------------------------------- #

from rdkit import Chem
import numpy as np
import random

# ------------------------------- Code ------------------------------- #

# Don't remember why but the if statement in the following row is needed,
# due to molecule_generation package usage
if __name__ == '__main__':

    # Specify the directory of the trained model
    model_dir = "..\\..\\model"

    # Specify the dataset path
    path1 = "..\\..\\dataset\\Commercial_MW\Commercial_MWlower330(changed).csv"
    path2 = "..\\..\\dataset\\Commercial_MW\Commercial_MW330-500(changed)1.csv"
    path3 = "..\\..\\dataset\\Commercial_MW\Commercial_MW330-500(changed)2.csv"
    path4 = "..\\..\\dataset\\Commercial_MW\Commercial_MWhigher500(changed).csv"

    # Open the files to read
    f1 = open(path1,"r")
    f2 = open(path2,"r")
    f3 = open(path3,"r")
    f4 = open(path4,"r")

    # List of input smiles strings
    input_smiles = []

    # Read all the smiles
    print("Let'start reading from files...")

    for x in f1:
        input_smiles.append(x.replace("\n",""))
    for x in f2:
        input_smiles.append(x.replace("\n",""))
    for x in f3:
        input_smiles.append(x.replace("\n",""))
    for x in f4:
        input_smiles.append(x.replace("\n",""))

    # Close the file
    f1.close()
    f2.close()
    f3.close()
    f4.close()

    # New line
    print("Here we go!")
    print()
    
    # Process latent vector for each input smiles string
    # Embeddings are not actually needed, just seeing if any error is triggered
    for i in range(len(input_smiles)):
        # Print index to manually observe at what index the error occurs
        print("{}% -> {} {}".format(round((i+1)/len(input_smiles)*100,2), i+1, input_smiles[i]))
        # Call to encode()
        try:
            Chem.MolFromSmiles(input_smiles[i])
        except Exception as e:
            # Empty ERR.txt
            f.truncate()

    # Empty ERR.txt
    f.truncate()
