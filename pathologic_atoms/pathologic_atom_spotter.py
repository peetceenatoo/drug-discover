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

import matplotlib.pyplot as plt
from molecule_generation import load_model_from_directory
import numpy as np
import random

# ------------------------------- Code ------------------------------- #

# Don't remember why but the if statement in the following row is needed,
# due to molecule_generation package usage
if __name__ == '__main__':

    # Specify the directory of the trained model
    model_dir = "..\\model"

    # Specify the dataset path
    path = "..\\dataset\\Commercial_MW\Commercial_MWhigher500(clean).csv"
    out1 = "chosen_molecules.csv"

    # Open the dataset file
    f = open(path,"r")

    # List of input smiles strings
    input_smiles = []

    # Read all the smiles
    for x in f:
        input_smiles.append(x.replace("\n",""))

    # Close the file
    f.close()

    # Pick num_of_molecules_to_plot random smiles
    num_of_molecules_to_plot = 100
    random.shuffle(input_smiles)
    chosen_smiles = input_smiles[0:num_of_molecules_to_plot]

    # Open the output file to store all the num_of_molecules_to_plot smiles
    fout = open(out1,"w")
    # Write the names of the CSV file columns
    fout.write("Entry,Smiles\n")

    # For each molecule in chosen_smiles, write index (Entry) and Smiles
    for i in range(num_of_molecules_to_plot):
        fout.write("{},{}\n".format(i,chosen_smiles[i]))

    # Close the output file
    fout.close()

    # Using the model at model_dir path
    with load_model_from_directory(model_dir) as model:
        print()

        # Process latent vector for each input smiles string
        # Embeddings are not actually needed, just seeing if any error is triggered
        for i in range(num_of_molecules_to_plot):
            # Print index to manually observe at what index the error occurs
            print(i)
            # Call to encode()
            try:
                embeddings = model.encode([chosen_smiles[i]])
            except Exception as e:
                # Empty ERR.txt
                f.truncate()

    # Empty ERR.txt
    f.truncate()
