# --------------------- Hide errors and warnings ---------------------- #

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

from molecule_generation import load_model_from_directory
import random
import numpy as np
import matplotlib.pyplot as plt

# ------------------------------- Code ------------------------------- #

# Don't remember why but the if statement in the following row is needed,
# due to molecule_generation package usage
if __name__ == '__main__':

    # Specify the directory of the trained model
    model_dir = "model"

    #specify the paths for the database
    path1 = "dataset\\Commercial_MW\Commercial_MWlower330(clean).csv"
    path2 = "dataset\\Commercial_MW\Commercial_MW330-500(clean)1.csv"
    path3 = "dataset\\Commercial_MW\Commercial_MW330-500(clean)2.csv"
    path4 = "dataset\\Commercial_MW\Commercial_MWhigher500(clean).csv"

    # List of input smiles strings
    input_smiles = []

    #open the files
    f1 = open(path1,"r")
    f2 = open(path2,"r")
    f3 = open(path3,"r")
    f4 = open(path4,"r")

    #take all the smiles
    for x in f1:
        input_smiles.append(x)

    for x in f2:
        input_smiles.append(x)

    for x in f3:
        input_smiles.append(x)

    for x in f4:
        input_smiles.append(x)

    #close the files
    f1.close()
    f2.close()
    f3.close()
    f4.close()

    #pick random smiles
    num_of_molecules_to_plot = 100
    random.shuffle(input_smiles)
    chosen_smiles = input_smiles[0:num_of_molecules_to_plot]

    # Using the model at model_dir path
    with load_model_from_directory(model_dir) as model:
        print()

        # Process latent vector for each input smiles string
        embeddings = model.encode(chosen_smiles)

        # Calculate num_of_molecules_to_generate diverse molecules for each input smiles string,
        # adding noise to the embeddings
        num_of_features = len(embeddings[0])

        for i in range(num_of_features):
            x = []
            y = []
            for j in range(num_of_molecules_to_plot):
                y.append(embeddings[j][i])
                x.append(j)
            plt.plot(x,y)
            plt.xlabel('chosen_molecules')
            plt.ylabel('features')
            plt.title("feature n.{}".format(i))
            plt.show()  

    # Empty ERR.txt
    f.truncate()
