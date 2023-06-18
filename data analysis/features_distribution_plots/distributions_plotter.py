# --------------------- Solve errors and warnings ---------------------- #

import os
import logging
from tensorflow.python.util import deprecation

# Disable logging output of tensorflow content [May be useless]
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'

# Do not print deprecation warnings of tensorflow content
deprecation._PRINT_DEPRECATION_WARNINGS = False

# Set tensorflow logging level to only print fatal errors
logging.getLogger('tensorflow').setLevel(logging.FATAL)

# ------------------------------- Imports ------------------------------- #

import matplotlib.pyplot as plt
from molecule_generation import load_model_from_directory
import random
import copy

# ------------------------------- Code ------------------------------- #

# Don't remember why but the if statement in the following row is needed,
# due to molecule_generation package usage
if __name__ == '__main__':

    # Define the number of molecule to analyse
<<<<<<< Updated upstream
    num_of_chosen_molecules = 100
=======
<<<<<<< HEAD
    num_of_chosen_molecules = 10000
=======
    num_of_chosen_molecules = 100
>>>>>>> 0d4704de37a845e7a271a1e4e3c05a08a80d77a5
>>>>>>> Stashed changes
    # Define the interval length
    interval = 0.02

    # Specify the directory of the trained model
    model_dir = "..\\..\\model"

    # Specify the paths for the database
    path1 = "..\\..\\dataset\\Commercial_MW\Commercial_MWlower330.csv"
    path2 = "..\\..\\dataset\\Commercial_MW\Commercial_MW330-500-1.csv"
    path3 = "..\\..\\dataset\\Commercial_MW\Commercial_MW330-500-2.csv"
    path4 = "..\\..\\dataset\\Commercial_MW\Commercial_MWhigher500.csv"
    out1 = "features_distributions.txt"

    # Open the files
    f1 = open(path1,"r")
    f2 = open(path2,"r")
    f3 = open(path3,"r")
    f4 = open(path4,"r")

    # List of input smiles strings
    input_smiles = []

    # Read all the smiles
    for x in f1:
        input_smiles.append(x.replace("\n",""))
    for x in f2:
        input_smiles.append(x.replace("\n",""))
    for x in f3:
        input_smiles.append(x.replace("\n",""))
    for x in f4:
        input_smiles.append(x.replace("\n",""))

    # Close the files
    f1.close()
    f2.close()
    f3.close()
    f4.close()

    # Pick num_of_chosen_molecules random smiles
    random.shuffle(input_smiles)
    chosen_smiles = input_smiles[0:num_of_chosen_molecules]

    # Using the model at model_dir path
    with load_model_from_directory(model_dir) as model:
        print()

        # Process latent vector for each input smiles string
        embeddings = model.encode(chosen_smiles)

    # Set number of embedding features
    num_of_features = len(embeddings[0])

    # Lists of min and max for all features
    min = copy.deepcopy(embeddings[0])
    max = copy.deepcopy(embeddings[0])

    for i in range(1,num_of_chosen_molecules):
        for j in range(num_of_features):
            if embeddings[i][j] < min[j]:
                min[j] = embeddings[i][j]
            if embeddings[i][j] > max[j]:
                max[j] = embeddings[i][j]

    # The module of the magnitude of interval + 1
    magnitude = 0
    factor = 1
    while interval*factor < 10 :
        magnitude += 1
        factor = factor * 10

    # Open the output file
    fout = open(out1,"w")
    fout.write("{}\n".format(interval))

    # For each feature
    for j in range(num_of_features):

        # Init x-axis data
        lower = round(min[j] - min[j]%interval + interval/2, magnitude)
        upper = round(max[j] - max[j]%interval + interval/2, magnitude)
        fout.write("{}:".format(lower - interval/2))

        # Define x-axis
        x = list( range( (int)(lower*factor), (int)(upper*factor+1), (int)(interval*factor) ) )
        for k in range(len(x)):
            x[k] = x[k]/factor

        # Init y-axis as list of len(x) elements
        y = [0]*len(x)

        for i in range(num_of_chosen_molecules):

            # Calculate the value on the y-axis
            k=0
            while k<len(y) and embeddings[i][j] >= x[k]-interval/2:
                k += 1
            y[k-1] += 1

        # Set the data in percentage
        for k in range(len(y)):
            y[k] = y[k]/num_of_chosen_molecules
            fout.write("{};".format(y[k]))
            y[k] = y[k]*100
        fout.write("\n")

        # Plot
        plt.plot(x,y)
        plt.xlabel('feature values')
        plt.ylabel('feature counts in percentage')
        plt.title("feature n.{}".format(j+1))
        plt.savefig("feature{}.png".format(j+1))
        plt.clf()

    fout.close()
