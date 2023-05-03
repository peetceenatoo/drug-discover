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
    out1 = "features_distribution_plots\chosen_molecules.csv"
    out2 = "features_distribution_plots\\features_max_min.csv"

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

    #create the output file
    fout1 = open(out1,"w")

    for c in chosen_smiles:
        fout1.write("{}".format(c))

    fout1.close()

    # Using the model at model_dir path
    with load_model_from_directory(model_dir) as model:
        print()

        # Process latent vector for each input smiles string
        embeddings = model.encode(chosen_smiles)

        num_of_features = len(embeddings[0])

        #prepare the x-axis
        lower = -1
        upper = 1
        jump = 0.05
        magnitude = 100 #max of magnitudes

        x = list(range((int)(lower*magnitude),(int)(upper*magnitude+1),(int)(jump*magnitude)))
        for i in range(len(x)):
            x[i] = x[i]/magnitude
        x.insert(0,lower-jump/2)
        x.append(upper+jump/2)

        #plotting the features
        fout2 = open(out2,"w")
        fout2.write("feature,max,min\n")
        for i in range(num_of_features):
            y = [0]*len(x)
            min = embeddings[0][i]
            max = embeddings[0][i]
            for j in range(num_of_molecules_to_plot):
                #updating max and min
                if embeddings[j][i] < min:
                    min = embeddings[j][i]
                if embeddings[j][i] > max:
                    max = embeddings[j][i]
                #upate data for the plot
                if embeddings[j][i] < x[0]:
                    y[0] += 1
                elif embeddings[j][i] >= x[len(x)-1]:
                    y[len(x)-1] += 1
                else:
                    k=1
                    while embeddings[j][i] >= x[k]-jump/2 and k<len(x)-1:
                        k += 1
                    y[k-1] +=1 
            #update max and min in the file
            fout2.write("{},{},{}\n".format(i,max,min))
            #update the plot
            plt.plot(x,y)
            plt.xlabel('feature value')
            plt.ylabel('feature count')
            plt.title("feature n.{}".format(i))
            plt.savefig("features_distribution_plots\\feature{}.png".format(i))
            plt.clf()
        fout2.close()

    # Empty ERR.txt
    f.truncate()
