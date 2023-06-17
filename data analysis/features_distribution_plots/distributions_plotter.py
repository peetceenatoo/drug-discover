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

# ------------------------------- Code ------------------------------- #

# Don't remember why but the if statement in the following row is needed,
# due to molecule_generation package usage
if __name__ == '__main__':

    # Specify the directory of the trained model
    model_dir = "model"

    # Specify the paths for the database
    path1 = "..\\dataset\\Commercial_MW\Commercial_MWlower330(clean).csv"
    path2 = "..\\dataset\\Commercial_MW\Commercial_MW330-500(clean)1.csv"
    path3 = "..\\dataset\\Commercial_MW\Commercial_MW330-500(clean)2.csv"
    path4 = "..\\dataset\\Commercial_MW\Commercial_MWhigher500(clean).csv"
    out1 = "chosen_molecules.csv"
    out2 = "features_max_min.csv"

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

    # Pick num_of_molecules_to_plot random smiles
    num_of_molecules_to_plot = 100
    random.shuffle(input_smiles)
    chosen_smiles = input_smiles[0:num_of_molecules_to_plot]

    # Open the output file to store all the num_of_molecules_to_plot smiles
    fout = open(out1,"w")

    # Write on the output file
    for c in chosen_smiles:
        fout.write("{}".format(c))

    # Close the output file
    fout.close()

    # Using the model at model_dir path
    with load_model_from_directory(model_dir) as model:
        print()

        # Process latent vector for each input smiles string
        for i in range(len(chosen_smiles)):
            embeddings = model.encode([chosen_smiles[i]])

        # Set number of embedding features
        num_of_features = len(embeddings[0])

        # Init x-axis data
        lower = -1
        upper = 1
        jump = 0.05
        # Magnitude must equal reciprocal( min( magnitude(jump), magnitude(lower), magnitude(upper) ) )4
        # Needed because range() only works with integer parameters
        magnitude = 100

        # Define x-axis
        x = list( range( (int)(lower*magnitude), (int)(upper*magnitude+1), (int)(jump*magnitude) ) )
        for i in range(len(x)):
            x[i] = x[i]/magnitude
        x.insert(0, lower - jump/2)
        x.append(upper + jump/2)

        # Write the names of the CSV file columns
        fout = open(out2,"w")
        fout.write("feature,max,min\n")

        # For every feature (from 0 to num_of_features-1)
        for i in range(num_of_features):

            # Init y-axis as list of len(x) elements
            y = [0]*len(x)
            # Init min and max
            min = embeddings[0][i]
            max = embeddings[0][i]

            # For each smiles, eventually save the current smiles (j-th) feature as min or max
            for j in range(num_of_molecules_to_plot):

                # Update max and min
                if embeddings[j][i] < min:
                    min = embeddings[j][i]
                if embeddings[j][i] > max:
                    max = embeddings[j][i]

                # If the current feature for the current smiles is to the left of the range on the x-axis
                if embeddings[j][i] < x[0]:
                    # Plot one more molecule on the left bond
                    y[0] += 1
                # If the current feature for the current smiles is to the right of the range on the x-axis
                elif embeddings[j][i] >= x[len(x)-1]:
                    # Plot one more molecule on the right bond
                    y[len(x)-1] += 1
                # Otherwise, plot one more molecule at the right value on the x-asis
                else:
                    # Calculate the value on the x-axis
                    k=1
                    while embeddings[j][i] >= x[k]-jump/2 and k<len(x)-1:
                        k += 1
                    # Plot one more
                    y[k-1] += 1 

            # Write max and min for the current feature (i-th)
            fout.write("{},{},{}\n".format(i,max,min))

            # Plot
            plt.plot(x,y)
            plt.xlabel('feature value')
            plt.ylabel('feature count')
            plt.title("feature n.{}".format(i))
            plt.savefig("features_distribution_plots\\feature{}.png".format(i))
            plt.clf()

        # Close the output file
        fout.close()
