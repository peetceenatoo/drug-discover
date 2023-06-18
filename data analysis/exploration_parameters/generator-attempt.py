# --------------------- Solve errors and warnings ---------------------- #

#     Currently this script produces a whole wall of errors and warnings
#     because molecule_generation module uses many deprecated functions.
#     Although, these all do not affect the correct functioning of
#     the script and will be fixed by Microsoft team before the deprecated
#     functions will have been deleted.

import os
import logging
from tensorflow.python.util import deprecation

# Disable logging output of tensorflow content [May be useless] 
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'

# Do not print deprecation warnings of tensorflow content
deprecation._PRINT_DEPRECATION_WARNINGS = False

# Set tensorflow logging level to only print fatal errors
logging.getLogger('tensorflow').setLevel(logging.FATAL) 

# ----------------------------------------------------------------------- #


# ------------------------------- Imports ------------------------------- #

from molecule_generation import load_model_from_directory
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import DataStructs
import numpy as np
import os
import random
import time
import copy

# ------------------------------ Global ------------------------------- #

# Compute module of the magnitude of the size of the intervals
interval_size = 0
magnitude = 0

# ------------------------------ addNoise ------------------------------- #

def addNoise(distributions, area, embedding, num=None):

    # If the number of elements to add disturb to is not explicitly defined,
    # then add noise to the whole embedding
    if num == None:
        num = len(embedding)

    # Calculate the indexes of the num elements to be modified
    all_indexes = [i for i in range(len(embedding))]
    random.shuffle(all_indexes)
    indexes_to_be_modified = all_indexes[0:num]

    # Initialize embedding to be returned
    modified_embedding = copy.deepcopy(embedding)

    # Set the seed for the np.random.normal(...) function
    np.random.seed(int(time.time()))

    # For each index previously calculated
    for i in range(len(indexes_to_be_modified)):
        current_feature = indexes_to_be_modified[i]
        current_feature_distribution = distributions[current_feature]
        current_interval = round(embedding[current_feature] - embedding[current_feature]%interval_size,magnitude)
        intervals = sorted(current_feature_distribution.keys())
        
        if current_interval < intervals[0] or current_interval > intervals[len(intervals)-1]:
            continue
        left_tail = 0
        right_tail = 0
        for j in range(len(intervals)):
            if left_tail < area/2:
                left_tail += intervals[j]
            if right_tail < 1 - area/2:
                right_tail += intervals[j]
        if current_interval < left_tail or current_interval > right_tail:
            continue

        base = 0
        if current_interval
        # Add noise to such element
        modified_embedding[indexes_to_be_modified[i]] += np.random.normal(0, scale/magnitude)

        # --- About how to choose the scale for np.random.normal(...) function --- 
        #   Using num = len(embedding): we noticed that much lower values than scale/magnitude
        #                               are likely to provide many duplicates. On the other hand,
        #                               much bigger values than scale/magnitude are likely to
        #                               provide molecules which are too different from the input.

    # Return the embedding with noise addition
    return modified_embedding

# ------------------------------ readDistributions ------------------------------- #

def readDistributions():
    # Specify the paths for the database
    path = "..\\features_distribution_plots\\features_distributions.txt"

    # Open the file
    f = open(path,"r")

    # Read the size of the discretization range
    interval_size = float(f.readline().strip())
    while interval_size*(10**magnitude) < 1:
        magnitude+=1

    # Read the ranges for each feature and put them in the rows map
    rows_list = []
    for row in f:
            # Remove leading/trailing whitespace and newline characters
            row = row.strip()
            
            # Split the row by colon
            parts = row.split(":")
            bottom = float(parts[0])
            probs_part = parts[1]
            
            # Split the second part by semicolon
            probs = probs_part.split(";")
            probs = probs[:len(probs)-1]

            # Compute the map of probabilities for the current row
            temp_map = {}
            for i in range(len(probs)):
                temp_map[round(bottom+i*interval_size,magnitude)] = float(probs[i])
                
            # Put current map of probabilities in the rows map
            rows_list.append(temp_map)
                
    # Close the file
    f.close()

    # Return the list
    return rows_list, interval_size


# ------------------------------- Code ------------------------------- #

# Don't remember why but the if statement in the following row is needed,
# due to molecule_generation package usage
if __name__ == '__main__':

    # Read from distributions file
    readDistributions()

    # Specify the directory of the trained model
    model_dir = "..\\..\\model"

    # List of input smiles strings
    input_smiles = ["CCOc1cc(ccc1OC)c2nnc(SCC(=O)Nc3cccc(Cl)c3)nc2c4ccc(OC)c(OCC)c4"]

    print(f"Encoded: {input_smiles}","\n")

    features_distributions, interval_size = readDistributions()

    # Using the model at model_dir path
    with load_model_from_directory(model_dir) as model:
        print()

        # Number of molecules to be generated from each input smiles string
        num_of_molecules_to_generate = 1

        print("Start encoding...")
        # Process latent vector for each input smiles string
        embeddings = model.encode(input_smiles)

        # Calculate num_of_molecules_to_generate diverse molecules for each input smiles string,
        # adding noise to the embeddings
        list_of_embeddings = []
        print("Start adding noise...")
        for i in range(num_of_molecules_to_generate):
            list_of_embeddings.append(addNoise(???))

        # Decode without a scaffold constraint
        print("Start decoding...")
        list_of_decoded_smiles = model.decode(list_of_embeddings)

    # Compute fingerprints
    print("Generating fingerprints...")
    in_fps = [AllChem.RDKFingerprint(Chem.MolFromSmiles(x)) for x in input_smiles]
    out_fps = [AllChem.RDKFingerprint(Chem.MolFromSmiles(x)) for x in list_of_decoded_smiles]

    # Calculate similarity
    for i in range(len(list_of_decoded_smiles)):
        print("Input molecule: {}\nOutput molecule: {}\nSimilarity: {}\n".format(input_smiles[0], list_of_decoded_smiles[i], DataStructs.TanimotoSimilarity(in_fps[0], out_fps[i])))

    print("Have a nice day :)")