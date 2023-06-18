# --------------------- Solve errors and warnings ---------------------- #

#     Currently this script produces a whole wall of errors and warnings
#     because molecule_generation module uses many deprecated functions.
#     Although, these all do not affect the correct functioning of
#     the script and will be fixed by Microsoft team before the deprecated
#     functions will have been deleted.

from math import sqrt
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
distributions = []

# ------------------------------ addNoise ------------------------------- #

def addNoise(embedding, area, num_of_molecules, num=None):

    # If the number of elements to add disturb to is not explicitly defined,
    # then add noise to the whole embedding
    if num == None:
        num = len(embedding)

    global interval_size
    global magnitude
    global distributions

    # For each feature calculates its range for noise
    noise_intervals_lengths = []
    for i in range(len(distributions)):
        current_feature_distribution = distributions[i]
        current_interval = round(embedding[i] - embedding[i]%interval_size,magnitude)
        
        # Fetures in distribution tails are not touched
        intervals = sorted(current_feature_distribution.keys())
        left_tail = 0
        count = 0
        for j in range(len(intervals)):
            if count < area/2:
                count += current_feature_distribution[intervals[j]]
            else:
                left_tail = intervals[j]
                break
        right_tail = 0
        count = 0
        for j in range(len(intervals)):
            if count < area/2:
                count += current_feature_distribution[intervals[len(intervals)-1-j]]
            else:
                right_tail = intervals[j]
                break
        if current_interval < left_tail or current_interval > right_tail:
            noise_intervals_lengths.append(0.0)
            continue

        # Find the interval in which noise must be added
        if embedding[i]%interval_size < interval_size-embedding[i]%interval_size:
            even = embedding[i]%interval_size
            odd = interval_size-embedding[i]%interval_size
            is_near_lower = True
        else:
            even = interval_size-embedding[i]%interval_size
            odd = embedding[i]%interval_size
            is_near_lower = False

        lower_bound = embedding[i]
        upper_bound = embedding[i]
        current_area = 0
        j = 0
        while True:
            if j%2 == 0:
                p_lower = current_feature_distribution[round(current_interval-interval_size*j/2,magnitude)]
                p_upper = current_feature_distribution[round(current_interval+interval_size*j/2,magnitude)]
                if current_area + even/interval_size*(p_lower + p_upper) >= area:
                    break
                lower_bound -= even
                upper_bound += even
                current_area += even/interval_size*(p_lower + p_upper)
            else:
                if is_near_lower:
                    p_lower = current_feature_distribution[round(current_interval-interval_size*(1+(j-j%2)/2),magnitude)]
                    p_upper = current_feature_distribution[round(current_interval+interval_size*(j-j%2)/2,magnitude)]
                else:
                    p_lower = current_feature_distribution[round(current_interval-interval_size*(j-j%2)/2,magnitude)]
                    p_upper = current_feature_distribution[round(current_interval+interval_size*(1+(j-j%2)/2),magnitude)]
                if current_area + odd/interval_size*(p_lower + p_upper) >= area:
                    break
                lower_bound -= odd
                upper_bound += odd
                current_area += odd/interval_size*(p_lower + p_upper)
            j += 1
        
        r = sqrt((area-current_area)*interval_size/(p_lower+p_upper))
        noise_intervals_lengths.append(upper_bound-lower_bound+2*r)
        
    modified_embeddings = []
    all_indexes = [j for j in range(len(embedding))]
    # Set the seed for the np.random.normal(...) function
    np.random.seed(int(time.time()))
    for i in range(num_of_molecules):

        # Calculate the indexes of the num elements to be modified
        random.shuffle(all_indexes)
        indexes_to_be_modified = all_indexes[0:num]

        # Initialize embedding to be returned
        modified_embedding = copy.deepcopy(embedding)

        for j in range(num):
            current_range = noise_intervals_lengths[indexes_to_be_modified[j]]
            modified_embedding[indexes_to_be_modified[j]] += random.random()*current_range-current_range/2
        modified_embeddings.append(modified_embedding)

    # Return the embeddings with noise addition
    return modified_embeddings

# ------------------------------ readDistributions ------------------------------- #

def readDistributions():
    # Specify the paths for the database
    path = "..\\features_distribution_plots\\features_distributions.txt"

    # Open the file
    f = open(path,"r")

    # Read the size of the discretization range
    global interval_size
    global magnitude
    global distributions
    interval_size = float(f.readline().strip())
    while interval_size*(10**magnitude) < 1:
        magnitude += 1

    # Read the ranges for each feature and put them in the rows map
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
            distributions.append(temp_map)
                
    # Close the file
    f.close()
    return


# ------------------------------- Code ------------------------------- #

# Don't remember why but the if statement in the following row is needed,
# due to molecule_generation package usage
if __name__ == '__main__':

    # Read from distributions file
    readDistributions()

    # Specify the directory of the trained model
    model_dir = "..\\..\\model"

    # List of input smiles strings
    input_smiles = ["COc1ccc(OS(=O)(=O)C(F)(F)F)c(Br)c1"]

    print(f"Encoded: {input_smiles}","\n")

    # Using the model at model_dir path
    with load_model_from_directory(model_dir) as model:
        print()

        # Number of molecules to be generated from each input smiles string
        num_of_molecules_to_generate = 10

        print("Start encoding...")
        # Process latent vector for each input smiles string
        embeddings = model.encode(input_smiles)

        # Calculate num_of_molecules_to_generate diverse molecules for each input smiles string,
        # adding noise to the embeddings
        list_of_embeddings = []
        print("Start adding noise...")
        list_of_embeddings = addNoise(embeddings[0],0.2,num_of_molecules_to_generate)

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