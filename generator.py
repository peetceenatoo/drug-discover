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
import numpy as np
import random
import time

# ------------------------------ Functions ------------------------------ #

def addNoise(embedding, scale=1, num=None):

    # If the number of elements to add disturb to is not explicitly defined,
    # then add noise to the whole embedding
    if num == None:
        num = len(embedding)

    # Calculate the indexes of the num elements to be modified
    all_indexes = [i for i in range(len(embedding))]
    random.shuffle(all_indexes)
    indexes_to_be_modified = all_indexes[0:num]

    # Initialize embedding to be returned
    modified_embedding = embedding.copy()

    # Set the seed for the np.random.normal(...) function
    np.random.seed(int(time.time()))

    # For each index previously calculated
    for i in range(len(indexes_to_be_modified)):

        # Calculate the magnitude of the element at i-th index
        magnitude = 1.0
        while abs(modified_embedding[indexes_to_be_modified[i]]*magnitude)<1.0 :
            magnitude = magnitude*10

        # Add noise to such element
        modified_embedding[indexes_to_be_modified[i]] += np.random.normal(0, scale/magnitude)

        # --- About how to choose the scale for np.random.normal(...) function --- 
        #   Using num = len(embedding): we noticed that much lower values than scale/magnitude
        #                               are likely to provide many duplicates. On the other hand,
        #                               much bigger values than scale/magnitude are likely to
        #                               provide molecules which are too different from the input.

    # Return the embedding with noise addition
    return modified_embedding

# ------------------------------- Code ------------------------------- #

# Don't remember why but the if statement in the following row is needed,
# due to molecule_generation package usage
if __name__ == '__main__':

    # Specify the directory of the trained model
    model_dir = "model"

    # List of input smiles strings
    input_smiles = ["CCC(C)NC(=O)COc1cccc(C)c1"]
    print(f"Encoded: {input_smiles}","\n")

    # Using the model at model_dir path
    with load_model_from_directory(model_dir) as model:
        print()

        # Number of molecules to be generated from each input smiles string
        num_of_molecules_to_generate = 20

        # Process latent vector for each input smiles string
        embeddings = model.encode(input_smiles)
        
        # Calculate num_of_molecules_to_generate diverse molecules for each input smiles string,
        # adding noise to the embeddings
        list_of_embeddings = []
        for i in range(num_of_molecules_to_generate):
            list_of_embeddings.append(addNoise(embeddings[0], 1))

        # Decode without a scaffold constraint
        list_of_decoded_smiles = model.decode(list_of_embeddings)

        # Print all generated smiles
        for i in range(num_of_molecules_to_generate):
            print("Decoded molecule n.{}: ".format(i+1), list_of_decoded_smiles[i],"\n")
        
        print("Have a nice day :)")
