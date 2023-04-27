# ------------------ Hide errors and warnings --------------------- #

import array
import os
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'
from tensorflow.python.util import deprecation
deprecation._PRINT_DEPRECATION_WARNINGS = False
import sys
f = open('ERR.txt', 'w')
sys.stderr = f
from contextlib import redirect_stdout
import logging
logging.getLogger('tensorflow').setLevel(logging.ERROR) 

# ------------------------------- Import ------------------------------- #

from molecule_generation import load_model_from_directory
import random
import numpy as np
import time

# ---------------------------- Functions ---------------------------- #

def disturb(array, scale=1, num=None):
    #np.random.seed(random.randint(-100000000,100000000))
    if num == None:
        num = len(array)
    numbers = [i for i in range(len(array))]
    random.shuffle(numbers)
    indexes = numbers[0:num]
    temp_embedding = array.copy()
    np.random.seed(int(time.time()))
    for i in range(num):
        #taking magnitude of the element
        magnitude = 1.0
        while abs(temp_embedding[indexes[i]]*magnitude)<1.0 :
            magnitude = magnitude*10
        temp_embedding[indexes[i]] += np.random.normal(0,scale/magnitude*10)  
    return temp_embedding

# ------------------------------- Code ------------------------------- #

# Don't remember why but the if statement in the following row is needed
if __name__ == '__main__':
    model_dir = "model"
    example_smiles = ["CCC(C)NC(=O)COc1cccc(C)c1"]
    print(f"Encoded: {example_smiles}")

    with load_model_from_directory(model_dir) as model:
        num_of_generated_molecules = 20

        # Process latent vector
        embeddings = model.encode(example_smiles)

        # Generate similar embeddings
        list_of_embeddings = []
        for i in range(num_of_generated_molecules):
            list_of_embeddings.append(disturb(embeddings[0],1))

        # Decode without a scaffold constraint.
        # DO NOT PRINT ANYTHING HERE !!
        with redirect_stdout(f):
            list_of_decoded_smiles = model.decode(list_of_embeddings)

        for i in range(num_of_generated_molecules):
            print("Decoded molecule n.{}: ".format(i+1), list_of_decoded_smiles[i])
