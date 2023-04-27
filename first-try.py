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

# ---------------------------- Functions ---------------------------- #

def foo(array, scale, num=None):
    #np.random.seed(random.randint(-100000000,100000000))
    if num == None:
        num = len(array)
    numbers = [i for i in range(len(array))]
    random.shuffle(numbers)
    indexes = numbers[0:num]
    for i in range(num):
        array[indexes[i]] += np.random.normal(0,scale)
    print("Stampo un array dentro a foo...\n", array)
    return array

# ------------------------------- Code ------------------------------- #

# Don't remember why but the if statement in the following row is needed
if __name__ == '__main__':
    model_dir = "C:/Users/picci/AppData/Local/Programs/Python/Python310/Lib/site-packages/molecule_generation/tmp/MoLeR_checkpoint"
    example_smiles = ["CCC(C)NC(=O)COc1cccc(C)c1"]

    with load_model_from_directory(model_dir) as model:

        # Process latent vector
        embeddings = model.encode(example_smiles)
        print(f"Encoded: {example_smiles}")

        # Generate 2 similar embeddings
        list_of_embeddings = []
        for i in range(2):
            list_of_embeddings.append(foo(embeddings,0.01))
            print("Stampo un array fuori da foo...\n", list_of_embeddings[i])

        # Decode without a scaffold constraint.
        # DO NOT PRINT ANYTHING HERE !!
        with redirect_stdout(f):
            list_of_decoded_smiles = []
            for i in range(2):
                list_of_decoded_smiles.append(model.decode(list_of_embeddings[i]))

        for i in range(2):
            print("Decoded molecule n.{}: ".format(i+1), list_of_decoded_smiles[i])
        


    







