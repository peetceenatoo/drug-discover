# --------------------- Solve errors and warnings ---------------------- #


#     Currently this script produces a whole wall of errors and warnings
#     because molecule_generation module uses many deprecated functions.
#     Although, these all do not affect the correct functioning of
#     the script and will be fixed by Microsoft team before the deprecated
#     functions will have been deleted.

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

# ----------------------------------------------------------------------- #


# ------------------------------- Imports ------------------------------- #

from molecule_generation import load_model_from_directory
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import DataStructs
import random
# ------------------------------- Code ------------------------------- #

# Don't remember why but the if statement in the following row is needed,
# due to molecule_generation package usage
if __name__ == '__main__':

    # Specify the various paths
    model_dir = "../../model"

    path1 = "../../dataset/Commercial_MW/Commercial_MWlower330.csv"
    path2 = "../../dataset/Commercial_MW/Commercial_MW330-500-1.csv"
    path3 = "../../dataset/Commercial_MW/Commercial_MW330-500-2.csv"
    path4 = "../../dataset/Commercial_MW/Commercial_MWhigher500.csv"
    out_path = "input-output_similarity.csv"

    # List of input smiles strings
    input_smiles = []

    # Open the files
    f1 = open(path1,"r")
    f2 = open(path2,"r")
    f3 = open(path3,"r")
    f4 = open(path4,"r")

    print("Start reading all the smiles...")

    # Read all the smiles from f1
    for x in f1:
        input_smiles.append(x.replace("\n",""))

    # Read all the smiles from f2
    for x in f2:
        input_smiles.append(x.replace("\n",""))

    # Read all the smiles from f3
    for x in f3:
        input_smiles.append(x.replace("\n",""))

    # Read all the smiles from f4
    for x in f4:
        input_smiles.append(x.replace("\n",""))

    #close the files
    f1.close()
    f2.close()
    f3.close()
    f4.close()

    #taking only num random smiles
    num = 10
    random.shuffle(input_smiles)
    chosen_smiles_in = input_smiles[0:num]

    # Using the model at model_dir path
    with load_model_from_directory(model_dir) as model:
        print()

        print("Start encoding...")
        # Process latent vector for each input smiles string
        try:
            embeddings = model.encode(chosen_smiles_in)
        except Exception as e:
            # Empty ERR.txt
            f.truncate()

        # Decode without a scaffold constraint.
        # DO NOT PRINT ANYTHING HERE !!
        print("Start decoding...")
        with redirect_stdout(f):
            chosen_smiles_out = model.decode(embeddings)

    #generate fingerprints
    print("Generating fingerprints...")
    in_fps = [AllChem.RDKFingerprint(Chem.MolFromSmiles(x)) for x in chosen_smiles_in]
    out_fps = [AllChem.RDKFingerprint(Chem.MolFromSmiles(x)) for x in chosen_smiles_out]

    out = open(out_path,"w")
    out.write("Input Molecule,OutputMolecule,Tanimoto Distance\n")

    #calculate similarity
    print("Printing output...")
    for i in range(len(chosen_smiles_in)):
        out.write("{},{},{}\n".format(chosen_smiles_in[i],chosen_smiles_out[i],DataStructs.TanimotoSimilarity(in_fps[i],out_fps[i])))

    out.close()

    print("Have a nice day :)")

    # Empty ERR.txt
    f.truncate()
