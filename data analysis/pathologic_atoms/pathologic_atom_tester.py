# ------------------------------- Imports ------------------------------- #

from molecule_generation import load_model_from_directory
import numpy as np
import random
import sys
import io
from rdkit import Chem

# ------------------------------- Code ------------------------------- #

# Don't remember why but the if statement in the following row is needed,
# due to molecule_generation package usage
if __name__ == '__main__':

    # Specify the paths for the database
    path1 = "..\\..\\dataset\\Commercial_MW\Commercial_MWlower330.csv"
    path2 = "..\\..\\dataset\\Commercial_MW\Commercial_MW330-500-1.csv"
    path3 = "..\\..\\dataset\\Commercial_MW\Commercial_MW330-500-2.csv"
    path4 = "..\\..\\dataset\\Commercial_MW\Commercial_MWhigher500.csv"

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

    chosen_smiles = []

    # Take those which contain [15N]
    for i in range(len(input_smiles)):
        if "[15N]" in input_smiles[i]:
            chosen_smiles.append(input_smiles[i])

    print("Nel dataset, {} contengono \"[15N]\" \n".format(len(chosen_smiles)))
    print("Ora proviamo a vedere se ci sono ancora errori sostituendo \"[15N]\" con \"[15N+]\" \n")

    changed_smiles = []

    # Compute changed
    print("I'm computing changed")
    for i in range(len(chosen_smiles)):
        changed_smiles.append(chosen_smiles[i].replace("[15N]","[15N+]"))     

    # Print chosen smiles
    print("################ Print Chosen #####################\n")
    for i in range(len(chosen_smiles)):
        # CHOSEN
        print(chosen_smiles[i])
        Chem.MolFromSmiles(chosen_smiles[i])
        print()

    print("################ Print Changed #####################\n")
    for i in range(len(changed_smiles)):
        # CHANGED
        print(changed_smiles[i])
        Chem.MolFromSmiles(changed_smiles[i])
        print()


