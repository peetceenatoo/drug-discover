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
    path1 = "..\\..\\dataset\\Commercial_MW\Commercial_MWlower330(clean).csv"
    path2 = "..\\..\\dataset\\Commercial_MW\Commercial_MW330-500(clean)1.csv"
    path3 = "..\\..\\dataset\\Commercial_MW\Commercial_MW330-500(clean)2.csv"
    path4 = "..\\..\\dataset\\Commercial_MW\Commercial_MWhigher500(clean).csv"

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

    print("Nel dataset ci sono {} molecole e ".format(len(input_smiles)))

    chosen_smiles_temp = []
    chosen_smiles = []

    # Take those which contain [O]
    for i in range(len(input_smiles)):
        if "[O]" in input_smiles[i]:
            chosen_smiles_temp.append(input_smiles[i])

    # Now only take the ones which contain a pattern
    for i in range(len(chosen_smiles_temp)):
        if "N([O])" in chosen_smiles_temp[i] or "N[O]" in chosen_smiles_temp[i] or "N1[O]" in chosen_smiles_temp[i] or "N2[O]" in chosen_smiles_temp[i] or "N3[O]" in chosen_smiles_temp[i] or "N4[O]" in chosen_smiles_temp[i] or "N5[O]" in chosen_smiles_temp[i]:
            chosen_smiles.append(chosen_smiles_temp[i])

    print("di queste {} contengono \"[O]\" \n".format(len(chosen_smiles)))
    print("Ora proviamo a vedere se ci sono ancora errori sostituendo \"[O]\" con \"[O+]\" \n")

    changed_smiles = []

    # Compute changed
    print("I'm computing changed")
    for i in range(len(chosen_smiles)):
        changed_smiles.append(chosen_smiles[i].replace("[O]","[O+]"))     

    # Print chosen smiles
    print("################ Print Chosen #####################\n")
    for i in range(len(chosen_smiles)):
        # CHOSEN
        print("{}% ".format(i/len(changed_smiles)*100))
        print(chosen_smiles[i])
        Chem.MolFromSmiles(chosen_smiles[i])
        print()

    print("################ Print Changed #####################\n")
    for i in range(len(changed_smiles)):
        # CHANGED
        print("{}% ".format(i/len(changed_smiles)*100))
        print(changed_smiles[i])
        Chem.MolFromSmiles(changed_smiles[i])
        print()


