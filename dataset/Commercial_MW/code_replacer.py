# ------------------------------- Imports ------------------------------- #
from rdkit import Chem

# ------------------------------- Code ------------------------------- #

# Don't remember why but the if statement in the following row is needed,
# due to molecule_generation package usage
if __name__ == '__main__':

    # Specify the paths for the database
    path1 = "Commercial_MWlower330.csv"
    path2 = "Commercial_MW330-500-1.csv"
    path3 = "Commercial_MW330-500-2.csv"
    path4 = "Commercial_MWhigher500.csv"

    # Open the input files
    f1 = open(path1,"r")
    f2 = open(path2,"r")
    f3 = open(path3,"r")
    f4 = open(path4,"r")

    # List of input smiles strings
    input_smiles1 = []
    input_smiles2 = []
    input_smiles3 = []
    input_smiles4 = []

    # Read all the smiles
    for x in f1:
        input_smiles1.append(x.replace("\n","").replace("[BH3]","[BH3-]"))
    for x in f2:
        input_smiles2.append(x.replace("\n","").replace("[BH3]","[BH3-]"))
    for x in f3:
        input_smiles3.append(x.replace("\n","").replace("[BH3]","[BH3-]"))
    for x in f4:
        input_smiles4.append(x.replace("\n","").replace("[BH3]","[BH3-]"))

    # Close the input files
    f1.close()
    f2.close()
    f3.close()
    f4.close()

    # Open the output files
    f1 = open(path1,"w")
    f2 = open(path2,"w")
    f3 = open(path3,"w")
    f4 = open(path4,"w")

    # Write on the output files
    for x in input_smiles1:
        f1.write(x+"\n")
    for x in input_smiles2:
        f2.write(x+"\n")
    for x in input_smiles3:
        f3.write(x+"\n")
    for x in input_smiles4:
        f4.write(x+"\n")

    f1.close()
    f2.close()
    f3.close()
    f4.close()

