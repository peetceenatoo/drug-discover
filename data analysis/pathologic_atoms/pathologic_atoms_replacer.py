# ------------------------------- Imports ------------------------------- #
from rdkit import Chem

# ------------------------------- Code ------------------------------- #

# Specify the paths for the database
path1 = "..\\..\\dataset\\Commercial_MW\\Commercial_MWlower330(clean).csv"
path2 = "..\\..\\dataset\\Commercial_MW\\Commercial_MW330-500(clean)1.csv"
path3 = "..\\..\\dataset\\Commercial_MW\\Commercial_MW330-500(clean)2.csv"
path4 = "..\\..\\dataset\\Commercial_MW\\Commercial_MWhigher500(clean).csv"
path5 = "..\\..\\dataset\\Commercial_MW\\Commercial_MWlower330.csv"
path6 = "..\\..\\dataset\\Commercial_MW\\Commercial_MW330-500-1.csv"
path7 = "..\\..\\dataset\\Commercial_MW\\Commercial_MW330-500-2.csv"
path8 = "..\\..\\dataset\\Commercial_MW\\Commercial_MWhigher500.csv"

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
    input_smiles1.append(x.replace("\n","").replace("[N]","[N+]").replace("[O]","[O+]").replace("[B]","[B-]"))
for x in f2:
    input_smiles2.append(x.replace("\n","").replace("[N]","[N+]").replace("[O]","[O+]").replace("[B]","[B-]"))
for x in f3:
    input_smiles3.append(x.replace("\n","").replace("[N]","[N+]").replace("[O]","[O+]").replace("[B]","[B-]"))
for x in f4:
    input_smiles4.append(x.replace("\n","").replace("[N]","[N+]").replace("[O]","[O+]").replace("[B]","[B-]"))

# Close the input files
f1.close()
f2.close()
f3.close()
f4.close()

# Open the output files
f1 = open(path5,"w")
f2 = open(path6,"w")
f3 = open(path7,"w")
f4 = open(path8,"w")

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

