from FPSim2 import FPSim2Engine

path1 = "Commercial_MWlower330.csv"
path2 = "Commercial_MW330-500-1.csv"
path3 = "Commercial_MW330-500-2.csv"
path4 = "Commercial_MWhigher500.csv"

# List of input smiles strings
input_smiles = []

# Open the files
f1 = open(path1,"r")
f2 = open(path2,"r")
f3 = open(path3,"r")
f4 = open(path4,"r")

print("Start reading all the smiles...")

# Read all the smiles from f1
i=0
for x in f1:
    input_smiles.append(x.replace("\n",""))
    i += 1

# Read all the smiles from f2
for x in f2:
    input_smiles.append(x.replace("\n",""))
    i += 1

# Read all the smiles from f3
for x in f3:
    input_smiles.append(x.replace("\n",""))
    i += 1

# Read all the smiles from f4
for x in f4:
    input_smiles.append(x.replace("\n",""))
    i += 1

#close the files
f1.close()
f2.close()
f3.close()
f4.close()

fp_filename = 'fingerprints.h5'
fpe = FPSim2Engine(fp_filename)

query = 'CC(C)(CC(=O)NCc1nc(CC(=O)O)cs1)NC(=O)OCC2c3ccccc3c4ccccc24'

results = fpe.similarity(query, 0.5, n_workers=2)

for i in range(len(results)):
    print("Molecule: {}, Index: {}, Similarity: {}".format(input_smiles[results[i][0]-1], results[i][0], results[i][1]))

print("Len: {}".format(len(results)))