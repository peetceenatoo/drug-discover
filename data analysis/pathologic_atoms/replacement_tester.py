# ------------------------------- Imports ------------------------------- #

from rdkit import Chem

# ------------------------------- Code ------------------------------- #

# Specify the paths for the database
path1 = "..\\..\\dataset\\Commercial_MW\Commercial_MWlower330(clean).csv"
path2 = "..\\..\\dataset\\Commercial_MW\Commercial_MW330-500(clean)1.csv"
path3 = "..\\..\\dataset\\Commercial_MW\Commercial_MW330-500(clean)2.csv"
path4 = "..\\..\\dataset\\Commercial_MW\Commercial_MWhigher500(clean).csv"

# Open the files
f1 = open(path1, "r")
f2 = open(path2, "r")
f3 = open(path3, "r")
f4 = open(path4, "r")

# List of input smiles strings
input_smiles = []

# Read all the smiles from the input files
cont = 0
for x in f1:
    cont = cont+1
    if "[O]" in x:
        input_smiles.append(x.replace("\n", ""))
for x in f2:
    if "[O]" in x:
        cont = cont+1
        input_smiles.append(x.replace("\n", ""))
for x in f3:
    if "[O]" in x:
        cont = cont+1
        input_smiles.append(x.replace("\n", ""))
for x in f4:
    if "[O]" in x:
        cont = cont+1
        input_smiles.append(x.replace("\n", ""))

# Close the files
f1.close()
f2.close()
f3.close()
f4.close()

# Print total number of molecules in the database
print("Overall there are {} molecules in the dataset and ".format(cont))

# Only take the ones which contain a pathologic pattern: the last "[O]"s seem to only be pathologic in the following cases
chosen_smiles = []
for i in range(len(input_smiles)):
    if "N([O])" in input_smiles[i] or "N[O]" in input_smiles[i] or "N1[O]" in input_smiles[i] or "N2[O]" in input_smiles[i] or "N3[O]" in input_smiles[i] or "N4[O]" in input_smiles[i] or "N5[O]" in input_smiles[i]:
        chosen_smiles.append(input_smiles[i])
# Print the number of chosen_smiles
print("{} contain [O] within a pathologic pattern.\n".format(len(chosen_smiles)))

# Replace pattern
changed_smiles = []
print("I'm computing the dynamic replacement...")
for i in range(len(chosen_smiles)):
    changed_smiles.append(chosen_smiles[i].replace("[O]","[O+]"))     

# Print chosen smiles and trigger an eventual error
print("################ Print Chosen #####################\n")
for i in range(len(chosen_smiles)):
    # CHOSEN
    print("{}% ".format(i/len(changed_smiles)*100))
    print(chosen_smiles[i])
    Chem.MolFromSmiles(chosen_smiles[i])
    print()

# Print changed smiles and trigger an eventual error (if it worked, no error will show up)
print("################ Print Changed #####################\n")
for i in range(len(changed_smiles)):
    # CHANGED
    print("{}% ".format(i/len(changed_smiles)*100))
    print(changed_smiles[i])
    Chem.MolFromSmiles(changed_smiles[i])
    print()

# Done
print("If no error popped up this way, it means you should replace this pathologic pattern!")


