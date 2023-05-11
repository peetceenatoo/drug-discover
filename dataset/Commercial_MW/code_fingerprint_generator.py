# ------------------------------- Imports ------------------------------- #

from FPSim2.io import create_db_file

# ------------------------------- Code ------------------------------- #

# Don't remember why but the if statement in the following row is needed,
# due to molecule_generation package usage
if __name__ == '__main__':

    # Specify the paths for the database
    path1 = "Commercial_MWlower330.csv"
    path2 = "Commercial_MW330-500-1.csv"
    path3 = "Commercial_MW330-500-2.csv"
    path4 = "Commercial_MWhigher500.csv"
    out1 = "fingerprints.h5"

    # List of input smiles strings
    input_smiles = []

    # Open the files
    f1 = open(path1,"r")
    f2 = open(path2,"r")
    f3 = open(path3,"r")
    f4 = open(path4,"r")

    print("Start reading all the smiles...")

    # Read all the smiles from f1
    i=1
    for x in f1:
        input_smiles.append([x.replace("\n",""),0+i])
        i += 1

    # Read all the smiles from f2
    for x in f2:
        input_smiles.append([x.replace("\n",""),0+i])
        i += 1

    # Read all the smiles from f3
    for x in f3:
        input_smiles.append([x.replace("\n",""),0+i])
        i += 1

    # Read all the smiles from f4
    for x in f4:
        input_smiles.append([x.replace("\n",""),0+i])
        i += 1

    #close the files
    f1.close()
    f2.close()
    f3.close()
    f4.close()

    print("Start creating the database... Good luck!")

    create_db_file(input_smiles,out1,'Morgan',{'radius': 3,'nBits': 2048})
