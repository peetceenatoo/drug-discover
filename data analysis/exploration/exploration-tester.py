# --------------------- Solve errors and warnings ---------------------- #

#     Currently this script produces a whole wall of errors and warnings
#     because molecule_generation module uses many deprecated functions.
#     Although, these all do not affect the correct functioning of
#     the script and will be fixed by Microsoft team before the deprecated
#     functions will have been deleted.

import os
import logging
from tensorflow.python.util import deprecation

# Disable logging output of tensorflow content [May be useless] 
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'

# Do not print deprecation warnings of tensorflow content
deprecation._PRINT_DEPRECATION_WARNINGS = False

# Set tensorflow logging level to only print fatal errors
logging.getLogger('tensorflow').setLevel(logging.FATAL) 

# ----------------------------------------------------------------------- #


# ------------------------------- Imports ------------------------------- #

from molecule_generation import load_model_from_directory
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import DataStructs
from FPSim2 import FPSim2Engine
import os
import random
import copy

# ------------------------------ Global ------------------------------- #

# Compute module of the magnitude of the size of the intervals
interval_size = 0
magnitude = 0
distributions = []

# ------------------------------ addNoise ------------------------------- #

def addNoise(embedding, area):

    global interval_size
    global magnitude
    global distributions

    # For each feature calculates its range for noise
    noise_intervals_lengths = []
    for i in range(len(distributions)):
        current_feature_distribution = distributions[i]
        current_interval = round(embedding[i] - embedding[i]%interval_size,magnitude)
        
        # Fetures in distribution tails are not touched
        intervals = sorted(current_feature_distribution.keys())
        left_tail = 0
        count = 0
        for j in range(len(intervals)):
            if count < area/2:
                count += current_feature_distribution[intervals[j]]
            else:
                left_tail = intervals[j]
                break
        right_tail = 0
        count = 0
        for j in range(len(intervals)):
            if count < area/2:
                count += current_feature_distribution[intervals[len(intervals)-1-j]]
            else:
                right_tail = intervals[j]
                break
        if current_interval < left_tail or current_interval > right_tail:
            noise_intervals_lengths.append(0.0)
            continue

        # Find the interval in which noise must be added
        if embedding[i]%interval_size < interval_size-embedding[i]%interval_size:
            even = copy.deepcopy(embedding)[i]%interval_size
            odd = interval_size - copy.deepcopy(embedding)[i]%interval_size
            is_near_lower = True
        else:
            even = interval_size - copy.deepcopy(embedding)[i]%interval_size
            odd = copy.deepcopy(embedding)[i]%interval_size
            is_near_lower = False

        lower_bound = copy.deepcopy(embedding)[i]
        upper_bound = copy.deepcopy(embedding)[i]
        current_area = 0
        j = 0
        while True:
            if j%2 == 0:
                p_lower = current_feature_distribution[round(current_interval-interval_size*j/2,magnitude)]
                p_upper = current_feature_distribution[round(current_interval+interval_size*j/2,magnitude)]
                if current_area + even/interval_size*(p_lower + p_upper) >= area:
                    break
                lower_bound -= even
                upper_bound += even
                current_area += even/interval_size*(p_lower + p_upper)
            else:
                if is_near_lower:
                    p_lower = current_feature_distribution[round(current_interval-interval_size*(1+(j-j%2)/2),magnitude)]
                    p_upper = current_feature_distribution[round(current_interval+interval_size*(j-j%2)/2,magnitude)]
                else:
                    p_lower = current_feature_distribution[round(current_interval-interval_size*(j-j%2)/2,magnitude)]
                    p_upper = current_feature_distribution[round(current_interval+interval_size*(1+(j-j%2)/2),magnitude)]
                if current_area + odd/interval_size*(p_lower + p_upper) >= area:
                    break
                lower_bound -= odd
                upper_bound += odd
                current_area += odd/interval_size*(p_lower + p_upper)
            j += 1
        
        r = (area-current_area)*interval_size/(p_lower+p_upper)
        noise_intervals_lengths.append(upper_bound-lower_bound+2*r)

    # Initialize embedding to be returned
    modified_embedding = copy.deepcopy(embedding)

    for j in range(len(embedding)):
        current_range = noise_intervals_lengths[j]
        modified_embedding[j] += random.random()*current_range-current_range/2

    # Return the embedding with noise addition
    return modified_embedding

# ------------------------------ readDistributions ------------------------------- #

def readDistributions():
    # Specify the paths for the database
    path = "..\\features_distribution_plots\\features_distributions.txt"

    # Open the file
    f = open(path,"r")

    # Read the size of the discretization range
    global interval_size
    global magnitude
    global distributions
    interval_size = float(f.readline().strip())
    while interval_size*(10**magnitude) < 1:
        magnitude += 1

    # Read the ranges for each feature and put them in the rows map
    for row in f:
            # Remove leading/trailing whitespace and newline characters
            row = row.strip()
            
            # Split the row by colon
            parts = row.split(":")
            bottom = float(parts[0])
            probs_part = parts[1]
            
            # Split the second part by semicolon
            probs = probs_part.split(";")
            probs = probs[:len(probs)-1]

            # Compute the map of probabilities for the current row
            temp_map = {}
            for i in range(len(probs)):
                temp_map[round(bottom+i*interval_size,magnitude)] = float(probs[i])
                
            # Put current map of probabilities in the rows map
            distributions.append(temp_map)
                
    # Close the file
    f.close()
    return

# ------------------------------ readDataset ------------------------------- #

def readDataset():

    dataset_paths = ["..\\..\\dataset\\Commercial_MW\\Commercial_MWlower330.csv", "..\\..\\dataset\\Commercial_MW\\Commercial_MW330-500-1.csv", "..\\..\\dataset\\Commercial_MW\\Commercial_MW330-500-2.csv", "..\\..\\dataset\\Commercial_MW\\Commercial_MWhigher500.csv"]

    # List of input smiles strings
    dataset_smiles = []

    for path in dataset_paths:

        # Open the files
        f = open(path,"r")

        # Read all the smiles from f
        for x in f:
            dataset_smiles.append(x.replace("\n",""))
            
        #close the files
        f.close()

    return dataset_smiles

# ------------------------------- Code ------------------------------- #

# Don't remember why but the if statement in the following row is needed,
# due to molecule_generation package usage
if __name__ == '__main__':

    # Read from distributions file
    readDistributions()

    # Read from dataset
    dataset_smiles = readDataset()

    # Specify the directory of the trained model
    model_dir = "..\\..\\model"

    # Specify the path for the file containing dataset fingerprints
    fp_filename = '..\\..\\fingerprints.h5'

    # Output files
    out1 = "noise_quality.txt"
    out2 = "project_quality.txt"

    # Input smiles
    input_smiles = ["COc1ccccc1N2CCN(CC2)C(=O)c3ccc(nc3)c4cccs4",
                    "CC(N1C(=S)S\C(=C/c2cccs2)\C1=O)C(=O)O",
                    "COC(=O)c1ccccc1NC(=O)NCCc2ccc(Cl)cc2",
                    "OC(CN1CCN(CC1)C(=O)c2ccc(nc2)n3ccnc3)c4ccccc4",
                    "Cc1cc(\C=C(/C#N)\c2ccccc2)c(C)n1c3ccccc3",
                    "O=C(C1CCCCN1S(=O)(=O)c2ccccc2)N3CCN(CC3)C(=O)c4ccccc4",
                    "Fc1ccc(cc1)S(=O)(=O)N2CCC(CC2)C(=O)NCCC(=O)NCc3cccnc3",
                    "COc1ccc(CNC(=O)N2CCCN(CC2)C(=O)NCc3ccc(OC)cc3)cc1",
                    "CCn1c(SCC(=O)c2cc(C)n(c2C)c3ccccc3)nnc1c4occc4",
                    "CC(=O)OC[C@H]1O[C@H](OC(=O)C)[C@@H](OC(=O)C)[C@@H](OC(=O)C)C1OC(=O)C",
                    "COc1ccc(NC(=O)C2CN(C)CC2c3ccccc3)cc1",
                    "CC(Nc1nc(Nc2ccc(F)cc2)nc(SCC(=O)Nc3ccc(F)cc3)n1)c4ccccc4",
                    "CCOc1ccc(CCNC(=O)COC(=O)c2ccncc2)cc1OCC",
                    "CCOC(=O)c1c(C)[nH]c(C)c1C(=O)CSc2ccccn2",
                    "O=C(Nc1ccccc1)Nc2cc(nn2c3ccccc3)C4CC4",
                    "CN(C)c1ccc(\C=N\\NC(=O)c2oc(Br)cc2)cc1",
                    "O[C@@H](c1ccnc2ccccc12)C(F)(F)F",
                    "CCOC(=O)N1CCN(CC1)C(=O)C2=CN(CC)c3ccc(cc3C2=O)S(=O)(=O)N4CCCCC4",
                    "CC(C)(C)c1ccc(cc1)C(=O)N2CCN(CC2)c3ncccc3C(F)(F)F",
                    "O[C@@H](c1ccc2ncccc2c1)C(F)(F)F",
                    "CCOc1ccc(\C=N\\NC(=O)c2cc(OC)c(OC)c(OC)c2)cc1OC",
                    "CN1CC(CC1=O)C(=O)NCc2ccc3CCCc3c2",
                    "O=C(CC1N(CCNC1=O)S(=O)(=O)c2ccccc2)Nc3ccccn3",
                    "Cc1c(sc2ccccc12)C(=O)NCCc3c(F)cccc3F",
                    "CC(C)(C)OC(=O)N1CCCC[C@H]1C(=O)NCC2CCCO2",
                    "CN(Cc1ccc(Cl)cc1)C(=O)c2oc3ccccc3c2",
                    "Fc1ccc2[nH]c(CCC(=O)NCC(N3CCOCC3)c4ccccn4)nc2c1",
                    "Clc1ccc(CC(=O)N2CCc3ccccc23)cc1",
                    "Cc1[nH]c(nc1C(=O)N(CCO)CCO)c2cccc(Cl)c2",
                    "Cc1cc(O)nc(SCC(=O)N2CCC(Cc3ccccc3)CC2)n1",
                    "O=C(N\\N=C\c1cccc2ccccc12)C(=O)Nc3ccccc3",
                    "CC(=O)N1CCC(C1)C(=O)Nc2n[nH]cc2Br",
                    "CC(N1CCN(Cc2ccc(NC(=O)C)cc2)CC1)c3ccccn3",
                    "O=C(Nc1ccccc1)Nc2cnn(CC3CCCCO3)c2",
                    "CCOC(=O)C1CCN(CC1)C(=O)c2csc(n2)c3ccc(C)cc3C",
                    "CC1=CN([C@@H]2O[C@H](CO)C(O)C2O)C(=O)NC1=O",
                    "CN(C)c1ccc(Nc2nc(cs2)c3ccc(Cl)cc3)cc1",
                    "CCn1c(SCc2cc(on2)c3ccccc3)nnc1c4ccoc4C",
                    "CC(C)CN1CCC(CC1)C(=O)NCCN2CCC(C)CC2"] 
    
    fout1 = open(out1, "w")
    fout2 = open(out2, "w")

    for input_smile in input_smiles:
        in_fps = AllChem.RDKFingerprint(Chem.MolFromSmiles(input_smile)) 

        print(f"Encoded: ", input_smile,"\n")
        fout1.write("I: {}\n".format(input_smile))
        fout2.write("I: {}\n".format(input_smile))

        # Using the model at model_dir path
        with load_model_from_directory(model_dir) as model:
            print()

            # Starting value for the area
            area = 0.0
            # Weights for finding/not finding duplicates
            present = 0.005
            not_present = -0.05

            print("Start encoding...")
            # Process latent vector for each input smiles string
            embedding = (model.encode([input_smile]))[0]

            # Adding noise to the embeddings
            list_of_decoded_smiles = []
            print("Start adding noise...")
            while area < 0.5:
                area = max(area,0.0)
                modified_molecule = (model.decode([addNoise(embedding, area)]))[0]
                if modified_molecule in list_of_decoded_smiles:
                    area += present
                else:
                    list_of_decoded_smiles.append(modified_molecule)
                    fout1.write("O: {}, {}\n".format(modified_molecule, DataStructs.TanimotoSimilarity(in_fps, AllChem.RDKFingerprint(Chem.MolFromSmiles(modified_molecule)))))
                    area += not_present

        # The first element is always the input molecule
        list_of_decoded_smiles = list_of_decoded_smiles[1:len(list_of_decoded_smiles)]

        # Similarity search
        print("Start similarity...")
        similarity_dataset = 0.5
        max_num_of_SMILES = 5
        output_smiles = []
        fpe = FPSim2Engine(fp_filename)
        for decoded_smile in list_of_decoded_smiles:
            results = fpe.similarity(decoded_smile, similarity_dataset, n_workers=2)
            results_smiles = []
            for res in results:
                if all(dataset_smiles[res[0]-1] != out_smi[0] for out_smi in output_smiles) and dataset_smiles[res[0]-1] != input_smile:
                    results_smiles.append(dataset_smiles[res[0]-1])
                    if len(results_smiles) >= max_num_of_SMILES:
                        break
            if len(results_smiles) == 0:
                continue
            res_fps = [AllChem.RDKFingerprint(Chem.MolFromSmiles(res_smi)) for res_smi in results_smiles]
            similarities = [DataStructs.TanimotoSimilarity(in_fps, out_fps) for out_fps in res_fps]
            output_smiles.append([results_smiles[similarities.index(max(similarities))], max(similarities)])
            
        output_smiles.sort(key=lambda x : x[1],reverse=True)
        # Print ouput
        for i in range(len(output_smiles)):
            print("Input molecule: {}\nOutput molecule: {}\nSimilarity: {}\n".format(input_smile, output_smiles[i][0], output_smiles[i][1]))
            fout2.write("O: {}, {}\n".format(output_smiles[i][0], output_smiles[i][1]))

        fout1.write("\n")
        fout2.write("\n")

    fout1.close()
    fout2.close()