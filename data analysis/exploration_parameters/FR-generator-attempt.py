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
            even = embedding[i]%interval_size
            odd = interval_size-embedding[i]%interval_size
            is_near_lower = True
        else:
            even = interval_size-embedding[i]%interval_size
            odd = embedding[i]%interval_size
            is_near_lower = False

        lower_bound = embedding[i]
        upper_bound = embedding[i]
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
    path1 = "..\\..\\dataset\\Commercial_MW\\Commercial_MWlower330.csv"
    path2 = "..\\..\\dataset\\Commercial_MW\\Commercial_MW330-500-1.csv"
    path3 = "..\\..\\dataset\\Commercial_MW\\Commercial_MW330-500-2.csv"
    path4 = "..\\..\\dataset\\Commercial_MW\\Commercial_MWhigher500.csv"

    # List of input smiles strings
    dataset_smiles = []

    # Open the files
    f1 = open(path1,"r")
    f2 = open(path2,"r")
    f3 = open(path3,"r")
    f4 = open(path4,"r")

    print("Start reading all the smiles...")

    # Read all the smiles from f1
    i=0
    for x in f1:
        dataset_smiles.append(x.replace("\n",""))
        i += 1

    # Read all the smiles from f2
    for x in f2:
        dataset_smiles.append(x.replace("\n",""))
        i += 1

    # Read all the smiles from f3
    for x in f3:
        dataset_smiles.append(x.replace("\n",""))
        i += 1

    # Read all the smiles from f4
    for x in f4:
        dataset_smiles.append(x.replace("\n",""))
        i += 1

    #close the files
    f1.close()
    f2.close()
    f3.close()
    f4.close()
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
    fp_filename = '..\\..\\dataset\\Commercial_MW\\fingerprints.h5'

    # Input smile
    input_smile = "COc1ccccc1N2CCN(CC2)C(=O)c3ccc(nc3)c4cccs4"
    in_fps = AllChem.RDKFingerprint(Chem.MolFromSmiles(input_smile)) 

    print(f"Encoded: ",input_smile,"\n")

    # Using the model at model_dir path
    with load_model_from_directory(model_dir) as model:
        print()

        # Number of molecules to be generated from each input smiles string
        num_of_molecules_to_generate = 20
        # Starting value for the area
        area = 0.01
        # Weights for finding/not finding duplicates
        present = 0.01
        not_present = -0.05

        print("Start encoding...")
        # Process latent vector for each input smiles string
        embedding = (model.encode([input_smile]))[0]

        # Calculate num_of_molecules_to_generate diverse molecules for each input smiles string,
        # adding noise to the embeddings
        list_of_decoded_smiles = []
        print("Start adding noise...")
        while len(list_of_decoded_smiles) < num_of_molecules_to_generate and area < 0.5:
            modified_molecule = (model.decode([addNoise(embedding,area)]))[0]
            if modified_molecule in list_of_decoded_smiles:
                area += present
            else:
                list_of_decoded_smiles.append(modified_molecule)
                area += not_present

    # The first element is always the input molecule
    list_of_decoded_smiles = list_of_decoded_smiles[1:len(list_of_decoded_smiles)]

    # Similarity search
    similarity_dataset = 0.5
    num_of_similar_smiles_per_decoded_smiles = 5
    output_smiles = []
    fpe = FPSim2Engine(fp_filename)
    for decoded_smile in list_of_decoded_smiles:
        results = fpe.similarity(decoded_smile, similarity_dataset, n_workers=2)
        results_smiles = []
        for res in results:
            if all(dataset_smiles[res[0]-1] != out_smi[0] for out_smi in output_smiles) and dataset_smiles[res[0]-1] != input_smile:
                results_smiles.append(dataset_smiles[res[0]-1])
                if len(results_smiles) >= num_of_similar_smiles_per_decoded_smiles:
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