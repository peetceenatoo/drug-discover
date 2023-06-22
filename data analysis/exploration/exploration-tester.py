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

# Global variables to compute the magnitude's module of the intervals' size
interval_size = 0
magnitude = 0
distributions = []

# ------------------------------ addNoise ------------------------------- #

def addNoise(embedding, area, tail):
    global interval_size
    global magnitude
    global distributions

    #     # For each feature calculate the range of the interval where to add noise
    noise_intervals_lengths = []
    for i in range(len(distributions)):
        current_feature_distribution = distributions[i]
        current_interval = round(embedding[i] - embedding[i]%interval_size, magnitude)
        
        # Get intervals in current feature's discrete domain
        intervals = sorted(current_feature_distribution.keys())
        
        left_tail = 0
        count = 0
        # Find interval in current feature's domain such that the probability from left accumulates to more than area*tail
        for j in range(len(intervals)):
            if count < area*tail:
                count += current_feature_distribution[intervals[j]]
            else:
                left_tail = intervals[j]
                break
        
        right_tail = 0
        count = 0
        # Find interval in current feature's domain such that the probability from right accumulates to more than area*tail
        for j in range(len(intervals)):
            if count < area*tail:
                count += current_feature_distribution[intervals[len(intervals)-1-j]]
            else:
                right_tail = intervals[j]
                break
        if current_interval < left_tail or current_interval > right_tail:
            noise_intervals_lengths.append(0.0)
            continue

        # Compute range of values, centered in embedding[i], such that overall probability (approximatively) amounts to area

        # Init even, odd (deltas to be added at even and odd iterations)
        # and is_near_lower (indicates whether embedding[i] is closer to the left or right bound of the domain's interval which contains it)
        if embedding[i]%interval_size < interval_size-embedding[i]%interval_size:
            even = copy.deepcopy(embedding)[i]%interval_size
            odd = interval_size - copy.deepcopy(embedding)[i]%interval_size
            is_near_lower = True
        else:
            even = interval_size - copy.deepcopy(embedding)[i]%interval_size
            odd = copy.deepcopy(embedding)[i]%interval_size
            is_near_lower = False

        # Starting from embedding[i], compute the range
        lower_bound = copy.deepcopy(embedding)[i]
        upper_bound = copy.deepcopy(embedding)[i]
        current_area = 0
        j = 0
        while True:
            # If current iteration is even, add even delta to left and right, and increase 
            if j%2 == 0:
                # Get probability of the interval which was reached to the left (j=0 => starting interval)
                p_lower = current_feature_distribution[round(current_interval-interval_size*j/2, magnitude)]
                # Get probability of the interval which was reached to the right (j=0 => starting interval)
                p_upper = current_feature_distribution[round(current_interval+interval_size*j/2, magnitude)]

                # Set reached area
                current_area += even/interval_size*(p_lower + p_upper)
                # If done, break
                if current_area >= area:
                    break

                # Set new bounds
                lower_bound -= even
                upper_bound += even
                  
            # If current iteration is even, add even delta to left and right 
            else:
                # If embedding[i] was closer to the left bound of its domain's interval, intervals to the left are visited one iteration ahead 
                if is_near_lower:
                    p_lower = current_feature_distribution[round(current_interval-interval_size*(1+(j-1)/2), magnitude)]
                    p_upper = current_feature_distribution[round(current_interval+interval_size*(j-1)/2, magnitude)]
                # Othwerwise, intervals to the right are visited one iteration ahead
                else:
                    p_lower = current_feature_distribution[round(current_interval-interval_size*(j-1)/2, magnitude)]
                    p_upper = current_feature_distribution[round(current_interval+interval_size*(1+(j-1)/2), magnitude)]
                    
                # Set reached area
                current_area += odd/interval_size*(p_lower + p_upper)
                # If done, break
                if current_area >= area:
                    break

                # Set new bounds
                lower_bound -= odd
                upper_bound += odd
                current_area += odd/interval_size*(p_lower + p_upper)
            j += 1
        
        # Compute delta to be added to the current range and store the final range in noise_interval_lengths
        r = (area-current_area)*interval_size/(p_lower+p_upper)
        noise_intervals_lengths.append(upper_bound-lower_bound+2*r)

    # Initialize embedding to be returned
    modified_embedding = copy.deepcopy(embedding)

    # For each feature, add random shift
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
    f = open(path, "r")

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
                temp_map[round(bottom+i*interval_size, magnitude)] = float(probs[i])
                
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
        f = open(path, "r")

        # Read all the smiles from f
        for x in f:
            dataset_smiles.append(x.replace("\n", ""))
            
        #close the files
        f.close()

    return dataset_smiles

# ------------------------------- Code ------------------------------- #

# Don't remember why but the if statement in the following row is needed,
# due to molecule_generation package usage
if __name__ == '__main__':

    print("Hi! The execution of this script may take some time. Please, get comfortable and wait!")

    # Read from distributions file
    readDistributions()

    # Read from dataset
    dataset_smiles = readDataset()

    # Specify the directory of the trained model
    model_dir = "..\\..\\model"

    # Specify the path for the file containing dataset fingerprints
    fp_filename = '..\\..\\fingerprints.h5'

    # Output files
    out = "project_quality.txt"

    # Input smiles
    input_smiles = ["COc1ccccc1N2CCN(CC2)C(=O)c3ccc(nc3)c4cccs4"] 
    
    fout = open(out, "w")

    # For each smile in input_smiles
    for input_smile in input_smiles:
        # Compute its fingerprint
        in_fps = AllChem.RDKFingerprint(Chem.MolFromSmiles(input_smile)) 

        print(f"Encoded: ", input_smile,"\n")
        fout.write("I: {}\n".format(input_smile))
        
        # Using the model at model_dir path
        with load_model_from_directory(model_dir) as model:
            print()

            # Starting value for the area
            area = 0.0
            # Bonus/Penalty for finding duplicates/no duplicates
            present = 0.005
            not_present = -0.05
            # Max percentage to stop the algorithm when reached
            p_max = 0.5

            print("Start encoding...")
            # Process latent vector for each input smiles string
            embedding = (model.encode([input_smile]))[0]

            # Generate molecules for each input smiles string, adding noise to the embeddings
            list_of_decoded_smiles = []
            print("Start adding noise...")
            # Until area reaches
            while area < p_max:
                # Prevent area from going below zero
                area = max(area, 0.0)

                # Compute two molecules, with different values for tail
                modified_molecule1 = (model.decode([addNoise(embedding, area, 1/2)]))[0]
                modified_molecule2 = (model.decode([addNoise(embedding, area, 2/3)]))[0]

                # If both had been found already, add present bonus
                if (modified_molecule1 in list_of_decoded_smiles) and (modified_molecule2 in list_of_decoded_smiles):
                    area += present
                # Otherwise, append the new molecule and add not_present penalty
                else:
                    if (modified_molecule1 not in list_of_decoded_smiles):
                        list_of_decoded_smiles.append(modified_molecule1)
                    if (modified_molecule2 not in list_of_decoded_smiles):
                        list_of_decoded_smiles.append(modified_molecule2)
                    area += not_present

        # The first element is always the input molecule
        list_of_decoded_smiles = list_of_decoded_smiles[1:len(list_of_decoded_smiles)]

        # Similarity search
        print("Start similarity...")
        # Threshold
        similarity_dataset = 0.5
        # Max num of similar SMILES to be considered for each generated molecule
        max_num_of_SMILES = 5

        # Compute similarity search for each generated molecule
        output_smiles = []
        fpe = FPSim2Engine(fp_filename)
        for decoded_smile in list_of_decoded_smiles:
            # Compute similarity using similarity_dataset as a threshold
            results = fpe.similarity(decoded_smile, similarity_dataset, n_workers=2)
            results_smiles = []
            for res in results:
                # Prevent duplicates from being appended
                if all(dataset_smiles[res[0]-1] != out_smi[0] for out_smi in output_smiles) and dataset_smiles[res[0]-1] != input_smile:
                    results_smiles.append(dataset_smiles[res[0]-1])
                    # If the max_num is reached, stop the algorithm
                    if len(results_smiles) >= max_num_of_SMILES:
                        break
            # If no one was found, go to the next generated molecule
            if len(results_smiles) == 0:
                continue
            # Compute similarity between each result and the input molecule
            res_fps = [AllChem.RDKFingerprint(Chem.MolFromSmiles(res_smi)) for res_smi in results_smiles]
            similarities = [DataStructs.TanimotoSimilarity(in_fps, out_fps) for out_fps in res_fps]
            # And append the most similar to the output_smiles
            output_smiles.append([results_smiles[similarities.index(max(similarities))], max(similarities)])
            
        # Sort by similarity
        output_smiles.sort(key=lambda x : x[1],reverse=True)
        
        # Print ouput
        for i in range(len(output_smiles)):
            print("Input molecule: {}\nOutput molecule: {}\nSimilarity: {}\n".format(input_smile, output_smiles[i][0], output_smiles[i][1]))
            fout.write("O: {}, {}\n".format(output_smiles[i][0], output_smiles[i][1]))
        
        fout.write("\n")

    fout.close()