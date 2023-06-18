# ------------------------------- Imports ------------------------------- #



# ------------------------------- Code ------------------------------- #

# Specify the paths for the database
path = "..\\features_distribution_plots\\features_distributions.txt"

# Open the file
f = open(path,"r")

# Read the size of the discretization range
range = float(f.readline().strip())

# Read the ranges for each feature and put them in the rows map
features_list = []
for row in f:
        # Remove leading/trailing whitespace and newline characters
        row = row.strip()
        
        # Split the row by colon
        parts = row.split(":")
        bottom = parts[0]
        probs_part = parts[1]
        
        # Split the second part by semicolon
        probs = probs_part.split(";")

        # Compute the map of probabilities for the current row
        temp_map = {}
        for i in range(len(probs)):
            temp_map[bottom+i*range] = float(probs[i])
            
        # Put current map of probabilities in the rows map
        features_list.append(temp_map)
            
# Close the file
f.close()
