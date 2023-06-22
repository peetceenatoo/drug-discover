
# Specify the paths for the database
path = "project_quality.txt"

print("# ----------- Let's compute data about the noise ----------- #")

# Open the first file
f = open(path,"r")

# To compute average number of outputs
input_count = 0
output_count = 0

# To compute counters
input_quality = {}      # temp   
counters = {}           # counters

# Read the ranges for each feature and put them in the rows map
for row in f:
    # Remove leading/trailing whitespace and newline characters
    row = row.strip()
    
    # New input string
    if row and row[0] == "I":
        current_string = row.split(" ")[1]
        print("Current string: {}".format(current_string))
        input_count += 1

    # New output string (save its quality)
    elif row and row[0] == "O":
        output_count += 1
        input_quality[row.split(" ")[1]] = float(row.split(" ")[2])

    # While, if done with the current input string, calculate quantities for current input string
    else:
        # Set counters to 0
        counters[current_string] = {}
        for x in [0.9, 0.7, 0.5]:
            counters[current_string][x] = 0
        # Count
        for output in input_quality:
            if input_quality[output] >= 0.9:
                counters[current_string][0.9] += 1
            elif input_quality[output] >= 0.7:
                counters[current_string][0.7] += 1
            elif input_quality[output] >= 0.5:
                counters[current_string][0.5] += 1
        # Print
        print("Above 0.9 there are {} mols".format(counters[current_string][0.9]))
        print("Above 0.7 there are {} mols".format(counters[current_string][0.7]))
        print("Above 0.5 there are {} mols".format(counters[current_string][0.5]))
        print()
        # Reset input qualities for current string
        input_quality = {}

# Print final overview
print("Average number of outputs for each input mol: {}".format(output_count/input_count))
# Set sum to zero
sum = {}
for x in [0.9, 0.7, 0.5]:
    sum[x] = 0
# Count
for input in counters:
    sum[0.5] += counters[input][0.5]
    sum[0.7] += counters[input][0.7]
    sum[0.9] += counters[input][0.9]
# Print average
print("Average above 0.9 for each input mol: {}".format(sum[0.9]/input_count))
print("Average above 0.7 for each input mol: {}".format(sum[0.7]/input_count))
print("Average above 0.5 for each input mol: {}".format(sum[0.5]/input_count))

# Close the file
f.close()