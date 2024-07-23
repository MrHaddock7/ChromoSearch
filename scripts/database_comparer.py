import csv

# Function that extracts the protein names from a fasta file, usually from a protein query on uniprot

def extract_protein_names(file_path):
    protein_names = []
    with open(file_path, 'r') as file:
        for line in file:
            if line.startswith('>'):
                # Extract everything from '>' until the first space
                protein_name = line[1:line.find(' ')]  # Use find to locate the first space
                protein_names.append(protein_name)
    return protein_names

# Example usage:
file_path = '/Users/klonk/Desktop/uniprotkb_carotenoids_2024_07_12.fasta'  # replace with your CSV file path
protein_names_from_database = extract_protein_names(file_path)

# Define the file name
file_name = '/Users/klonk/Desktop/Chromophore/Chromoproteins_2024/Output/Run_1/output_PT298_sorter.csv'
protein_names_from_sorter2 = []

# Read the CSV file (sorter2.csv type file)
with open(file_name, 'r') as file:
    reader = csv.reader(file)
    for row in reader:
        if len(row) > 1:  # Ensure there is a second column
            protein_names_from_sorter2.append(row[1])

# print(protein_names_from_sorter2[+])

# Function that compares the query to the outputs in the sorter2.csv type file

def list1inlist2(list1_from_waterman, list2_from_query):
    matches = []
    for item in list1_from_waterman:
        if item in list2_from_query:
            matches.append(item)
    return matches

output_list = list1inlist2(protein_names_from_sorter2, protein_names_from_database)
print(len(output_list))
if len(output_list) > 0:
    print(output_list)