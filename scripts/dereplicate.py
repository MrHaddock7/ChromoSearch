import pandas as pd
from copy import deepcopy
import glob
import os
import sys


def dereplicate_highest_score(df):
    # Sort dataframe by Name1 and Scores in descending order
    dataframe = deepcopy(df)
    sorted_df = dataframe.sort_values(by='Score', ascending=False)
    
    # Drop duplicates keeping the first occurrence (which has the highest score due to sorting)
    deduplicated_df = sorted_df.drop_duplicates(subset='Name1', keep='first').reset_index(drop=True)
    
    return deduplicated_df


def load_csv_files(directory):
    # Construct the search pattern
    pattern = os.path.join(directory, '*sorter2.csv')
    print(pattern)
    # Find all files matching the pattern
    csv_files = glob.glob(pattern)
    
    # Load each CSV file into a DataFrame and store them in a list of tuples (filename, DataFrame)
    dataframes = {}
    for file in csv_files:
        # Extract the base filename without the "sorter2.csv" part
        base_filename = os.path.basename(file).replace('sorter2.csv', '')
        df = pd.read_csv(file)
        dataframes[base_filename] = df
    
    return dataframes


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python dereplicate.py <directory> <output_dir>")
        sys.exit(1)

    directory = sys.argv[1]
    output_dir = sys.argv[2]

    if not os.path.isdir(directory):
        print(f"The directory {directory} does not exist.")
        sys.exit(1)

    dataframes = load_csv_files(directory)
    dereplicated_dfs = {}

    # dereplicate all dataframes, and save in output dir
    for key, value in dataframes.items():
        print(key)
        dereplicated_df = dereplicate_highest_score(value)
        output_loc = os.path.join(output_dir, f"{key}_dereplicated_sorted.csv")
        dereplicated_df.to_csv(output_loc, index=False)

    sys.exit(0)
