import pandas as pd
from copy import deepcopy
import os
from Bio import SeqIO
from Bio.SeqUtils import molecular_weight

## First dereplicating function


def dereplicate_highest_score(df):

    df = pd.read_csv(df)

    # """
    # Dereplicate the DataFrame by keeping the highest score for each unique 'Name1'.

    # Parameters:
    # df (DataFrame): Input DataFrame with columns 'Name1' and 'Score'.

    # Returns:
    # DataFrame: Dereplicated DataFrame with the highest score for each 'Name1'.
    # """
    dataframe = deepcopy(df)
    sorted_df = dataframe.sort_values(by="Score", ascending=False)
    dereplicated_df = dataframe.drop_duplicates(
        subset="Name1", keep="first"
    ).reset_index(drop=True)
    return dereplicated_df


## Second, calculate statistics


def calculate_mass_length(fasta_loc, df_entry, pblast_file_path):

    # Read the FASTA file
    sequences = list(SeqIO.parse(fasta_loc, "fasta"))
    sequence_dict = {record.id: record for record in sequences}

    # Initialize lists to store the statistics
    lengths = []
    masses = []

    # Loop through the DataFrame and calculate statistics
    for name in df_entry["Name1"]:
        # Match the entry in the FASTA file
        if name in sequence_dict:
            seq_record = sequence_dict[name]
            # Remove asterisks from the protein sequence
            cleaned_sequence = str(seq_record.seq).replace("*", "")
            sequence_length = len(cleaned_sequence)
            sequence_mass = molecular_weight(cleaned_sequence, seq_type="protein")

            lengths.append(sequence_length)
            masses.append(sequence_mass)
        else:
            raise Exception(
                "Entry in results does not match any entry in protein fasta file."
            )

    # Add the statistics to the DataFrame
    df_entry["Length"] = lengths
    df_entry["Normalized_score"] = df_entry["Score"] / df_entry["Length"]
    df_entry["Mass"] = masses

    pblast_df = pd.read_csv(pblast_file_path)

    # Rename columns in the pBLAST DataFrame for merging
    pblast_df.rename(columns={"qseqid": "Name1", "sseqid": "Name2"}, inplace=True)

    # Merge the DataFrames on 'Name1' and 'Name2'
    df_entry = pd.merge(
        df_entry,
        pblast_df[["Name1", "Name2", "evalue"]],
        on=["Name1", "Name2"],
        how="left",
    )

    df_entry = df_entry.sort_values(by="Normalized_score", ascending=False)

    # df_entry.to_csv(f'{output}/final_results_{genome}.csv')
    return df_entry
