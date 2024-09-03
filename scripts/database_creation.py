from Bio.Blast.Applications import NcbimakeblastdbCommandline
import os
import shutil

def create_blast_db(fasta_file, db_name, db_output):
    print("enter")
    """
    Convert a FASTA file into a BLAST database.

    Parameters:
    fasta_file (str): Path to the input FASTA file.
    db_name (str): name of BLAST database, such as 'positive_control'.
    db_output (str): directory to the output database, such as '/Users/klonk/Desktop/positive_n_negative_control/actual_database'

    Returns:
    str: Path to the created BLAST database folder.
    """

    # Check if the FASTA file exists
    if not os.path.exists(fasta_file):
        raise FileNotFoundError(f"FASTA file '{fasta_file}' not found.")

    # Construct the command for creating the BLAST database
    makeblastdb_cmd = NcbimakeblastdbCommandline(
        cmd='makeblastdb', 
        input_file=fasta_file, 
        dbtype="prot", 
        out=f"{db_output}/{db_name}"
    )

    # Execute the command
    stdout, stderr = makeblastdb_cmd()

    # Check for errors
    if stderr:
        raise Exception(f"Error occurred while creating BLAST database: {stderr}")
    
    shutil.copy(fasta_file, db_output)

    # Return the path to the database folder
    return f"{db_name}.pin"

# Example usage
fasta_file = "/Users/klonk/Desktop/positive_n_negative_control/positive_control_database.fasta"
db_name = "positive_control"
db_output = f'/Users/klonk/Desktop/positive_n_negative_control/actual_database'

try:
    db_path = create_blast_db(fasta_file, db_name, db_output)
    print(f"BLAST database created successfully: {db_path}")
except Exception as e:
    print(str(e))
