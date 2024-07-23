import subprocess

# Function to run Prodigal and return the protein sequences 
# Input: genome (FASTA format)
# Output: Protein candidates (also fasta format)

def run_prodigal(input_file, output_prot_file, gene):
    # Command to run Prodigal and output protein sequences
    command = ['prodigal', '-i', input_file, '-a', f'{output_prot_file}/output_{gene}_DNAtoProtein.fasta']

    try:
        subprocess.run(command, check=True)
        print(f"Prodigal finished successfully. Protein sequences saved to {output_prot_file}/output_{gene}_DNAtoProtein.fasta")
    except subprocess.CalledProcessError as e:
        print("Error running Prodigal:", e)

# # Example usage:
# input_file = r'/Users/klonk/Desktop/Chromophore/Chromoproteins_2024/secret/genome_data/e_coli_k12_genome.fasta'
# output_prot = 'output_proteins_k12.faa'
# output_cdf = "output_cdf.gff"
