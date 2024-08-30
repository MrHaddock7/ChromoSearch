import subprocess

# Function to run Prodigal and return the protein sequences 
# Input: genome (FASTA format)
# Output: Protein candidates (also fasta format)

def run_prodigal(input_file, output_prot_file, gene):
    # Command to run Prodigal and output protein sequences
    command = ['prodigal', '-i', input_file, '-a', f'{output_prot_file}/output_{gene}_DNAtoProtein.fasta', '-q']

    try:
        subprocess.run(command, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        # print(f"Prodigal finished successfully. Protein sequences saved to {output_prot_file}/output_{gene}_DNAtoProtein.fasta")
    except subprocess.CalledProcessError as e:
        print("Error running Prodigal:", e)
