import os
import logging
import tempfile
import subprocess
import csv


logger = logging.getLogger(__name__)



def make_blast_protein_database(input_database):
    """
    Creates the protein database for BLASTP - protein_blastp_search() in a temporary location.

    Args:
        input_database (str): Location of the input FASTA file protein sequences.

    Raises:
        FileNotFoundError: Input FASTA cannot be accessed.

    Returns:
        str: Location + prefix of BLAST protein database.
    """

    logger.debug('Creating BLAST protein database.')
    # Ensure the input file exists
    if not os.path.isfile(input_database):
        raise FileNotFoundError(f"The specified database FASTA file does not exist: {input_database}")


    # Set database prefix
    output_db_name = os.path.splitext(os.path.basename(input_database))[0]

    # Create a temporary directory for the database
    temp_dir = tempfile.mkdtemp()
    db_path = os.path.join(temp_dir, output_db_name)

    # Construct the makeblastdb command
    command = [
        'makeblastdb',
        '-in', input_database,
        '-dbtype', 'prot',
        '-out', db_path
    ]

    try:
        # Run the command and capture the output
        result = subprocess.run(command, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        
        # Log the output and error messages
        logger.info(result.stdout)
        
        logger.info(f"Protein BLAST database created successfully in temporary directory + prefix: {db_path}")

        # Return the database path, without the suffix
        return db_path
    except subprocess.CalledProcessError as e:
        logger.error(f"Error creating BLAST database: {e}")
        logger.error(e.stderr)
        raise
    except Exception as e:
        logger.error(f"An unexpected error occurred: {e}")
        raise
        



def protein_blastp_search(input_sequence, genome, output, input_database, threads):
    """Run BLASTP of the putative proteins against the predefined database.

    Args:
        input_sequence (str): Predicted proteins from the genome
        genome (str): Naming convention
        output (str): Output directory
        input_database (str): Location of the input FASTA file protein sequences.
        threads (_type_): Number of threads for the search to use
    """

    logger.debug('Entering protein_blastp_search function')

    # Get protein database
    protein_database = make_blast_protein_database(input_database)


    # Generate the BLASTP command
    blastp_command = [
        'blastp',
        '-db', protein_database,
        '-outfmt', '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore',
        '-query', input_sequence,
        '-num_threads', str(threads)
    ]

    output_csv_file = f'{output}/output_{genome}_protein_search.csv'

    ## Unclear what the subprocess.PIPE does

    try:
        blast_results = subprocess.run(blastp_command, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

        # Split the output into lines and write to CSV file
        with open(output_csv_file, 'w', newline='') as csvfile:
            csvwriter = csv.writer(csvfile, delimiter=',')
            # Write header
            csvwriter.writerow(["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore"])
            # Write rows
            for line in blast_results.stdout.strip().split('\n'):
                csvwriter.writerow(line.split('\t'))
        
        ## Include in the -q flag

        # print(f"P-blast results saved to {output_csv_file}")
        logger.info(f'pblast result saved in file: {output_csv_file}')
    except subprocess.CalledProcessError as e:
        print("P-blast failed with the following error message:\n", e.stderr)
        logger.error(f'Error in rotein_blastp_search (subprocess): {e}')
    except Exception as ex:
        print("An error occurred:", ex)
        print("BLASTP Output:\n", blast_results.stdout)
        logger.error(f'Error in rotein_blastp_search: {ex}')
    except KeyboardInterrupt:
        logger.warning("Data processing interrupted by user")
    logger.debug('Exiting protein_blastp_search function')