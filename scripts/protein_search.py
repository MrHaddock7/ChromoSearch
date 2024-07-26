import subprocess
import csv
import pandas as pd
import logging

logger = logging.getLogger(__name__)

def protein_blastp_search(input_sequence, genome, output, input_database):
    logger.debug('Entering protein_blastp_search function')
    # Generate the BLASTP command
    blastp_command = [
        'blastp',
        '-db', input_database,
        '-outfmt', '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore',
        '-query', input_sequence
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