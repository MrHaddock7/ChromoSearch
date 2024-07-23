import pandas as pd
import logging

logger = logging.getLogger(__name__)

def name_and_sequence_pair(input_genome_fasta, alignment_references, input_database_fasta):
    logger.debug('Entering pname_and_sequence_pair function')
    def parse_fasta_file(file_path):
        logger.debug('Entering parse_fasta_file function')
        try:
            sequences = []
            with open(file_path, 'r') as f:
                lines = f.readlines()
                current_name = None
                current_sequence = []
                for line in lines:
                    line = line.strip()
                    if line.startswith('>'):
                        if current_name is not None:
                            sequences.append((current_name, ''.join(current_sequence)))
                        current_name = line[1:].split()[0]  # Extracting the protein name
                        current_sequence = []
                    else:
                        current_sequence.append(line)
                # Append the last sequence
                if current_name is not None and current_sequence:
                    sequences.append((current_name, ''.join(current_sequence)))
        except Exception as e:
            logger.error(f'Error in parse_fasta_file function: {e}')
        except KeyboardInterrupt:
            logger.warning("Data processing interrupted by user")
        logger.debug('Exiting parse_fasta_file function')
        return sequences
    try:
        fasta_file = input_genome_fasta
        fasta_file2 = input_database_fasta

        sequences = parse_fasta_file(fasta_file)
        sequences2 = parse_fasta_file(fasta_file2)

        df = pd.DataFrame(sequences, columns=['Protein Name', 'Sequence'])
        df2 = pd.DataFrame(sequences2, columns=['Protein Name', 'Sequence'])
        df.set_index('Protein Name', inplace=True)
        df2.set_index('Protein Name', inplace=True)

        sorted_output_df = pd.read_csv(alignment_references)

        final_list = []
        for i in range(len(sorted_output_df)):
            final_list.append([(sorted_output_df.iloc[i, 0], df.loc[sorted_output_df.iloc[i, 0], 'Sequence']), (sorted_output_df.iloc[i, 1], df2.loc[sorted_output_df.iloc[i, 1], 'Sequence'])])
    except Exception as ex:
        logging.error(f'Error in name_and_sequence_pair function: {ex}')
    logger.debug('Exiting pname_and_sequence_pair function')
    return final_list