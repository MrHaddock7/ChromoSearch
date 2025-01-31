import pandas as pd
import logging

logger = logging.getLogger(__name__)

def name_and_sequence_pair(input_genome_fasta, alignment_references, input_database_fasta, blastpsw=True):
    
    logger.debug('Entering pname_and_sequence_pair function')
    def parse_fasta_file(file_path): # Input, path to FASTA file you want to parse
        logger.debug('Entering parse_fasta_file function')
        try:
            sequences = []
            with open(file_path, 'r') as f:
                lines = f.readlines()
                current_name = None 
                current_sequence = []
                for line in lines:
                    line = line.strip()
                    if line.startswith('>'): # Adding the name of the protein
                        if current_name is not None:
                            sequences.append((current_name, ''.join(current_sequence)))
                        current_name = line[1:].split()[0]  # Extracting the protein name
                        current_sequence = []
                    else:
                        current_sequence.append(line) # Adding the protein sequence
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

        sequences = parse_fasta_file(fasta_file) # Parsing the genome 
        sequences2 = parse_fasta_file(fasta_file2) # Parsing the database

        df = pd.DataFrame(sequences, columns=['Protein Name', 'Sequence'])
        df2 = pd.DataFrame(sequences2, columns=['Protein Name', 'Sequence'])
        df.set_index('Protein Name', inplace=True)
        df2.set_index('Protein Name', inplace=True)

        final_list = [] 

        sorted_output_df = pd.read_csv(alignment_references)

        ## Output pairs based on the highest scoring results from BlastP

        if blastpsw:
            logger.info('Paring from BlastP results')

            try:
                for i in range(len(sorted_output_df)):
                    index_1 = sorted_output_df.iloc[i, 0]
                    index_2 = sorted_output_df.iloc[i, 1]
 
                    seq_1 = df.loc[index_1, 'Sequence']
                    seq_2 = df2.loc[index_2, 'Sequence']
                    
                    final_list.append([(index_1, str(seq_1)), (index_2, str(seq_2))])
                

            except KeyError as e:
                print(f"KeyError encountered: {e}")
            except IndexError as e:
                print(f"IndexError encountered: {e}")
            except Exception as e:
                print(f"An unexpected error occurred: {e}")


        ## Outputs ALL possible pairs of genes and entries in the database.
        ## NOTE! This may be VERY 
        ## computationally intensive for larger databases and/or larger protein
        ## candidate families.

        else:
            logger.info('Paring ALL possible combinations of genes in query and database')

            final_list = []
            try:
                count = 0
                for index_1, row_1 in df.iterrows():
                    for index_2, row_2 in df2.iterrows():
                        count += 1
                        seq_1 = row_1['Sequence']
                        seq_2 = row_2['Sequence']
                        final_list.append([(index_1, str(seq_1)), (index_2, str(seq_2))])

            except Exception as e:
                print(f"An error occurred: {e}")
    
            

    except Exception as ex:
        logging.error(f'Error in name_and_sequence_pair function: {ex}')
    logger.debug('Exiting pname_and_sequence_pair function')
    return final_list
