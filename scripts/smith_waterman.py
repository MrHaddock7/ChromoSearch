import csv
import logging
import concurrent.futures as futures
from functools import partial
from itertools import islice
from Bio.Align import PairwiseAligner
from Bio.Align import substitution_matrices

logger = logging.getLogger(__name__)

def sequence_pairs_smith_waterman(match, mismatch, gap_open, gap_extend, matrix, sequence_pairs):
    """Basic funtion that takes a list of sequence_pairs and returns their names along with their scores. Used by smith_waterman_alignment().

    Args:
        match (int): Score for match
        mismatch (_type_): Score for mismatch
        gap_open (_type_): Penalty for gap opening
        gap_extend (_type_): Penalty for gap extension
        sequence_pairs (_type_): List of sequence pairs, format [[(name1, seq1), (name2, seq_2)], ...]
    Returns:
        list: List of dictionaries, each with the format {'Name1': name1, 'Name2': name2, 'Score': score}.
    """   

    # Basic alignment setup
    aligner = PairwiseAligner()
    aligner.mode = 'local'
    aligner.match_score = match        
    aligner.mismatch_score = mismatch    
    aligner.open_gap_score = gap_open   
    aligner.extend_gap_score = gap_extend 

    # TO IMPLEMENT: SUPPORT OTHER MATRICES 
    if matrix is not None:
        aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")

    # List of dirs to return
    results_list = []

    for pair in sequence_pairs:
        (name1, seq1), (name2, seq2) = pair
        try:
            best_alignment = aligner.align(seq1, seq2)[0]
            score = best_alignment.score

        except Exception as e:
            logger.error(f'Error in sequence_pair_smith_waterman(), for sequences {name1}, {name2}: {e}')

        except KeyboardInterrupt:
            logger.warning("Data processing interrupted by user.")

        # This is the most inefficient data structure ever created
        results_list.append({'Name1': name1, 'Name2': name2, 'Score': score})

    return results_list


# NOTE: Changed from one pair = one process to batch processing
def batch_sequence_pairs(sequence_pairs, threads):
    """Creates the sequence pairs for the batch process

    Args:
        sequence_pairs (list): List of sequence pairs, format [[(name1, seq1), (name2, seq_2)], ...]
        threads (int): Number of threads.

    Returns:
        list: List of sequence pairs batches.
    """

    CHUNK_SIZE = len(sequence_pairs) // threads
    
    # Uses generator for memory efficiency
    it = iter(sequence_pairs)
    while True:
        batch = list(islice(it, CHUNK_SIZE))
        if not batch:
            return
        yield batch

def write_smith_waterman_results(output, gene_name, results_dictonaries):
    """Function to write the results of the Waterman-Smith alignment into a .csv file.

    Args:
        output (str): Output file location.
        gene_name (str): Naming prefix for the results.
        results_dir (): List of dictionaries to write out.
    """

    with open(f'{output}/output_{gene_name}_smith_waterman.csv', 'w', newline='') as csvfile:
        
        try:
            fieldnames = ['Name1', 'Name2', 'Score']
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
            writer.writeheader()

            # Write all results
            for result in results_dictonaries:
                writer.writerow(result)
        
        except Exception as e:
            logger.error(f'Error in function write_smith_waterman_results(): {e}')

    return

def smith_waterman_alignment(output, sequence_pairs, gene_name, match=3, mismatch=-1, gap_open=-10, gap_extend=-4, matrix=True, threads= 1):

    logger.debug('Entering smith_waterman_alignment function')

    result_to_write = []
    try:
        # # Single-threaded mode
        # if threads == 1:
        #     logger.debug('Entered single-threaded mode...')
        #     result_to_write = sequence_pairs_smith_waterman(match, mismatch, gap_open, gap_extend, sequence_pairs, matrix)
            
        # else:
        # logger.debug('Entered multi-threaded mode...')

        # Set the first arguments for the function as static, and map to the batch sequence pairs
        partial_sequence_pair_smith_waterman = partial(sequence_pairs_smith_waterman, match, mismatch, gap_open, gap_extend, matrix)
        with futures.ProcessPoolExecutor() as ex:
            result = list(ex.map(partial_sequence_pair_smith_waterman, batch_sequence_pairs(sequence_pairs, threads)))
        
        # flatten result list of lists of dictionaries
        # Could retain use of a generator if memory is a problem - Unlikely
        result_to_write = [item for sublist in result for item in sublist]

    except Exception as e:
        logger.error(f'Error in smith_waterman_alignment function: {e}')

    except KeyboardInterrupt:
        logger.warning("Data processing interrupted by user")

    logger.debug("Waterman-Smith finished, writing results...")

    # Write results to file
    write_smith_waterman_results(output, gene_name, result_to_write)

    return
