import csv
import logging
from Bio.Align import substitution_matrices
import concurrent.futures as futures
from functools import partial
from Bio.Align import PairwiseAligner

logger = logging.getLogger(__name__)

def matrix_sw_parall(match, mismatch, gap_open, gap_extend, sequence_pairs):

    ## Used when you want to parallelize the alignment with BLOSUM62 matrix

    aligner = PairwiseAligner()
    aligner.mode = 'local'
    aligner.match_score = match        # Match score
    aligner.mismatch_score = mismatch    # Mismatch penalty
    aligner.open_gap_score = gap_open   # Gap opening penalty
    aligner.extend_gap_score = gap_extend # Gap extension penalty
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")

    (name1, seq1), (name2, seq2) = sequence_pairs
    try:
        best_alignment = aligner.align(seq1, seq2)[0]
        score = best_alignment.score

    except Exception as e:
        logger.error(f'Error in matrix_sw_parall function: {e}')

    except KeyboardInterrupt:
        logger.warning("Data processing interrupted by user")

    return {'Name1': name1, 'Name2': name2, 'Score': score}


def raw_sw_parall(match, mismatch, gap_open, gap_extend, sequence_pairs):

    ## Used when you want to parallelize the alignment using scoring system

    aligner = PairwiseAligner()
    aligner.mode = 'local'
    aligner.match_score = match        # Match score
    aligner.mismatch_score = mismatch    # Mismatch penalty
    aligner.open_gap_score = gap_open   # Gap opening penalty
    aligner.extend_gap_score = gap_extend # Gap extension penalty

    (name1, seq1), (name2, seq2) = sequence_pairs
    try:
        best_alignment = aligner.align(seq1, seq2)[0]
        score = best_alignment.score

    except Exception as e:
        logger.error(f'Error in raw_sw_parall function: {e}')

    except KeyboardInterrupt:
        logger.warning("Data processing interrupted by user")

    return {'Name1': name1, 'Name2': name2, 'Score': score}


def smith_waterman_alignment(output, sequence_pairs, gene_name, parallel=True, match=3, mismatch=-1, gap_open=-10, gap_extend=-4, matrix=True):
    logger.debug('Entering smith_waterman_alignment function')
    try:
        with open(f'{output}/output_{gene_name}_smith_waterman.csv', 'w', newline='') as csvfile:

            fieldnames = ['Name1', 'Name2', 'Score']
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
            writer.writeheader()

            if matrix and parallel:
                logger.debug('Entered matrix and parallel')
                partial_matrix_sw_parall = partial(matrix_sw_parall, match, mismatch, gap_open, gap_extend)
                with futures.ProcessPoolExecutor() as ex:
                    result = list(ex.map(partial_matrix_sw_parall, sequence_pairs))
                    for i in result:
                        logger.debug('Writing matrix and parallel')
                        writer.writerow(i)

            elif matrix and not parallel:
                logger.debug('Entered matrix and NOT parallel')
                
                aligner = PairwiseAligner()
                aligner.mode = 'local'
                aligner.match_score = match        # Match score
                aligner.mismatch_score = mismatch    # Mismatch penalty
                aligner.open_gap_score = gap_open   # Gap opening penalty
                aligner.extend_gap_score = gap_extend # Gap extension penalty
                aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")

                for (name1, seq1), (name2, seq2) in sequence_pairs:
                    try:
                        best_alignment = aligner.align(seq1, seq2)[0]
                        score = best_alignment.score
                    except Exception as e:
                        logger.error(f'Error in smith_waterman_alignment function: {e}')
                    except KeyboardInterrupt:
                        logger.warning("Data processing interrupted by user")
                    
                    writer.writerow({'Name1': name1, 'Name2': name2, 'Score': score})

            elif parallel and not matrix:
                logger.debug('Entered parallel and NOT matrix')
                partial_raw_sw_parall = partial(raw_sw_parall, match, mismatch, gap_open, gap_extend)
                with futures.ProcessPoolExecutor() as ex:
                    result = list(ex.map(partial_raw_sw_parall, sequence_pairs))
                for i in result:
                    print(i)
                    writer.writerow(i)
            else:
                logger.debug('Entered NOT matrix and NOT parallel')
                
                aligner = PairwiseAligner()
                aligner.mode = 'local'
                aligner.match_score = match        # Match score
                aligner.mismatch_score = mismatch    # Mismatch penalty
                aligner.open_gap_score = gap_open   # Gap opening penalty
                aligner.extend_gap_score = gap_extend # Gap extension penalty

                for (name1, seq1), (name2, seq2) in sequence_pairs:
                    try:
                        best_alignment = aligner.align(seq1, seq2)[0]
                        score = best_alignment.score
                    except Exception as e:
                        logger.error(f'Error in smith_waterman_alignment function: {e}')
                    except KeyboardInterrupt:
                        logger.warning("Data processing interrupted by user")
                    
                    writer.writerow({'Name1': name1, 'Name2': name2, 'Score': score})
                    
        logger.info(f'Result from smith_waterman_alignment function saved in: output_{gene_name}_smith_waterman.csv')

    except Exception as e:
        logger.error(f'Error in smith_waterman_alignment function: {e}')

    except KeyboardInterrupt:
        logger.warning("Data processing interrupted by user")

    logger.debug('Exiting smith_waterman_alignment function')