import csv
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
import logging
from Bio.Align import substitution_matrices
import concurrent.futures as futures
from functools import partial

# from protein_sequence_obtainer import name_and_sequence_pair as nm

logger = logging.getLogger(__name__)

def matrix_sw_parall(sequence_pairs):
    (name1, seq1), (name2, seq2) = sequence_pairs
    blosum62 = substitution_matrices.load("BLOSUM62")
    try:
        alignments = pairwise2.align.localdx(seq1, seq2, blosum62)
        best_alignment = alignments[0]
        score = best_alignment[2]

    except Exception as e:
        logger.error(f'Error in matrix_sw_parall function: {e}')

    except KeyboardInterrupt:
        logger.warning("Data processing interrupted by user")

    return {'Name1': name1, 'Name2': name2, 'Score': score}


def raw_sw_parall(match, mismatch, gap_open, gap_extend, sequence_pairs):
    (name1, seq1), (name2, seq2) = sequence_pairs
    try:
        alignments = pairwise2.align.localms(seq1, seq2, match, mismatch, gap_open, gap_extend)
        best_alignment = alignments[0]
        score = best_alignment[2]

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
                with futures.ProcessPoolExecutor() as ex:
                    result = list(ex.map(matrix_sw_parall, sequence_pairs))
                    for i in result:
                        writer.writerow(i)

            elif matrix and not parallel:
                logger.debug('Entered matrix and NOT parallel')
                blosum62 = substitution_matrices.load("BLOSUM62")
                for (name1, seq1), (name2, seq2) in sequence_pairs:
                    try:
                        alignments = pairwise2.align.localdx(seq1, seq2, blosum62)
                        best_alignment = alignments[0]
                        score = best_alignment[2]
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
                for (name1, seq1), (name2, seq2) in sequence_pairs:
                    try:
                        alignments = pairwise2.align.localms(seq1, seq2, match, mismatch, gap_open, gap_extend)
                        best_alignment = alignments[0]
                        score = best_alignment[2]
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

# if __name__ == "__main__":
#     pairs = nm('/Users/william/Documents/Github/Chromoproteins_2024/secret/outputs/PT5/output_PT5_DNAtoProtein.fasta', 'secret/outputs/PT5/output_PT5_sorter.csv', 'data/chromoproteins_uniprot/uniprotkb_chromophore_keyword_KW_0157_AND_reviewed_2024_06_24.fasta')
    
#     output_file = "alignment_scores.csv"
#     smith_waterman_alignment('secret/outputs/test', pairs, 'test', matrix=True)
#     print(f"Alignment scores saved to '{output_file}'")