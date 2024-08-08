## NOTE! This script is meant to be run through the terminal

import argparse
import logging
import os
import tempfile

from scripts.protein_sequence_obtainer import name_and_sequence_pair as nm
from scripts.smith_waterman import smith_waterman_alignment as sm
from scripts.protein_search import protein_blastp_search as pbs
from scripts.sorter import csv_sorter_final
from scripts.sorter import csv_sorter
from scripts.sorter import csv_sorter2
from scripts.DNAtoProtein_prodigal import run_prodigal as DNAtoProtein

## Thanos' code

from scripts.characterize_proteins import dereplicate_highest_score
from scripts.characterize_proteins import calculate_mass_length

def main(fasta_path, 
         output_path, 
         gene, 
         database, 
         process=True, 
         save_intermediates=True, 
         parallel=True, 
         matrix=True, 
         match=3, 
         mismatch=-1, 
         gap_open=-10, 
         gap_extend=-4, 
         blastpnsw=True, 
         quiet_mode=False,
         mass_n_length=True):
    
    logging.basicConfig(
        level=logging.DEBUG,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        handlers=[logging.FileHandler("project.log")]
    )

    logger = logging.getLogger(__name__)   
   
    logger.info('Starting main function')
    logger.info(f'Processing the {gene} gene')
    logger.info(f'{database}')

    output_dir = f'{output_path}/{gene}'
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        logger.info(f'Created directory: {output_dir}')


    # if save_intermediates:
    temp_protein_search = os.path.join(output_dir)
    temp_SW_csv = os.path.join(output_dir)

    def print_quiet_mode(message):
        if not quiet_mode:
            print(message)

    try:

        ## Suggestion: Include all print messages in -q flag

        if process:
            print_quiet_mode(f'DNAtoProtein: started...')
            DNAtoProtein(fasta_path, output_dir, gene)
            print_quiet_mode(f'DNAtoProtein: finished')


        print_quiet_mode('pbs: started...')
        pbs(f'{output_dir}/output_{gene}_DNAtoProtein.fasta', gene, temp_protein_search, input_database=f'{database}')
        print_quiet_mode(f'pbs: finished')

        print_quiet_mode(f'csv: started...')
        # csv_sorter(os.path.join(temp_protein_search, f'output_{gene}_protein_search.csv'), gene, temp_SW_csv)
        csv_sorter_final(f'{output_dir}/output_{gene}_protein_search.csv', gene, temp_SW_csv, 'evalue', cut_off_value=0.05, name_output='sorted_pBLAST')
        print_quiet_mode(f'csv: finished')

        print_quiet_mode(f'smith waterman + name_and_sequence_pair started...')
        sequence_pairs = nm(f'{output_dir}/output_{gene}_DNAtoProtein.fasta', os.path.join(temp_SW_csv, f'output_{gene}_sorted_pBLAST.csv'), input_database_fasta=f'{database}.fasta', blastpsw=blastpnsw)
        
        print_quiet_mode(f'Performing the Smith-Waterman algorithm on {len(sequence_pairs)} sequence pairs...')

        sm(output_dir, gene_name=gene, sequence_pairs=sequence_pairs, parallel=parallel, matrix=matrix, match=match, mismatch=mismatch, gap_open=gap_open, gap_extend=gap_extend)
        print_quiet_mode(f'smith waterman + name_and_sequence_pair finished')

        # csv_sorter2(f'{output_dir}/output_{gene}_smith_waterman.csv', gene, output_dir)
        csv_sorter_final(f'{output_dir}/output_{gene}_smith_waterman.csv', gene, temp_SW_csv, only_sort=True, sort_value_metric='Score', name_output='sorted_alignment')

        ## Implementation of Thanos' code

        if mass_n_length:
            print('mass_n_length started...')
            dereplicated_results = dereplicate_highest_score(f'{output_dir}/output_{gene}_sorted_alignment.csv')
            calculate_mass_length(f'{output_dir}/output_{gene}_DNAtoProtein.fasta', dereplicated_results, gene, output_dir)
            print('mass_n_length finished')
    finally:
        if not save_intermediates:
            os.remove(os.path.join(temp_protein_search, f'output_{gene}_protein_search.csv'))
            os.remove(os.path.join(temp_SW_csv, f'output_{gene}_smith_waterman.csv'))

    logger.info(f'Finished processing the {gene} gene')

if __name__ == "__main__":
    ## Suggestion: We should rewrite some of the helper strings

    parser = argparse.ArgumentParser(description="Process a single genomic data file and perform various bioinformatics tasks.")
    parser.add_argument("fasta_file", help="Path to the fasta file with the whole genome for the strain")
    parser.add_argument("output_path", help="Path to where to save the output files")
    parser.add_argument("gene", help="Name of the gene to process")
    parser.add_argument("-db", "--database", default='databases/chromoproteins_uniprot/uniprotkb_chromophore_keyword_KW_0157_AND_reviewed_2024_06_24', help="Path to the chromoprotein database")
    parser.add_argument("-p", "--parallel", action="store_false", help="If you want to disable the parallelization")
    parser.add_argument("-M", "--matrix", action="store_false", help="If you want to disable BLOSUM62 matrix and use standard scores")
    parser.add_argument("--match", type=int, default=3, help="Score for a match")
    parser.add_argument("--mismatch", type=int, default=-1, help="Penalty for a mismatch")
    parser.add_argument("--gap_open", type=int, default=-10, help="Gap opening penalty")
    parser.add_argument("--gap_extend", type=int, default=-4, help="Gap extension penalty")
    parser.add_argument("-s", "--save_intermediates", action="store_true", help="Flag to save intermediate files")
    parser.add_argument("-P", "--process", action="store_false", help="Flag to skip the DNA to Protein processing step")
    parser.add_argument("-bpsw", "--blastpandsmithwaterman", action="store_false", help="Flag to turn off Smith-Waterman algorithm based on blastp results, and instead perform Smith-Waterman on all possible candidate protein-database protein combinations. NOTE! This can be VERY computationally intensive for larger sequences and/or databases.")
    parser.add_argument("-q", "--quiet", action="store_true", help="Quiets the text-outputs of the ChromoSearch")

    args = parser.parse_args()

    fasta_path_argument = args.fasta_file
    output_path_argument = args.output_path
    gene_argument = args.gene
    database_argument = args.database
    save_intermediates_argument = args.save_intermediates
    parallel_argument = args.parallel
    matrix_argument = args.matrix
    match_argument = args.match
    mismatch_argument = args.mismatch
    gap_open_argument = args.gap_open
    gap_extend_argument = args.gap_extend
    process_argument = args.process
    blastpnsw_argument = args.blastpandsmithwaterman
    quiet_argument = args.quiet



    main(fasta_path=fasta_path_argument, 
         output_path=output_path_argument, 
         gene=gene_argument, 
         database=database_argument, 
         process=process_argument, save_intermediates=save_intermediates_argument, 
         parallel=parallel_argument, matrix=matrix_argument, 
         match=match_argument, 
         mismatch=mismatch_argument, 
         gap_open=gap_open_argument, 
         gap_extend=gap_extend_argument, 
         blastpnsw=blastpnsw_argument,
         quiet_mode=quiet_argument)