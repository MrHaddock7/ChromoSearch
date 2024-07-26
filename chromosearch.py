import argparse
import logging
import os
import tempfile

from scripts.protein_sequence_obtainer import name_and_sequence_pair as nm
from scripts.smith_waterman import smith_waterman_alignment as sm
from scripts.protein_search import protein_blastp_search as pbs
from scripts.sorter import csv_sorter
from scripts.sorter import csv_sorter2
from scripts.DNAtoProtein_prodigal import run_prodigal as DNAtoProtein

def main(fasta_path, output_path, gene, database='/Users/klonk/Desktop/uniprotkb_e_coli_photo_AND_reviewed_tru_2024_07_26.fasta', process=True, save_intermediates=True, parallel=True, matrix=True, match=3, mismatch=-1, gap_open=-10, gap_extend=-4, blastpnsw=True):
    
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

    if process:
        print(f'DNAtoProtein: started...')
        DNAtoProtein(fasta_path, output_dir, gene)
        print(f'DNAtoProtein: finished')

    # if save_intermediates:
    temp_protein_search = os.path.join(output_dir)
    temp_SW_csv = os.path.join(output_dir)

    try:
        # print(f'pbs: started...')
        # pbs(f'{output_dir}/output_{gene}_DNAtoProtein.fasta', gene, temp_protein_search)
        # print(f'pbs: finished')

        print(f'csv: started...')
        csv_sorter(os.path.join(temp_protein_search, f'output_{gene}_protein_search.csv'), gene, temp_SW_csv)
        print(f'csv: finished')

        print(f'smith waterman + name_and_sequence_pair started...')
        sequence_pairs = nm(f'{output_dir}/output_{gene}_DNAtoProtein.fasta', os.path.join(temp_SW_csv, f'output_{gene}_sorted_pBLAST.csv'), input_database_fasta=f'{database}', blastpsw=blastpnsw)
        print(len(sequence_pairs))
        sm(output_dir, gene_name=gene, sequence_pairs=sequence_pairs, parallel=parallel, matrix=matrix, match=match, mismatch=mismatch, gap_open=gap_open, gap_extend=gap_extend)
        print(f'smith waterman + name_and_sequence_pair finished')

        csv_sorter2(f'{output_dir}/output_{gene}_smith_waterman.csv', gene, output_dir)

    finally:
        if not save_intermediates:
            os.remove(os.path.join(temp_protein_search, f'output_{gene}_protein_search.csv'))
            os.remove(os.path.join(temp_SW_csv, f'output_{gene}_smith_waterman.csv'))

    logger.info(f'Finished processing the {gene} gene')

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process a single genomic data file and perform various bioinformatics tasks.")
    parser.add_argument("fasta_file", help="Path to the fasta file with the whole genome for the strain")
    parser.add_argument("output_path", help="Path to where to save the output files")
    parser.add_argument("gene", help="Name of the gene to process")
    parser.add_argument("-db", "--database", default='/Users/klonk/Desktop/uniprotkb_e_coli_photo_AND_reviewed_tru_2024_07_26.fasta', help="Path to the chromoprotein database")
    parser.add_argument("-p", "--parallel", action="store_false", help="If you want to disable the parallelization")
    parser.add_argument("-M", "--matrix", action="store_false", help="If you want to disable BLOSUM62 matrix and use standard scores")
    parser.add_argument("--match", type=int, default=3, help="Score for a match")
    parser.add_argument("--mismatch", type=int, default=-1, help="Penalty for a mismatch")
    parser.add_argument("--gap_open", type=int, default=-10, help="Gap opening penalty")
    parser.add_argument("--gap_extend", type=int, default=-4, help="Gap extension penalty")
    parser.add_argument("-s", "--save_intermediates", action="store_true", help="Flag to save intermediate files")
    parser.add_argument("-P", "--process", action="store_false", help="Flag to skip the DNA to Protein processing step")

    args = parser.parse_args()

    fasta_file = args.fasta_file
    output_path = args.output_path
    gene = args.gene
    database = args.database
    save_intermediates = args.save_intermediates
    parallel = args.parallel
    matrix = args.matrix
    match = args.match
    mismatch = args.mismatch
    gap_open = args.gap_open
    gap_extend = args.gap_extend
    process = args.process

    main(fasta_file, output_path, gene, database, process, save_intermediates, parallel, matrix, match, mismatch, gap_open, gap_extend)