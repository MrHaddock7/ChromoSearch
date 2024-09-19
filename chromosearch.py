## NOTE! This script is meant to be run through the terminal

import argparse
import logging
import os
import tempfile
import pandas as pd

from scripts.initialization_scripts import check_requirements
from scripts.protein_sequence_obtainer import name_and_sequence_pair as nm
from scripts.smith_waterman import smith_waterman_alignment as sm
from scripts.protein_search import protein_blastp_search as pbs
from scripts.sorter import csv_sorter
from scripts.DNAtoProtein_prodigal import run_prodigal as DNAtoProtein
from scripts.statistical_analysis import statistics_calculation
from scripts.initialization_scripts import suppress_output

## Thanos' code

from scripts.characterize_proteins import dereplicate_highest_score
from scripts.characterize_proteins import calculate_mass_length

## main function


def main(
    fasta_path,
    output_path,
    gene,
    database,
    process=True,
    save_intermediates=True,
    threads=1,
    matrix=True,
    match=3,
    mismatch=-1,
    gap_open=-10,
    gap_extend=-4,
    blastpnsw=True,
    mass_n_length=True,
    multiple_test_correction="fdr_bh",
):

    logging.basicConfig(
        level=logging.DEBUG,
        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
        handlers=[logging.FileHandler("project.log")],
    )

    logger = logging.getLogger(__name__)

    logger.info("Starting main function")
    logger.info(f"Processing the {gene} gene")
    logger.info(f"{database}")

    output_dir = f"{output_path}/{gene}"
    temp_output = f"temp/{gene}"
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        logger.info(f"Created directory: {output_dir}")

    if not os.path.exists(temp_output):
        os.makedirs(temp_output)
        logger.info(f"Created directory: {temp_output}")

    # =================================================================
    # Check requirements before running the pipeline
    # =================================================================

    # Dictionary used for checking requirements
    # packages are used by
    requirements = {
        "packages": [
            ("Bio", "1.84"),  # biopython
            ("matplotlib", "3.9.2"),
            ("numpy", "2.1.0"),
            ("pandas", "2.2.2"),
            ("scipy", "1.14.1"),
            ("seaborn", "0.13.2"),
            ("statsmodels", "0.14.2"),
        ],
        "python_version": "3.12.5",
        "prodigal_version": "2.6.3",
        "blast_version": "2.16.0",
    }

    check_requirements(requirements)

    # if save_intermediates:
    temp_protein_search = os.path.join(temp_output)
    temp_SW_csv = os.path.join(temp_output)

    try:

        ## Initiation of the pipeline

        if process:
            print(f"Identifying candidate proteins in DNA: started...")
            DNAtoProtein(fasta_path, output_dir, gene)
            print(f"Identifying candidate proteins in DNA: complete")

        print("Running blastP search: started...")
        pbs(
            f"{output_dir}/output_{gene}_DNAtoProtein.fasta",
            gene,
            temp_protein_search,
            input_database=f"{database}",
            threads=threads,
        )

        print(f"Running blastP search: complete")

        print(f"Removing hits with high E-values: started")
        csv_sorter(
            input_csv=f"{temp_protein_search}/output_{gene}_protein_search.csv",
            genome=gene,
            output=temp_protein_search,
            sort_value_metric="evalue",
            cut_off_value=float(0.05),
            name_output="sorted_pBLAST",
        )
        print(f"Removing hits with high E-values: complete")

        print(f"smith waterman + name_and_sequence_pair started...")
        sequence_pairs = nm(
            f"{output_dir}/output_{gene}_DNAtoProtein.fasta",
            f"{temp_protein_search}/output_{gene}_sorted_pBLAST.csv",
            input_database_fasta=f"{database}",
            blastpsw=blastpnsw,
        )
        print(
            f"Performing the Smith-Waterman algorithm on {len(sequence_pairs)} sequence pairs..."
        )

        sm(
            temp_SW_csv,
            gene_name=gene,
            sequence_pairs=sequence_pairs,
            threads=threads,
            matrix=matrix,
            match=match,
            mismatch=mismatch,
            gap_open=gap_open,
            gap_extend=gap_extend,
        )
        print(f"smith waterman + name_and_sequence_pair finished")

        csv_sorter(
            f"{temp_SW_csv}/output_{gene}_smith_waterman.csv",
            gene,
            temp_SW_csv,
            only_sort=True,
            sort_value_metric="Score",
            name_output="sorted_alignment",
        )

        ## Implementation of normalization code
        results_with_mass_and_length = 0
        if mass_n_length:
            print("Calculation of mass and length of candidate proteins: started...")
            dereplicated_results = dereplicate_highest_score(
                f"{temp_SW_csv}/output_{gene}_sorted_alignment.csv"
            )
            results_with_mass_and_length = calculate_mass_length(
                f"{output_dir}/output_{gene}_DNAtoProtein.fasta",
                dereplicated_results,
                f"{temp_protein_search}/output_{gene}_sorted_pBLAST.csv",
            )

            print("Calculation of mass and length of candidate proteins: finished")

        # Statistical analysis - thanos
        # ==================================================================================================================

        # Create directory for outputs of statistical analysis (plots)
        statistics_directory = f"{output_dir}/Statistical_analysis/"

        os.makedirs(statistics_directory, exist_ok=True)

        # TODO: add support for changing plot_dpi through the command line
        final_results_dataframe = statistics_calculation(
            results_with_mass_and_length,
            statistics_directory,
            multiple_test_correction,
        )

        # final corrections to the dataframe
        final_results_dataframe.rename(
            columns={"Name1": "Genome_entry_id", "Name2": "Database_hit_id"},
            inplace=True,
        )

        # only save the final results, with statistics
        final_results_dataframe.to_csv(
            f"{output_dir}/chromosearch_{gene}_final_results.csv"
        )

    finally:
        if not save_intermediates:
            os.remove(
                os.path.join(temp_protein_search, f"output_{gene}_protein_search.csv")
            )
            os.remove(os.path.join(temp_SW_csv, f"output_{gene}_smith_waterman.csv"))

    logger.info(f"Finished processing the {gene} gene")


if __name__ == "__main__":

    ## The code below is only valid if the chromosearch-function is run via the terminal.

    parser = argparse.ArgumentParser(
        description="Process a single genomic data file and perform various bioinformatics tasks."
    )
    parser.add_argument(
        "-i",
        "--input-genome",
        help="Path to the fasta file for the bacterial genome.",
        required=True,
    )
    parser.add_argument(
        "-p", "--prefix", required=True, help="Naming prefix for output files."
    )
    parser.add_argument(
        "-db",
        "--database",
        default="databases/chromoproteins_uniprot/uniprotkb_chromophore_keyword_KW_0157_AND_reviewed_2024_06_24.fasta",
        help="Path to the chromoprotein database",
    )
    parser.add_argument(
        "-t",
        "--threads",
        type=int,
        default=1,
        help="Number of threads available to the pipeline. Set to 0 or negative numbers to use all available cores",
    )
    parser.add_argument(
        "-M",
        "--matrix",
        action="store_false",
        help="If you want to disable BLOSUM62 matrix and use standard scores",
    )
    parser.add_argument("--match", type=int, default=3, help="Score for a match")
    parser.add_argument(
        "--mismatch", type=int, default=-1, help="Penalty for a mismatch"
    )
    parser.add_argument(
        "--gap_open",
        type=int,
        default=-10,
        help="Gap opening penalty for the Waterman-Smith alignments.",
    )
    parser.add_argument(
        "--gap_extend",
        type=int,
        default=-4,
        help="Gap extension penalty for the Waterman-Smith alignments.",
    )
    parser.add_argument(
        "-s",
        "--save_intermediates",
        action="store_true",
        help="Flag to save intermediate files",
    )
    parser.add_argument(
        "-P",
        "--process",
        action="store_false",
        help="Flag to skip the DNA to Protein processing step",
    )
    parser.add_argument(
        "-bpsw",
        "--blastpandsmithwaterman",
        action="store_false",
        help="Flag to turn off Smith-Waterman algorithm based on blastp results, and instead perform Smith-Waterman on all possible candidate protein-database protein combinations. NOTE! This can be VERY computationally intensive for larger sequences and/or databases.",
    )
    parser.add_argument(
        "-q",
        "--quiet",
        action="store_true",
        help="Quiets the text-outputs of the ChromoSearch.",
    )
    # Argument for the multiple test correction method used by statistics.py
    parser.add_argument(
        "--mutliple-correction",
        type=str,
        default="fdr_bh",
        choices=[
            "bonferroni",
            "sidak",
            "holm-sidak",
            "holm",
            "simes-hochberg",
            "hommel",
            "fdr_bh",
            "fdr_by",
            "fdr_tsbh",
            "fdr_tsbky",
        ],
        help='Method used to correct p-values for multiple testing using the statsmodels.stats.multitest module. Available methods: "bonferroni", "sidak", "holm-sidak", "holm", "simes-hochberg", "hommel", "fdr_bh", "fdr_by", "fdr_tsbh", "fdr_tsbky" Default:  Benjamini-Hochberg.',
    )
    parser.add_argument("output_path", help="Path to where to save the output files")

    args = parser.parse_args()

    fasta_path_argument = args.input_genome
    output_path_argument = args.output_path
    gene_argument = args.prefix
    database_argument = args.database
    save_intermediates_argument = args.save_intermediates
    matrix_argument = args.matrix
    match_argument = args.match
    mismatch_argument = args.mismatch
    gap_open_argument = args.gap_open
    gap_extend_argument = args.gap_extend
    process_argument = args.process
    blastpnsw_argument = args.blastpandsmithwaterman
    mutliple_test_correction = args.mutliple_correction

    # check options for threads arg
    # assign all available cpus
    if args.threads <= 0:
        threads = os.cpu_count()
    elif args.threads > os.cpu_count():
        raise ValueError(
            "Number of threads specified exceeds those available in the system. Please specify a lower count, or run in single-threaded mode."
        )
    else:
        threads = args.threads

    if threads == 1:
        print(
            "You have chosen to run the pipeline using only 1 thread. This might take some time...\n"
        )

    # Use decorated main to suppress output, other than errors
    decorated_with_suppress_main = suppress_output(args.quiet)(main)

    decorated_with_suppress_main(
        fasta_path=fasta_path_argument,
        output_path=output_path_argument,
        gene=gene_argument,
        database=database_argument,
        process=process_argument,
        save_intermediates=save_intermediates_argument,
        threads=threads,
        matrix=matrix_argument,
        match=match_argument,
        mismatch=mismatch_argument,
        gap_open=gap_open_argument,
        gap_extend=gap_extend_argument,
        blastpnsw=blastpnsw_argument,
        multiple_test_correction=mutliple_test_correction,
    )
