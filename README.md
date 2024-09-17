# ChromoSearch

Our pipeline offers a streamlined approach for identifying chromoproteins by analyzing genomic data. It starts with sequenced genomic data in .fasta format and uses Prodigal to identify potential protein-coding sequences, which are translated into putative proteins. These are then compared against databases of known chromoproteins or pigment-producing enzymes using pBLAST, followed by a more detailed Smith-Waterman alignment for high-scoring matches. The molecular weight, length, and a normalized score are calculated for each protein to assess their plausibility, and a Gumbel distribution is used to model rare hits. The pipeline automates these steps, providing an efficient way to process next-generation sequencing data and identify potential chromoproteins for further validation.

This approach offers a faster, cost-effective alternative to traditional protein identification methods by combining genomic data analysis with automated protein sequence prediction and comparison.

# Installation

## Conda

The recommended method of installation is through the environment.yml file found in the base directory. Conda will automatically solve all requirements as necessary.

To install:

```
git clone https://github.com/MrHaddock7/ChromoSearch.git
conda env create -f environment.yml
```

## Manual

If the necessary packages are installed and available in the working environment, then the repository can be cloned as is.

ChromoSearch has been tested up to now with the following package versions:

```
python: 3.12.5
pandas: 2.2.2
biopython: 1.84
blast: 2.16.0
prodigal: 2.6.3
```


## How to use

If you want to, you could start a virtual environment to run the pipeline on.

```
python3 -m venv venv
. venv/bin/activate # on Linux, MacOS; or
. venv\Scripts\activate # on Windows
pip install -r requirements.txt
```

### Example use:

```
git clone https://github.com/MrHaddock7/ChromoSearch.git
cd chromosearch
python3 chromosearch.py path/to/genome.fasta path/for/output_folder name_of_gene_or_job
```

### For help flags:

```
python3 chromosearch.py -h
```

## How it works

The input for the pipeline is a .fasta file consiting of the genome you have sequenced. The pipeline will take this and find all protein coding sequences and translate them into protein sequences.

We will then run a pBLAST search of of all the proteins and compare them to a database of either known chromoproteins, or pigment creating enzymes. We will append these results to a csv file and sort and filter the hits according to e-value < 0.05.

We will use these results and do a smith-waterman alignment on all good hits from the pBLAST run. You can then compare the results from both csv files and see if you have any potential candidates.

![Visualisation of pipeline](pictures/pipeline4.drawio.svg)

## How to interpret the output data

The output data is provided as a csv file found in the output directory as specified by the user, ordered by the normalized score for each hit, in descending order. Below, each of columns are explained.

-(prodigal_protein), Contains the identified candidate protein ID as provided by prodigal (change Name1) 

-(database_protein), Contains the protein ID for the protein as provided by from the local database, based on a relevant UniProt query. (change Name2).

-(score), Contains a score as provided by the Smith-Waterman local alignment.

-(length) Contains the calculated length of the identified candidate protein.

-(normalized_score) Contains a normalized score for the identified candidate protein, simply defined as the score divided by the length of the protein. 

-(mass) Contains the estimated mass of the protein based on the protein sequence.

## An interactive interface for ChromoSearch

For a more user-friendly experience with the ChromoSearch pipeline, you can use our interactive interface by running the main_interface.py script. While this interface offers slightly less options and many settings are fixed compared to the command-line version, it is perfect for users who are less familiar with command-line operations or simply prefer a more intuitive and easy-to-navigate option. 
