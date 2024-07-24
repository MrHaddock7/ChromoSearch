# ChromoSearch

This is a pipeline that can be used to find matches of chromoproteins or pigment creating enzymes in a genome. This could be used as a first step when looking for the source of a color in a organism.

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

![Visualisation of pipeline](pictures/pipeline.drawio.svg)
