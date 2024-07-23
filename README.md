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
