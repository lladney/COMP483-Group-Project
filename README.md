# Calculating Amino Acid Composition for SARS-CoV-2 Omicron Variant Proteomes 

## Functionality
A Python pipeline to query NCBI with a user-input search information and statistically process the retrieved protein FASTA files. The pipeline provides a user interface to NCBI to pull FASTA sequences by asking the user for what they want to pull. The resulting FASTA sequences are stored in an input file, which is then used to conduct statistical analyses. Amino acid percentages are calculated individually for all FASTA sequences. Average amino acid percentages are calculated between all protein sequences and displayed in a table in an output csv file. These percentages are sorted by rank, so the user can request an amino acid at a specific rank (i.e, most frequent, or 7th most frequent amino acids). Amino acid frequencies are displayed on a barplot, with minimum and maximum frequencies displayed clearly.

Input
- user input (a search term, number of protein sequences, and an optional date range) to retrieve NCBI protein FASTA sequences

Output:
- csv file containing average amino acid prevalence across the protein sequences
- barplots displaying AAC, with labeled minimum and maximum frequencies.


## Dependencies: 

Built-in Python modules to import:

- sys (manipulates Python runtime environment)
- csv (reads and writes to csv files)
- re
- datetime (imports current date)

Can be installed using pip in the command line:

- numpy (mathematically alters arrays)
- biopython (parses FASTA files and queries NCBI)
- matplotlib (creates 2D graphs and plots)
- seaborn (creates statistical graphics)
- scipy (computes statistics)

```
pip install numpy
pip install biopython
pip install matplotlib
pip install seaborn
pip install scipy
```

## Test Data
Can be found in "Sequences" folder
SARS-CoV-2 Proteomes from Omicron Variants (9)
<img width="872" alt="Screen Shot 2023-03-24 at 5 00 33 PM" src="https://user-images.githubusercontent.com/125703033/227651451-e03a6295-3a8c-4f9e-b9b1-4f5002cbcc7d.png">

Useful Links:
- https://covariants.org/variants/21L.Omicron
- https://www.ncbi.nlm.nih.gov/gene - RefSeq for 11 Omicron sars-cov 2 genomes
- https://protfasta.readthedocs.io/en/latest/read_fasta.html (remove invalid sequences from FASTA)
- https://github.com/Wytamma/GISAIDR (retrieve files from GISAID)
- https://www.ncbi.nlm.nih.gov/nuccore/NC_045512.2?from=266&to=21555&report=fasta
- https://github.com/biopython/biopython/blob/master/Bio/SeqUtils/ProtParam.py
