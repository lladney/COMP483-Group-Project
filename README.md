# COMP483-Group-Project 

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

##Test Data
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
