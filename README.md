# Calculating Amino Acid Composition for User-Specified Proteomes 

## Functionality
A Python pipeline to query NCBI with a user-input search information and statistically process the retrieved protein FASTA files. The pipeline provides a user interface to NCBI to pull FASTA sequences by asking the user for what they want to pull. The resulting FASTA sequences are stored in an input file, which is then used to conduct statistical analyses. Amino acid percentages are calculated individually for all FASTA sequences. Average amino acid percentages are calculated between all protein sequences and displayed in a table in an output csv file. These percentages are sorted by rank, so the user can request an amino acid at a specific rank (i.e, most frequent, or 7th most frequent amino acids). Amino acid frequencies are displayed on a barplot, with minimum and maximum frequencies displayed clearly. There is also R code to generate boxplots indicating the frequencies of amino acids and a heat map showcasing a matrix of correlations between proteins.


## Dependencies: 
Built-in Python modules to import:

- sys (manipulates Python runtime environment)
- csv (reads and writes to csv files)
- re (checks if strings match regular expressions)
- datetime (imports current date)
  - date
- math (for calculations)
  - log10
  - floor

Can be installed using pip in the command line:

- numpy (mathematically alters arrays)
- biopython (parses FASTA files, queries NCBI, performs protein analysis)
  - SeqIO
  - Entrez
  - Bio.SeqUtils.ProtParam import ProteinAnalysis
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

R packages to import:

- corrplot


## Inputs
The only file needed to execute this code is **AAPipeline.py**. The R code can be run after the Python code using one of the generated output files (protein.csv). Once the .py file has been downloaded and the user runs the module, the user will be prompted to enter several inputs:

**a.  _Enter email:_**

The user is asked to enter an email address to identify themselves to NCBI in order to access protein sequences from the database.
  
**b.  _Enter NCBI search term:_**

The user is asked to enter the NCBI search term corresponding to the proteins they want to extract. The NCBI database for protein sequence extraction is set to NCBI's Protein database. The user may enter a search term as general as the disease name, or they can list out multiple parameters corresponding to the protein sequences they want to extract, which include but are not limited to: author, accession number, assembly, bioproject, cultivar, division, EC/RN number, feature key, gene name, isolate, journal, organism, strain, exc. A link to NCBI's Protein Advanced Search Builder is included below if the user wants to look at all available search filters:
      https://www.ncbi.nlm.nih.gov/protein/advanced
 
 **c.  _How many protein sequences would you like to extract?_**

The user is asked to enter the number of protein sequences they would like to pull FASTA records for. The first (user-entered integer)  protein sequences will be extracted in FASTA form when the user-specified NCBI search term is executed. 
  
**d.  _Would you like to extract protein sequences from a specified date range? Enter (Y/N):_**

The user can either enter Y (Yes) to extract protein sequences from a specified date range, which will result in further prompting to enter a start date and end date in the form YYYY/MM/DD, or N (No) if they want to work with the default date range, where the start date is set to 2000/01/01 and the end date is set to the current date.


## Outputs
The following files will be outputted once the AAPipeline.py code has been executed:

**1. average_amino_acid_frequencies.csv**

Displays amino acids and corresponding frequencies, averaged across all queried FASTA sequences. The "Rank" column corresponds to the nth most frequent amino acid across sequences. For example, an amino acid with Rank = 4 is the 4th most common amino acid across all queried FASTA sequences.

**2. Figure 1**

Amino acids are displayed on a barplot with average frequencies across all queried FASTA sequences. Percentages for minimum and maximum amino acids are displayed.

**3. proteinSearch.txt**

This text file contains the protein sequences obtained from NCBI using the user-specified parameters in FASTA form. 

**4. protein.csv**

This is a matrix of amino acid frequencies that should be uploaded to R and used to run the code found in CorrelationMxs.R

**5. test.out.csv**



The following images will be outputted once the CorrelationsMx.R code has been executed:

**1. boxplots**

A boxplot will be created for each amino acid, showing the distribution of frequencies for each one. The boxplots can be compared to determine any significant differences in the frequency distributions.

**2. heatmap**

A heatmap denoting the correlations between proteins will be created.


## Test Data

Test data can be generated by running the Python pipeline code. The following responses were entered as user input to generate the files in the Outputs folder.

__Enter email:__ [enter user email]

__Enter NCBI search term:__ sars cov 2

__How many protein sequences would you like to extract?__ 10

__Default date range for protein sequence extraction is from 01/01/2000 -- Current.__
__Would you like to extract protein sequences from a specified date range?__
__Enter (Y/N):__ N

