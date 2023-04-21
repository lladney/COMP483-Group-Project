import sys
import csv
import numpy as np
from Bio import SeqIO
from Bio import Entrez
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import matplotlib.pyplot as plt
import seaborn as sns
import scipy.stats
import re
from datetime import date
today = date.today()

# Define input and output files
#input_file = "C:/Users/param/Downloads/sars-cov-2-prot1.txt"
input_file = open('proteinSearch.txt', 'w')  # open and write sequences to proteins text file 
output_file = "/Users/laraladney/Documents/sars-cov2_amino_acid_frequencies-prot1.csv"
#matrix_output_file = "/Users/laraladney/Downloads/matrix_out.png"
#matrix_values_file = "/Users/laraladney/Downloads/matrix_values.csv"
#matrix_output_file = "C:/Users/param/Downloads/matrix_out.png"
#matrix_values_file = "C:/Users/param/Downloads/matrix_values.csv"

# Protein Sequence Retrieval from NCBI based on search term
Entrez.email = input("Enter email: ")   # user prompted to enter email (tell NCBI who you are to access sequences)

protTerm = input("Enter NCBI search term: ")     # user prompted to enter protein sequence ID
numSeqs = input("How many protein sequences would you like to extract? ") # user prompted to enter # seqs

print("Default date range for protein sequence extraction is from 01/01/2000 - Current.")
print("Would you like to extract protein sequences from a specified date range?")
dateY_N = input("Enter (Y/N): ")
startDate = ""
endDate = ""
if dateY_N == "N" or "n":
    startDate == "2000/01/01" and endDate == today.strftime("%d/%m/%Y")
if dateY_N == "Y" or "y":
    startDate == input("Using format YYYY/MM/DD, enter start date: ")
    endDate == input("Using format YYYY/MM/DD, enter end date: ")

searchResultHandle = Entrez.esearch(db = "protein", term = protTerm, retmax = numSeqs, mindate = startDate, maxdate = endDate)
searchResult = Entrez.read(searchResultHandle)
ids = searchResult["IdList"]

handle = Entrez.efetch(db="protein", id=ids, rettype="fasta", retmode="text")
record = handle.read()

input_file = open('proteinSearch.txt', 'w')
input_file.write(record.rstrip('\n'))
input_file.close()

'''
# CODE TO RETREVE PROTEIN SEQUENCES BASED ON URL
print("Analysis of protein sequences will be performed on search results from uniprot.org.")
url = input("Please enter URL of search results: ")     # user prompted to enter protein sequence ID
all_fastas = requests.get(url).text
fasta_list = re.split(r'\n(?=>)', all_fastas)
[fasta for fasta in fasta_list if '>' in fasta]
fastas = '\n'.join(fasta_list)
input_file.write(fastas)

# CODE TO RETRIEVE SINGLE PROTEIN SEQUENCE FROM NCBI
Entrez.email = input("Enter email: ")   # user prompted to enter email (tell NCBI who you are to access sequences)
seqID = input("Enter protein ID: ")     # user prompted to enter protein sequence ID
handle = Entrez.efetch(db="protein", id=seqID, rettype="gb", retmode="fasta") # handle to fetch protein sequence (user input)
record = SeqIO.read(handle, "genbank")  # use SeqIO parser to read handle (in genbank format)
sequence = str(record.seq)              # save sequence from record as a string
input_file.write(sequence)              # write protein sequence to input file
'''



# need to write to outfile:
# amino acids with the same/most similar frequency across the variants (least likely to mutate)
# amino acids with larges ranges of frequencies across the variants (most likely to mutate)
# table containing shapiro wilk test, skewness, kurtosis

input_file = open('proteinSearch.txt', 'r')

# Define function to extract amino acid sequences from a FASTA file
def extract_amino_acids(input_file):
    amino_acids = []
    for record in SeqIO.parse(input_file, "fasta"):
        sequence = str(record.seq)
        protein = ProteinAnalysis(sequence)
        aa_percentages = protein.get_amino_acids_percent()
        amino_acids.append(aa_percentages)
    return amino_acids

# Define function to write amino acid frequencies and min/max percentages to a CSV file
def write_csv(output_file, data):
    with open(output_file, "w", newline="") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["Rank", "Amino Acid", "Frequency"])
        rank = 1
        for aa, freq in data.items():
#            min_freq = min(amino_acids, key=lambda x: x[aa])[aa]
#            max_freq = max(amino_acids, key=lambda x: x[aa])[aa]
            writer.writerow([rank, aa, freq])
            rank += 1

# Extract amino acid sequences from input file
amino_acids = extract_amino_acids(input_file)

# Close input file
input_file.close()                     

# Calculate total frequencies of each amino acid
total_aa_freqs = {}
for aa in amino_acids:
    for k, v in aa.items():
        if k in total_aa_freqs:
            total_aa_freqs[k] += v
        else:
            total_aa_freqs[k] = v
for k, v in total_aa_freqs.items():
    total_aa_freqs[k] = v/int(numSeqs)
    
total_aa_freqs_sorted = dict(sorted(total_aa_freqs.items(), key = lambda item: item[1], reverse = True))
#rankdict = dict(zip(total_aa_freqs_sorted.keys(), rankdata([-i for i in total_aa_freqs_sorted.values()], method='ordinal')))
            
# Calculate Min/Max %'s of Frequencies
min_percent = min(total_aa_freqs.values())
max_percent = max(total_aa_freqs.values())

# Write amino acid frequencies to output file
write_csv(output_file, total_aa_freqs_sorted)

# Plot amino acid frequencies using seaborn
sns.set_style("whitegrid")
plt.figure(figsize=(10, 6))
sns.barplot(x=list(total_aa_freqs.keys()), y=list(total_aa_freqs.values()), palette="Blues")
plt.title("Amino Acid Frequencies in SARS-CoV-2 Proteins")
plt.xlabel("Amino Acid")
plt.ylabel("Frequency")
plt.text(-0.6, min_percent-0.005, f"Min: {min_percent:.2%}")
plt.text(19.2, max_percent-0.005, f"Max: {max_percent:.2%}")
'''
# Create a matrix of amino acid frequencies
index_table = np.zeros((len(total_aa_freqs), len(amino_acids)))
aa_indices = {aa: i for i, aa in enumerate(total_aa_freqs)}
for aa_freq in amino_acids:
    for aa, freq in aa_freq.items():
        index_table[aa_indices[aa],] += freq

# Write matrix values to output file
with open(matrix_values_file, "w", newline="") as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(["Amino Acid"] + list(range(1, len(amino_acids)+1)))
    for aa, row in zip(total_aa_freqs.keys(), index_table):
        writer.writerow([aa] + list(row))

# Plot amino acid frequency matrix using seaborn
plt.figure(figsize=(10, 8))
sns.heatmap(index_table, cmap="Reds", xticklabels=list(range(1, len(amino_acids)+1)), yticklabels=list(total_aa_freqs.keys()))
plt.title("Amino Acid Frequency Matrix in SARS-CoV-2 Proteins")
plt.xlabel("Protein Sequence Position")
plt.ylabel("Amino Acid")
plt.savefig(matrix_output_file)
'''
plt.show()
