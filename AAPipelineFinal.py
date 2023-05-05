# Import necessary packages
import sys                                              # module contains methods and variables for modifying Python's runtime environment
import csv                                              # module implements classes to read and write tabular data in csv format
import numpy as np                                      # module for array creation and working with numerical data
from Bio import SeqIO                                   # module functioning as an interface to input and output fasta format files
from Bio import Entrez                                  # module to search NCBI for protein sequences with user-specified parameters
from Bio.SeqUtils.ProtParam import ProteinAnalysis      # module for analysis of protein sequences
import matplotlib.pyplot as plt                         # module to plot numerical data
import seaborn as sns                                   # module for visualization and exploratory analysis, based on matplotlib
import scipy.stats                                      # module to conduct statistical tests
import re                                               # module for regular expression matching operations
import pandas as pd                                     # module for building dataframes
from math import log10, floor                           # module used to output user-specified number of significant figures in output data
from datetime import date                               # module to pull current dates
                                   

# Define input and output files
input_file = open('proteinSearch.txt', 'w')             # open and write sequences to proteins text file 
output_file = "/Users/laraladney/Documents/sars-cov2_amino_acid_frequencies-prot1.csv" # path to output file 1
output_file2 = "/Users/laraladney/Documents/aa_freqs_per_prot.csv" # path to output file 2

# Protein Sequence Retrieval from NCBI based on search term
Entrez.email = input("Enter email: ")                   # user prompted to enter email (tell NCBI who you are to access sequences)

protTerm = input("Enter NCBI search term: ")            # user prompted to enter protein sequence ID
numSeqs = input("How many protein sequences would you like to extract? ") # user prompted to enter # seqs to retrieve

print("Default date range for protein sequence extraction is from 01/01/2000 - Current.") 
print("Would you like to extract protein sequences from a specified date range?") 
startDate = "2000/01/01"                                # start date for default date range
today = date.today()                                    # current date pulled using datetime module for end date for default date range
endDate = today.strftime("%d/%m/%Y")                    # end date for default date range
dateY_N = ""                                            # initialize user-specified option to opt with default date range or enter their own
while dateY_N != "N" and dateY_N != "n" and dateY_N != "Y" and dateY_N != "y": # while loop to loop through user-entered options until "Y" (yes) or "N" (no) entered
    dateY_N = input("Enter (Y/N): ")                    # user prompted to enter "Y" (yes) or "N" (no) in regard to setting their own date range
    if dateY_N != "N" and dateY_N != "n" and dateY_N != "Y" and dateY_N != "y": # if the user enters a value other than "Y" or "N" they will be asked to enter one of those options
        print("Please choose yes (Y) or no (N).")
if dateY_N == "Y" or "y" and dateY_N != "N" and dateY_N != "n": # if the user enters "Y", then they will be prompted for further inputs
    startDate == input("Using format YYYY/MM/DD, enter start date: ") # user prompted to enter start date
    endDate == input("Using format YYYY/MM/DD, enter end date: ") # user prompted to enter end date
    
searchResultHandle = Entrez.esearch(db = "protein", term = protTerm, retmax = numSeqs, idtype = "protein", datetype = "pdat", mindate = startDate, maxdate = endDate) # entrez search handle with user-specified options
searchResult = Entrez.read(searchResultHandle)          # read in handle with parameters set to user inputs 
ids = searchResult["IdList"]                            # list of IDs created from protein sequences retrieved

handle = Entrez.efetch(db="protein", id=ids, rettype="fasta", retmode="text") # protein sequences retrieved by IDs in fasta format
record = handle.read()                                  # record created reading in the handle containing the fasta format protein sequences           

input_file = open('proteinSearch.txt', 'w')             # input file created in write format to write protein sequences in fasta format
input_file.write(record.rstrip('\n'))                   # each fasta format protein sequence is stripped of the new line character
input_file.close()                                      # close the input file containing the fasta format protein sequences

input_file = open('proteinSearch.txt', 'r')             # reopen input file but in read format

protein_ids = []                                        # make a list of protein IDs

# Define function to extract amino acid sequences from a FASTA file
def extract_amino_acids(input_file):                    # function to extract amino acids from the input file
    amino_acids = []                                    # create a list of amino acids
    for record in SeqIO.parse(input_file, "fasta"):     # parse through each record in the input file using the SeqIO parser to detect fasta formatted sequences
        sequence = str(record.seq)                      # sequence from each record obtained
        protein = ProteinAnalysis(sequence)             # protein analysis set to be performed on each individual sequence
        aa_percentages = protein.get_amino_acids_percent() # percentages of amino acids obtained from protein sequence
        protein_ids.append(record.id)                   # each protein ID appended to the list of protein IDs
        amino_acids.append(aa_percentages)              # each amino acid percentage appended to list of amino acid percentages
    return amino_acids                                  # list of amino acid percentages returned
  
# Define function to write amino acid frequencies and min/max percentages to a CSV file
def write_csv(output_file, data):
    with open(output_file, "w", newline="") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["Rank", "Amino Acid", "Frequency"])
        rank = 1
        for aa, freq in data.items():

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
            
sig = int(input("How many significant figures would you like to preserve for amino acid frequencies? "))
for k, v in total_aa_freqs.items():
    total_aa_freqs[k] = round((v/int(numSeqs)),sig-int(floor(log10(abs(v/int(numSeqs)))))-1)
    
total_aa_freqs_sorted = dict(sorted(total_aa_freqs.items(), key = lambda item: item[1], reverse = True))
          
# Calculate Min/Max %'s of Frequencies
min_percent = min(total_aa_freqs.values())
max_percent = max(total_aa_freqs.values())

# Write amino acid frequencies to output file
write_csv(output_file, total_aa_freqs_sorted)

aa_outfile = open("test.out.csv", "w", newline="") # prints protein frequencies 
writer = csv.writer(aa_outfile, delimiter=",")

protein_outfile = open("protein.csv","w")

# Creates and writes headers
header = ["Protein ID"] + list(amino_acids[0].keys())
protein_outfile.write(",".join(header)+'\n')
writer.writerow(header)

# Writes amino acid frequencies
for i in range(len(protein_ids)):
    aa_dict = amino_acids[i]
    aa_freqlist = list(aa_dict.values())
    row = protein_ids[i] + "," + ",".join([str(x) for x in aa_freqlist])
    protein_outfile.write(row + "\n")
    writer.writerow(row)

aa_outfile.close()
protein_outfile.close()

# Plot amino acid frequencies using seaborn
sns.set_style("whitegrid")
plt.figure(figsize=(10, 6))
sns.barplot(x=list(total_aa_freqs.keys()), y=list(total_aa_freqs.values()), palette="Blues")
plt.title("Amino Acid Frequencies")
plt.xlabel("Amino Acid")
plt.ylabel("Frequency")
plt.text(-0.6, min_percent-0.005, f"Min: {min_percent:.2%}")
plt.text(19.2, max_percent-0.005, f"Max: {max_percent:.2%}")
plt.show()
