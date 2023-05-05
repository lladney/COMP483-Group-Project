# Import necessary packages
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
import pandas as pd
from datetime import date
today = date.today()
from math import log10, floor

# Define input and output files
input_file = open('proteinSearch.txt', 'w')  # open and write sequences to proteins text file 
output_file = "/Users/laraladney/Documents/sars-cov2_amino_acid_frequencies-prot1.csv"
output_file2 = "/Users/laraladney/Documents/aa_freqs_per_prot.csv"

# Protein Sequence Retrieval from NCBI based on search term
Entrez.email = input("Enter email: ")   # user prompted to enter email (tell NCBI who you are to access sequences)

protTerm = input("Enter NCBI search term: ")     # user prompted to enter protein sequence ID
numSeqs = input("How many protein sequences would you like to extract? ") # user prompted to enter # seqs

print("Default date range for protein sequence extraction is from 01/01/2000 - Current.")
print("Would you like to extract protein sequences from a specified date range?")
startDate = "2000/01/01"
endDate = today.strftime("%d/%m/%Y")
dateY_N = ""
while dateY_N != "N" and dateY_N != "n" and dateY_N != "Y" and dateY_N != "y":
    dateY_N = input("Enter (Y/N): ")
    if dateY_N != "N" and dateY_N != "n" and dateY_N != "Y" and dateY_N != "y":
        print("Please choose yes (Y) or no (N).")
if dateY_N == "Y" or "y" and dateY_N != "N" and dateY_N != "n":
    startDate == input("Using format YYYY/MM/DD, enter start date: ")
    endDate == input("Using format YYYY/MM/DD, enter end date: ")
    
searchResultHandle = Entrez.esearch(db = "protein", term = protTerm, retmax = numSeqs, idtype = "protein", datetype = "pdat", mindate = startDate, maxdate = endDate)
searchResult = Entrez.read(searchResultHandle)
ids = searchResult["IdList"]

handle = Entrez.efetch(db="protein", id=ids, rettype="fasta", retmode="text")
record = handle.read()

input_file = open('proteinSearch.txt', 'w')
input_file.write(record.rstrip('\n'))
input_file.close()

input_file = open('proteinSearch.txt', 'r')

protein_ids = []

# Define function to extract amino acid sequences from a FASTA file
def extract_amino_acids(input_file):
    amino_acids = []
    for record in SeqIO.parse(input_file, "fasta"):
        sequence = str(record.seq)
        protein = ProteinAnalysis(sequence)
        aa_percentages = protein.get_amino_acids_percent()
        protein_ids.append(record.id)
        amino_acids.append(aa_percentages)
        
    return amino_acids

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

header = ["Protein ID"] + list(amino_acids[0].keys())
protein_outfile.write(",".join(header)+'\n')
writer.writerow(header)

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
