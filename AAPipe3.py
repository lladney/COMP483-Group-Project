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

# Define input and output files
input_file = open('proteinSearch.txt', 'w')  # open and write sequences to proteins text file 
output_file = "/Users/param/Downloads/sars-cov2_amino_acid_frequencies-prot1.csv"
output_file2 = "/Users/param/Downloads/aa_freqs_per_prot.csv"

# Protein Sequence Retrieval from NCBI based on search term
Entrez.email = "wmccain@luc.edu" # user prompted to enter email (tell NCBI who you are to access sequences)

protTerm = "sars cov 2" # user prompted to enter protein sequence ID
numSeqs = 10 # user prompted to enter # seqs

dateY_N = "N"
startDate = ""
endDate = ""
if dateY_N == "N" or "n":
    startDate == "2000/01/01" and endDate == today.strftime("%d/%m/%Y")

searchResultHandle = Entrez.esearch(db = "protein", term = protTerm, retmax = numSeqs, mindate = startDate, maxdate = endDate)
searchResult = Entrez.read(searchResultHandle)
ids = searchResult["IdList"]

handle = Entrez.efetch(db="protein", id=ids, rettype="fasta", retmode="text")
record = handle.read()

input_file = open('proteinSearch.txt', 'w')
input_file.write(record.rstrip('\n'))
input_file.close()


input_file = open('proteinSearch.txt', 'r')

# Define function to extract amino acid sequences and protein IDs from a FASTA file
def extract_amino_acids(input_file):
    records = SeqIO.parse(input_file, "fasta")
    amino_acids = []
    for record in records:
        protein_id = record.description.split("|")
        if len(protein_id) > 1:
            protein_id = protein_id[1]
        else:
            protein_id = record.id
        sequence = str(record.seq)
        protein = ProteinAnalysis(sequence)
        aa_percentages = protein.get_amino_acids_percent()
        amino_acids.append((protein_id, sequence, aa_percentages))
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
    for dict_item in aa:
        if isinstance(dict_item, dict):
            for k, v in dict_item.items():
                if k in total_aa_freqs:
                    total_aa_freqs[k] += v
                else:
                    total_aa_freqs[k] = v

for k, v in total_aa_freqs.items():
    total_aa_freqs[k] = v/int(numSeqs)

total_aa_freqs_sorted = dict(sorted(total_aa_freqs.items(), key = lambda item: item[1], reverse = True))

# Calculate Min/Max %'s of Frequencies
if len(total_aa_freqs) > 0:
    min_percent = min(total_aa_freqs.values())
    max_percent = max(total_aa_freqs.values())
else:
    min_percent = 0
    max_percent = 0

# Write amino acid frequencies to output file
write_csv(output_file, total_aa_freqs_sorted)

# Define function to write amino acid frequencies and protein IDs to a CSV file
def write_csv2(output_file2, data):
    with open(output_file2, "w", newline="") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["Protein ID", "Amino Acid", "Frequency"])
        #for protein_id, aa in data:
        for protein_id, aa, freq in data.items():
            writer.writerow([protein_id, aa, freq])

aa_outfile = open("test.out2.csv", "w", newline="")
writer = csv.writer(aa_outfile, delimiter=",")

header = ["Protein ID"]
if isinstance(amino_acids[0], dict):
    header += list(amino_acids[0].keys())
    writer.writerow(header)

for aa_dict in amino_acids:
    row = [aa_dict[0]]
    if isinstance(aa_dict[1], dict):
        row += list(aa_dict[1].values())
    writer.writerow(row)

aa_outfile.close()
