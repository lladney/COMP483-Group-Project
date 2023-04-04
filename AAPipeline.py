import sys
import csv
from Bio import SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import matplotlib.pyplot as plt
import seaborn as sns

# Define input and output files
input_file = "C:/Users/param/Downloads/sars-cov-2-prot1.txt"
output_file = "C:/Users/param/Downloads/sars-cov2_amino_acid_frequencies-prot1.csv"

# Define function to extract amino acid sequences from a FASTA file
def extract_amino_acids(input_file):
    amino_acids = []
    for record in SeqIO.parse(input_file, "fasta"):
        sequence = str(record.seq)
        protein = ProteinAnalysis(sequence)
        aa_percentages = protein.get_amino_acids_percent()
        aa_freqs = {}
        for aa, percentage in aa_percentages.items():
            aa_freqs[aa] = percentage / 100
        amino_acids.append(aa_freqs)
    return amino_acids

# Define function to write amino acid frequencies and min/max percentages to a CSV file
def write_csv(output_file, data):
    with open(output_file, "w", newline="") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["Amino Acid", "Frequency", "Min %", "Max %"])
        for aa, freq in data.items():
            min_freq = min(amino_acids, key=lambda x: x[aa])[aa]
            max_freq = max(amino_acids, key=lambda x: x[aa])[aa]
            writer.writerow([aa, freq, min_freq, max_freq])

# Extract amino acid sequences from input file
amino_acids = extract_amino_acids(input_file)

# Calculate total frequencies of each amino acid
total_aa_freqs = {}
for aa in amino_acids:
    for k, v in aa.items():
        if k in total_aa_freqs:
            total_aa_freqs[k] += v
        else:
            total_aa_freqs[k] = v
            
# Calculate Min/Max %'s of Frequencies
min_percent = min(total_aa_freqs.values())
max_percent = max(total_aa_freqs.values())

# Write amino acid frequencies to output file
write_csv(output_file, total_aa_freqs)

# Plot amino acid frequencies using seaborn
sns.set_style("whitegrid")
plt.figure(figsize=(10, 6))
sns.barplot(x=list(total_aa_freqs.keys()), y=list(total_aa_freqs.values()), palette="Blues")
plt.title("Amino Acid Frequencies in SARS-CoV-2 Proteins")
plt.xlabel("Amino Acid")
plt.ylabel("Frequency")
plt.text(-0.6, min_percent-0.005, f"Min: {min_percent:.2%}")
plt.text(19.2, max_percent-0.005, f"Max: {max_percent:.2%}")
plt.show()
