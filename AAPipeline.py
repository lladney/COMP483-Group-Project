"""from Bio.SeqUtils.ProtParam import ProteinAnalysis
X = ProteinAnalysis("MAEGEITTFTALTEKFNLPPGNYKKPKLLYCSNGGHFLRILPDGTVDGT"
                    "RDRSDQHIQLQLSAESVGEVYIKSTETGQYLAMDTSGLLYGSQTPSEEC"
                    "LFLERLEENHYNTYTSKKHAEKNWFVGLKKNGSCKRGPRTHYGQKAILF"
                    "LPLPV")
print(X.count_amino_acids()['A'])
6
print(X.count_amino_acids()['E'])
12
print("%0.2f" % X.get_amino_acids_percent()['A'])
0.04
print("%0.2f" % X.get_amino_acids_percent()['L'])
0.12
print("%0.2f" % X.molecular_weight())
17103.16
print("%0.2f" % X.aromaticity())
0.10
print("%0.2f" % X.instability_index())
41.98
print("%0.2f" % X.isoelectric_point())
7.72
sec_struc = X.secondary_structure_fraction()  # [helix, turn, sheet]
print("%0.2f" % sec_struc[0])  # helix
0.28
epsilon_prot = X.molar_extinction_coefficient()  # [reduced, oxidized]
print(epsilon_prot[0])  # with reduced cysteines
17420
print(epsilon_prot[1])  # with disulfid bridges
17545
"""
import sys
import csv
from Bio import SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import matplotlib.pyplot as plt
import seaborn as sns

# Define input and output files
input_file = "C:/Users/param/Downloads/sars-cov-2-om1.txt"
output_file = "C:/Users/param/Downloads/sars-cov2_amino_acid_frequencies-om1.csv"

# Define function to extract amino acid sequences from a FASTA file
def extract_amino_acids(input_file):
    amino_acids = []
    for record in SeqIO.parse(input_file, "fasta"):
        sequence = str(record.seq)
        protein = ProteinAnalysis(sequence)
        amino_acids.append(protein.get_amino_acids_percent())
    return amino_acids

# Define function to write amino acid frequencies to a CSV file
def write_csv(output_file, data):
    with open(output_file, "w", newline="") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["Amino Acid", "Frequency"])
        for aa, freq in data.items():
            writer.writerow([aa, freq])

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

# Write amino acid frequencies to output file
write_csv(output_file, total_aa_freqs)

# Plot amino acid frequencies using seaborn
sns.set_style("whitegrid")
plt.figure(figsize=(10, 6))
sns.barplot(x=list(total_aa_freqs.keys()), y=list(total_aa_freqs.values()), palette="Blues")
plt.title("Amino Acid Frequencies in SARS-CoV-2 Proteins")
plt.xlabel("Amino Acid")
plt.ylabel("Frequency")
plt.show()
