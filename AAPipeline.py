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
