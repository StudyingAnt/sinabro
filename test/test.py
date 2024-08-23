import sys
#sys.path.insert(0, "/mnt/c/Users/CEEL-PC-005/Desktop/Joon/sinabro")
#sys.path.insert(0, "/home/augustine/Desktop/Labs/0_Dry/projSETA/sinabro")
import sinabro.sinabro as snbr
from sinabro._helper._helper import _reverse_complement_mut_type
from sinabro._helper._helper import _unpack_complex_mut_type
from sinabro._helper._helper import generate_random_sequence
from sinabro._mutate._mutate_seq_with_mut_type import mutate_seq_with_mut_type

import os
from pathlib import Path, PurePath

import numpy as np

import pandas as pd

from Bio import SeqIO

#new_mut_type = _reverse_complement_mut_type("T[C>T]C")
#print(new_mut_type)

#mut_types = _unpack_complex_mut_type("TY[C>T]")
#print(mut_types)

# 35 bp 11 codon 10 aa
# cnt  12345678911234567892123456789312345
# idx  01234567891123456789212345678931234
# c.  -1123456789112345678921234567893123+1  
#seq = "TATGCTGACTCGGTCATCGATCGGTTCTCATTGAT"
seq = "TATG"+str(generate_random_sequence(3*8))+"TGAT"
#print(seq)
#seq = "TATGAGTTGGGACCCATCGTAAGTAACCTGAT"

#print(np.array([7,4,3])*np.array([1/4,1/3,5/12]))

#print(2*[1/4,1/3,5/12])

traj = snbr.Trajectory(id="Test", data=seq)
traj.autofill(condition="max_length", method="mut_types", strand="both", start = 1, mut_types=["AA[C>T]", "AC[C>T]", "TT[C>T]"], mut_type_probs = [1/4,1/3,5/12])
traj.show()

##
bat = "molMol2"

exp_dir = PurePath("/mnt/c/Users/CEEL-PC-005/Desktop/Joon/Final_scripts/Simple_robustness_analysis_of_cds_files/Bat_genome")
in_dir = PurePath(exp_dir, "Bat1KPilotProject/HLmolMol2")
out_dir = PurePath(exp_dir, "results")

genes_fasta_file = PurePath(in_dir, "HL"+bat+".all.pc_transcripts.cdsplus.fa")

seq_records = list(SeqIO.parse(genes_fasta_file, "fasta"))

gene_idx = 1000
transcript_name = seq_records[gene_idx].id.split(".")[0]
seq = seq_records[gene_idx].seq

apobec_prob = pd.read_csv(PurePath(in_dir, bat+"_apobec.csv"))

mut_types = apobec_prob["mut_type"]
apobecs = ["APOBEC3A", "APOBEC3C", "APOBEC1"]

rlt_dict = {}
for apobec in apobecs:
    rlt_dict[apobec] = []

n_traj = 1000
for apobec in apobecs:
    for _ in range(n_traj):
        traj = snbr.Trajectory(transcript_name, seq)
        traj.autofill(condition="nonsynonymous", method="mut_types", strand="both",
                      start = 1, end = traj._seq_length-1,
                      mut_types=mut_types, mut_type_probs=apobec_prob[apobec])
        #traj.show()
        print(" ")

        rlt_dict[apobec].append(traj._length-1)

out_data = pd.DataFrame(rlt_dict)
out_file = PurePath(out_dir, f"{transcript_name}_rus_apobec.csv")
out_data.to_csv(out_file, index=False)







"""
# 
seq = "TATGCTGACTCGGTCATCGATCGGTTCTCATTGAT"
traj = snbr.Trajectory(id="Test", data=seq)
traj.autofill(condition="max_length", method="random")
traj.show(show_idx=True)
"""
# 
"""
seq = "TATGCTGACTCGGTCATCGATCGGTTCTCATTGAT"
traj = snbr.Trajectory(id="Test", data=seq)
traj.autofill(condition="evaluation", method="signature", mutational_signature="SBS2", measure="blastp", threshold=10)
traj.show(show_idx=True)
"""
"""
# 
seq = "TATGCTGACTCGGTCATCGATCGGTTCTCATTGAT"
traj = snbr.Trajectory(id="Test", data=seq)
traj.autofill(condition="max_length", method="signature", mutational_signature="SBS2")
traj.show(show_idx=True)
"""
"""
# 
seq = "TATGCTGACTCGGTCATCGATCGGTTCTCATTGAT"
traj = snbr.Trajectory(id="Test", data=seq)
traj.autofill(condition="nonsynonymous", method="random")
traj.show()

# 
seq = "TATGCTGACTCGGTCATCGATCGGTTCTCATTGAT"
traj = snbr.Trajectory(id="Test", data=seq)
traj.autofill(condition="nonsynonymous", method="mut_type", mut_type="T[C>T]")
traj.show()

# 
seq = "TATGCTGACTCGGTCATCGATCGGTTCTCATTGAT"
traj = snbr.Trajectory(id="Test", data=seq)
traj.autofill(condition="nonsynonymous", method="signature", mutational_signature="SBS2")
traj.show()

# 
#seq = "TATGCTGACTCGGTCATCGATCGGTTCTCATTGAT"
#traj = snbr.Trajectory(id="Test", data=seq)
#traj.autofill(condition="max_length", method="mut_type", mut_type="T[G>A]A")
#traj.show(show_idx=True)





traj1 = snbr.Trajectory("test", seq)

traj1.show()
#traj1.autofill(condition="max_length", method="random")
traj1.autofill(condition="max_length", method="mut_type", mut_type="T[C>T]")
#traj1.autofill(condition="max_length", method="signature", mutational_signature="SBS2")
#traj1.autofill(condition="nonsynonymous", method="random")
#traj1.autofill(condition="nonsynonymous", method="mut_type", mut_type="T[C>T]")
#traj1.autofill(condition="nonsynonymous", method="signature", mutational_signature="SBS2")
#traj1.pop()
#traj1.autofill(condition="nonsynonymous", method="signature", mutational_signature="SBS1")
#traj1.add_mutated_sequence(method="signature", mutational_signature="SBS2")
#traj1.add_mutated_sequence()
#traj1.add_mutated_sequence(method="mut_type", mut_type="T[C>T]")

traj1.show()

#print(traj1._hgvs_aa)
#home = "/mnt/c/Users/CEEL-PC-005"
#output_path = PurePath(home, "Desktop/Joon")
#traj1.save(output_path, prefix="Test", suffix=str(1))
"""