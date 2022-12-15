import sys
#sys.path.insert(0, "/mnt/c/Users/CEEL-PC-005/Desktop/Joon/sinabro")
#sys.path.insert(0, "/home/augustine/Desktop/Labs/0_Dry/projSETA/sinabro")
import sinabro.sinabro as snbr

import os
from pathlib import Path, PurePath


# 35 bp 11 codon 10 aa
# cnt  12345678911234567892123456789312345
# idx  01234567891123456789212345678931234
# c.  -1123456789112345678921234567893123+1  
seq = "TATGCTGACTCGGTCATCGATCGGTTCTCATTGAT"

"""
# 
seq = "TATGCTGACTCGGTCATCGATCGGTTCTCATTGAT"
traj = snbr.Trajectory(id="Test", data=seq)
traj.autofill(condition="max_length", method="random")
traj.show(show_idx=True)
"""
# 
seq = "TATGCTGACTCGGTCATCGATCGGTTCTCATTGAT"
traj = snbr.Trajectory(id="Test", data=seq)
traj.autofill(condition="max_length", method="mut_type", mut_type="T[G>A]A")
traj.show(show_idx=True)

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