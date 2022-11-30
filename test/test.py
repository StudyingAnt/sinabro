import sys

sys.path.insert(0, "/mnt/c/Users/CEEL-PC-005/Desktop/Joon/sinabro")
import sinabro.sinabro as snbr

import os
from pathlib import Path, PurePath

seq = "TATGCTGACTCGGATCGGTTCTTGAT"

traj1 = snbr.Trajectory("test", seq)

traj1.autofill(condition="max_length", method="random")
#traj1.autofill(condition="nonsynonymous", method="random")
#traj1.autofill(condition="nonsynonymous", method="mut_type", mut_type="T[C>T]")
#raj1.autofill(condition="nonsynonymous", method="signature", mutational_signature="SBS2")
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
