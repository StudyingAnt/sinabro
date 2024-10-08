import os
import math
import sys

import numpy as np
import pandas as pd
import multiprocessing as mp
from pathlib import Path, PurePath

from Bio.Seq import Seq, MutableSeq
from typing import Union

sys.path.insert(0, "/mnt/c/Users/CEEL-PC-005/Desktop/Joon/sinabro")
#sys.path.insert(0, "/home/jhsong/sinabro")
from sinabro._helper._helper import (
    MutInfo,
    _get_target_of_mut_type,
    _get_result_of_mut_type,
    _get_motif_indices,
    _get_idx_offset_of_mut_type,
    _get_idx_target_from_idx_motif,
    _get_idx_motif_from_idx_target,
    generate_random_sequence
)

from sinabro._mutate._mutate_seq_with_mutational_signature import (
    make_full_mutational_signature,
    compute_mut_prob_matrix,
    mutate_seq_with_mutational_signature
)

signatures = [
    "SBS1", "SBS2", "SBS3", "SBS4", "SBS5", "SBS6", "SBS7a", "SBS7b", 
    "SBS7c", "SBS7d", "SBS8", "SBS9", "SBS10a", "SBS10b", "SBS10c", 
    "SBS10d", "SBS11", "SBS12", "SBS13", "SBS14", "SBS15", "SBS16", 
    "SBS17a", "SBS17b", "SBS18", "SBS19", "SBS20", "SBS21", "SBS22", 
    "SBS23", "SBS24", "SBS25", "SBS26", "SBS27", "SBS28", "SBS29", 
    "SBS30", "SBS31", "SBS32", "SBS33", "SBS34", "SBS35", "SBS36", 
    "SBS37", "SBS38", "SBS39", "SBS40", "SBS41", "SBS42", "SBS43", 
    "SBS44", "SBS45", "SBS46", "SBS47", "SBS48", "SBS49", "SBS50", 
    "SBS51", "SBS52", "SBS53", "SBS54", "SBS55", "SBS56", "SBS57",
    "SBS58", "SBS59", "SBS60", "SBS84", "SBS85", "SBS86", "SBS87", 
    "SBS88", "SBS89", "SBS90", "SBS91", "SBS92", "SBS93", "SBS94", 
    "SBS95"
]


def simulation(signature):
    n_sims = 10000
    for run in [1]:
        for length in [2200]:
            for gc in [0.4]:
                mut_type_cnts = {}

                for mut in ["[C>A]","[C>G]","[C>T]","[T>A]","[T>C]","[T>G]"]:
                    for nt1 in ["A", "C", "G", "T"]:
                        for nt2 in ["A", "C", "G", "T"]:
                            mut_type_cnts[f"{nt1}{mut}{nt2}"] = 0

                for i in range(n_sims):
                    seq = generate_random_sequence(length, gc=gc)
                    mutinfo = mutate_seq_with_mutational_signature(seq, mutational_signature=signature)
                    mut_type_cnts[mutinfo.mut_type] = mut_type_cnts[mutinfo.mut_type]+1

                rlt = {key: value/n_sims for key, value in mut_type_cnts.items()}

                data = pd.DataFrame(
                    {
                        "mut_type": rlt.keys(),
                        "simulation": rlt.values()
                    }
                )
                output_path = Path(f"/mnt/c/Users/CEEL-PC-005/Desktop/Joon/signature_interaction/exp/20221207-validation/out/run_{run}/{signature.lower()}")
                output_path.mkdir(parents=True, exist_ok=True)
                #data.to_csv(f"/home/jhsong/sinabro_validation/run_{run}/{signature.lower()}/{signature}_length_{str(length)}_gc_{str(gc)}_sim_{run}.csv", index=False)
                data.to_csv(f"/mnt/c/Users/CEEL-PC-005/Desktop/Joon/signature_interaction/exp/20221207-validation/out/run_{run}/{signature.lower()}/{signature}_length_{str(length)}_gc_{str(gc)}_sim_{run}.csv", index=False)

with mp.Pool(processes=4) as pool:
    pool.map(simulation, signatures)

