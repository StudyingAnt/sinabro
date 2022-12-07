import time
import numpy as np

from Bio.Seq import Seq, MutableSeq
from typing import Union

from .._helper._helper import (
    MutInfo,
    _get_target_of_mut_type,
    _get_result_of_mut_type,
    _get_motif_indices,
    _get_idx_offset_of_mut_type,
    _get_idx_target_from_idx_motif,
    _get_idx_motif_from_idx_target
)


def mutate_seq_with_mut_type(seq, mut_type, start=None, end=None):
    if isinstance(seq, str):
        seq = MutableSeq(seq)
        data_type = "str"
    elif isinstance(seq, Seq):
        seq = MutableSeq(seq)
        data_type = "Seq"
    elif isinstance(seq, MutableSeq):
        data_type = "MutableSeq"
    else:
        raise TypeError(
            "seq should be a string, Seq, or MutableSeq object"
        )

    if start is not None:
        if isinstance(start, int):
            pass                
        else:
            raise TypeError(
                "start should be an integer object"
            )
    else:
        start = 0
        
    if end is not None:
        if isinstance(end, int):
            pass                
        else:
            raise TypeError(
                "end should be an integer object"
            )
    else:
        end = len(seq)-1

    idx_motifs = _get_motif_indices(seq, mut_type)
    idx_motifs = [idx for idx in idx_motifs if idx >= start and idx <= end]

    if not idx_motifs:
        return MutInfo(MutableSeq(""), -1, "", "", 1)
    
    # select specific motif to mutate
    seed1 = time.time()
    seed2 = (seed1 - int(seed1))*(10**7)
    np.random.seed(int(seed1+seed2))
    idx_motif = np.random.choice(idx_motifs, 1)[0] 
    idx_target = _get_idx_target_from_idx_motif(idx_motif, mut_type)

    nt_before = seq[idx_target:idx_target+1]
    seq[idx_motif:idx_target+1] = _get_result_of_mut_type(mut_type)
    nt_after = seq[idx_target:idx_target+1]
    
    hgvs_mrna = f"c.{str(idx_target)}{nt_before}>{nt_after}"

    if data_type == "str":
        return MutInfo(str(seq), idx_target, hgvs_mrna, mut_type, 0)
    elif data_type == "Seq":
        return MutInfo(Seq(seq), idx_target, hgvs_mrna, mut_type, 0)
    elif data_type == "MutableSeq":
        return MutInfo(seq, idx_target, hgvs_mrna, mut_type, 0)