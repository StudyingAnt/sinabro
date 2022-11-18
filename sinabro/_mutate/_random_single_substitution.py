import re
import os
import math
import random
import sys

import numpy as np
import pandas as pd

from Bio import SeqIO
from Bio.Seq import Seq, MutableSeq
from Bio.SeqRecord import SeqRecord
from Bio.Data import CodonTable

from typing import Union

from ._helper import (
    _get_target_of_mut_type,
    _get_result_of_mut_type,
    _get_motif_indices,
    _get_idx_offset_of_mut_type,
    _get_idx_target_from_idx_motif,
    _get_idx_motif_from_idx_target
)


def random_single_substitution(seq, start=None, end=None):
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

    idx = np.random.choice(range(start, end+1), 1)[0]    
    nt_before = seq[idx:idx+1]
    A = ["A", "C", "G", "T"]
    A.remove(nt_before)
    nt_after = np.random.choice(A, 1,
                                p=np.array([1/3, 1/3, 1/3])
                                )[0]
    seq[idx:idx+1] = nt_after
    description = f"c.{str(idx)}{nt_before}>{nt_after}"

    if data_type == "str":
        return str(seq), idx, description
    elif data_type == "Seq":
        return Seq(seq), idx, description
    elif data_type == "MutableSeq":
        return seq, idx, description