import re
from typing import Union
from Bio.Seq import Seq, MutableSeq

from ..types.types import MutInfo, MutationRecord
from ..evaluate.algn_scoring import compute_align_distance

def dna_to_aa(seq: str, table: int = 1, to_stop: bool = False) -> str:
    """
    Translate a DNA string to an amino-acid string.
    - table: NCBI translation table ID (1 = standard)
    - to_stop=False keeps '*' for internal/terminal stops
    """
    # Create Seq object and translate; stop_symbol defaults to '*'
    aa = Seq(seq).translate(table=table, to_stop=to_stop, stop_symbol='*')
    return str(aa)


def hamming(a: str, b: str) -> int:
    """
    Return Hamming distance between equal-length strings a and b.
    """
    a = dna_to_aa(a)
    b = dna_to_aa(b)
    if len(a) != len(b):
        raise ValueError("Sequences must be the same length for Hamming distance")
    return sum(ch1 != ch2 for ch1, ch2 in zip(a, b))


def compute_distance(records, **kwargs):
    methods = {
        "hamming": hamming,
        "blosum": compute_align_distance
    }

    method = kwargs.get("method", "hamming")
    if method not in methods:
        raise ValueError(
            f"Unsupported method '{method}'. "
            f"Choose one of: {', '.join(methods.keys())}"
        )

    threshold_hit = kwargs.get("threshold_hit", False)
    if threshold_hit:
        last_idx = -1
    else:
        last_idx = -2

    ref_seq = str(records[0].sequence[1:-1])    
    mut_seq = str(records[last_idx].sequence[1:-1])   
        
    distance = methods[method](ref_seq, mut_seq)

    return distance     



