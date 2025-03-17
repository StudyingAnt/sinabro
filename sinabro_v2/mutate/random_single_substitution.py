import numpy as np

from typing import Union
from Bio.Seq import Seq, MutableSeq

from ..types.types import MutInfo, MutationRecord
from ._helper import preserve_seq_type

@preserve_seq_type
def random_single_substitution(seq, **kwargs):
    """
    Performs a random single nucleotide substitution on the given sequence.
    
    The function processes the sequence as a MutableSeq and returns a MutInfo object
    containing the modified sequence and mutation details.
    """
    start = 0
    end = len(seq) - 1

    # Select a random index between start and end (inclusive)
    idx_target = np.random.choice(range(start, end + 1), 1)[0]
    # Get the nucleotide at the target index
    nt_before = seq[idx_target:idx_target + 1]
    # Create a list of nucleotides and remove the nucleotide at the target index
    A = ["A", "C", "G", "T"]
    A.remove(nt_before)
        
    # Using the global random generator without reseeding for improved randomness
    nt_after = np.random.choice(A, 1, p=[1/3, 1/3, 1/3])[0]
    
    # Replace the nucleotide at the target index with the new nucleotide
    seq[idx_target:idx_target + 1] = nt_after
    hgvs_mrna = f"c.{idx_target}{nt_before}>{nt_after}"
    mut_type = f"[{nt_before}>{nt_after}]"
    
    # Always work with MutableSeq type inside the function.
    return MutInfo(seq, idx_target, hgvs_mrna, mut_type, 0)