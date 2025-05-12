import re
from typing import Union
from Bio.Seq import Seq, MutableSeq

from ..types.types import MutInfo, MutationRecord

from Bio.codonalign.codonseq import CodonSeq, cal_dn_ds

import numpy as np 
import warnings

# constant: standard DNA stop codons
STOP_CODONS = {"TAA", "TAG", "TGA"}
STOP_RE = re.compile(r"TAA|TAG|TGA") 

def _has_inframe_stop(seq_without_terminal: str) -> bool:
    """
    Return True if any in-frame stop codon (TAA/TAG/TGA) exists
    **before** the terminal codon.
    """
    m = STOP_RE.search(seq_without_terminal)
    # no stop found → False
    if m is None:
        return False
    # check reading frame (index % 3 == 0)
    return (m.start() % 3) == 0

def compute_dnds(records, **kwargs):
    """
    Compute pairwise dN, dS, and ω (dN/dS) between the first and last
    sequence in a list of MutationRecord object.
    """
    method = kwargs.get("method", "NG86")

    threshold_hit = kwargs.get("threshold_hit", False)
    if threshold_hit:
        last_idx = -1
    else:
        last_idx = -2

    # build trimmed strings
    ref_seq = str(records[0].sequence[1:-4])    # remove first bp + stop
    mut_seq = str(records[last_idx].sequence[1:-4])   # remove first bp + stop

    # check stop loss
    last_codon = str(records[last_idx].sequence[-4:-1])
    stop_loss = last_codon not in STOP_CODONS

    # detect premature in-frame stops in mutant (ignore terminal codon)
    has_premature_stop = _has_inframe_stop(mut_seq)

    # run dN/dS; handle possible stop-codon errors gracefully
    try:
        ref_codon = CodonSeq(ref_seq)
        mut_codon = CodonSeq(mut_seq)
        dN, dS = cal_dn_ds(ref_codon, mut_codon, method=method)
        omega = (dN / dS) if dS else np.inf
    except KeyError:
        # Biopython raises KeyError if a stop codon sneaks in
        warnings.warn("Premature stop codon encountered; "
                      "dN/dS set to NaN.", RuntimeWarning)
        dN = dS = omega = np.nan

    return dN, dS, omega, has_premature_stop, stop_loss