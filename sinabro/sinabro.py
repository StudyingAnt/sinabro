#########1#########2#########3#########4#########5#########6#########7$$#######
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


def is_codon_synonymous(codon1, codon2):
    codon_table = CodonTable.standard_dna_table.forward_table
    codon_table["TAA"] = "-"
    codon_table["TAG"] = "-"
    codon_table["TGA"] = "-"
    
    if codon_table[codon1] == codon_table[codon2]:
        return True
    else:
        return False


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
        return str(seq), description
    elif data_type == "Seq":
        return Seq(seq), description
    elif data_type == "MutableSeq":
        return seq, description
    
    
class Trajectory():    
    def __init__(self, data):
        self._data = []
        if data is None:
            raise ValueError("data must not be None")
        elif isinstance(data, str):
            self._original_seq = Seq(data)
        elif isinstance(data, Seq):
            self._original_seq = data
        elif isinstance(data, MutableSeq):
            self._original_seq = Seq(data)
        else:
            raise TypeError(
                "data should be a string, Seq, or MutableSeq object"
            )
        self._data.append(self._original_seq)
        self._length = len(self._original_seq)
        self._descriptions = ["original"]

    def append(self, seq, description=""):
        if seq is None:
            raise ValueError("seq must not be None")
        elif isinstance(seq, str):
            self._data.append(Seq(seq))
        elif isinstance(seq, Seq):
            self._data.append(seq)
        elif isinstance(seq, MutableSeq):
            self._data.append(Seq(seq))
        else:
            raise TypeError(
                "data should be a string, Seq, or MutableSeq object"
            )
        self._descriptions.append(description)

        return self._data
        
    def add_mutated_sequence(self, method="random", motif=None, P=None):
        if method == "random":
            last_seq = self._data[-1]
            new_seq, description = random_single_substitution(last_seq,
                                                              start=1,
                                                              end=self._length-2
                                                              )
            self.append(new_seq, description=description)
        elif method == "motif":
            pass
        elif method == "signature":
            pass
        else:
            raise ValueError(
                "mode must be either random, motif, or signature"
            )

        return self._data
                
