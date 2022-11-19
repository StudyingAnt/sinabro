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

from typing import Union

from _mutate._random_single_substitution import (
    random_single_substitution
)
from _mutate._mutate_seq_with_mut_type import (
    mutate_seq_with_mut_type
)
from _mutate._mutate_seq_with_mutational_signature import (
    mutate_seq_with_mutational_signature
)

from _mutate._helper import MutInfo


def is_codon_synonymous(codon1, codon2):
    codon_table = CodonTable.standard_dna_table.forward_table
    codon_table["TAA"] = "-"
    codon_table["TAG"] = "-"
    codon_table["TGA"] = "-"
    
    if codon_table[codon1] == codon_table[codon2]:
        return True
    else:
        return False

    
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
        self._seq_length = len(self._original_seq)
        self._hgvss = ["."]
        self._mut_types = ["."]
        self._length = 0
    
    def show(self):
        print(f"Sequence\tHGVS\tMutation Type")
        for i in range(self._length):
            line_tokens = [
                str(self._data[i]),
                self._hgvss[i],
                self._mut_types[i] 
            ]
            line = "\t".join(line_tokens)
            print(line)

    def append(self, seq, hgvs="", mut_type=""):
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
        self._hgvss.append(hgvs)
        self._mut_types.append(mut_type)

        return self._data
        
    def add_mutated_sequence(
        self,
        method: str = "random",
        mut_type: str = None,
        mutational_signature: str = None,
        cosmic_version: float = 3.3,
        genome_ref: str = "GRCh37",
        custom_signature_path: str = None,
        column: str = None,
        strand_bias: float = 0.5
    ) -> list:
        if method == "random":
            last_seq = self._data[-1]
            mutinfo = random_single_substitution(
                last_seq,
                start=1,
                end=self._seq_length-2
                )
            if not mutinfo.e:
                self.append(mutinfo.new_seq, 
                            hgvs=mutinfo.hgvs, 
                            mut_type=mutinfo.mut_type)
                self._length = self._length+1
        elif method == "mut_type":
            last_seq = self._data[-1]
            mutinfo = mutate_seq_with_mut_type(
                last_seq,
                mut_type=mut_type,
                start=1,
                end=self._seq_length-2
                )
            if not mutinfo.e:
                self.append(mutinfo.new_seq, 
                            hgvs=mutinfo.hgvs, 
                            mut_type=mutinfo.mut_type)
                self._length = self._length+1
            else:
                pass
        elif method == "signature":
            last_seq = self._data[-1]
            mutinfo = mutate_seq_with_mutational_signature(
                last_seq,
                mutational_signature=mutational_signature,
                cosmic_version=cosmic_version,
                genome_ref=genome_ref,
                custom_signature_path=custom_signature_path,
                column=column,
                strand_bias=strand_bias
            )
            if not mutinfo.e:
                self.append(mutinfo.new_seq, 
                            hgvs=mutinfo.hgvs, 
                            mut_type=mutinfo.mut_type)
                self._length = self._length+1
        else:
            raise ValueError(
                "mode must be either random, motif, or signature"
            )

        return self._data

    def autofill(self, codition="max_length", max_length=10, method=None, *kwargs):
        