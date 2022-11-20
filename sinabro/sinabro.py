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

from _autofill._helper import (
    _is_codon_synonymous,
    _get_codon,
    _get_idx_from_hgvs
)

    
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
        for i in range(self._length+1):
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
        
    def add_mutated_sequence(self, **kwargs) -> list:
        method = kwargs["method"]
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
            else:
                pass
        elif method == "mut_type":
            mut_type = kwargs["mut_type"]
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
            mutational_signature = kwargs["mutational_signature"]
            if "cosmic_version" in kwargs.keys():
                cosmic_version = kwargs["cosmic_version"]
            else:
                cosmic_version = 3.3
            if "genome_ref" in kwargs.keys():
                genome_ref = kwargs["genome_ref"]
            else:
                genome_ref = "GRCh37"
            if "custom_signature_path" in kwargs.keys():
                custom_signature_path = kwargs["custom_signature_path"]
            else:
                custom_signature_path =None
            if "column" in kwargs.keys():
                column = kwargs["column"]
            else:
                column = None
            if "strand_bias" in kwargs.keys():
                strand_bias = kwargs["strand_bias"]
            else:
                strand_bias = 0.5

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
                pass
        else:
            raise ValueError(
                "method must be either random, motif, or signature"
            )

        return self._data

    def autofill(
        self, 
        condition: str = "max_length", 
        max_length: int = 10,  
        **kwargs
        ):
        method = kwargs['method']
        if condition == "max_length":
            pass
        elif condition == "nonsynonymous":
            if method == "random":
                last_seq = self._data[-1]
                self.add_mutated_sequence(last_seq)
            elif method == "mut_type":
                last_seq = self._data[-1]
                self.add_mutated_sequence(last_seq,
                    method="mut_type",
                    mut_type=mut_type
                    )
            elif method == "signature":
                mutational_signature = kwargs["mutational_signature"]
                
                stop_flag = False
                while not stop_flag:
                    last_seq = self._data[-1]
                    self.add_mutated_sequence(**kwargs)

                    curr_seq = self._data[-1]
                    prev_seq = self._data[-2]

                    idx_target = _get_idx_from_hgvs(self._hgvss[-1])

                    new_codon = _get_codon(curr_seq, idx_target)
                    old_codon = _get_codon(prev_seq, idx_target)

                    if new_codon == old_codon:
                        pass
                    else:
                        stop_flag = True
        else:
            raise ValueError(
                "condition must be either max_length or nonsynonymous"
            )

        return self._data