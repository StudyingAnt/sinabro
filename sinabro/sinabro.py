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

from ._mutate._random_single_substitution import (
    random_single_substitution
)
from ._mutate._mutate_seq_with_mut_type import (
    mutate_seq_with_mut_type
)
from ._mutate._mutate_seq_with_mutational_signature import (
    mutate_seq_with_mutational_signature
)

from ._helper._helper import (
    MutInfo,
    _is_codon_synonymous,
    _get_codon,
    _get_pos_aa,
    _get_idx_from_hgvs_mrna,
    _get_amino_acid_from_codon
)

    
class Trajectory():    
    def __init__(self, id, data, note="."):
        self._id = id
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
        self._hgvs_mrnas = ["."]
        self._hgvs_aa = ["."]
        self._mut_types = ["."]
        self._notes = [note]        
        self._length = 0
    
    def show(self, verbose= True):
        #print(f"Sequence\tHGVS_mRNA\tHGVS_Protein\tMutation_Type")
        print(f"ID: {self._id}")
        print(f"Trajectory length: {self._length}")

        if verbose:
            print(f"")
            print(f"Sequence\tHGVS_mRNA\tHGVS_Protein\tMutation_Type")
            #print(f"Sequence\tHGVS_mRNA\tMutation_Type")
            for i in range(self._length+1):
                line_tokens = [
                    str(self._data[i]),
                    self._hgvs_mrnas[i],
                    self._hgvs_aa[i],
                    self._mut_types[i] 
                ]
                line = "\t".join(line_tokens)
                print(line)
        else:
            pass

    def append(self, seq, hgvs_mrna=".", mut_type=".", note="."):
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
        self._hgvs_mrnas.append(hgvs_mrna)
        self._mut_types.append(mut_type)
        self._notes.append(note)

        return self._data

    def pop(self):
        last_seq = self._data.pop()
        last_hgvs_mrna = self._hgvs_mrnas.pop()
        last_hgvs_aa = self._hgvs_aa.pop()
        last_mut_type = self._mut_types.pop()
        self._length = self._length-1

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
                            hgvs_mrna=mutinfo.hgvs_mrna, 
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
                            hgvs_mrna=mutinfo.hgvs_mrna, 
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
        if condition == "max_length": # have to write
            pass
        elif condition == "nonsynonymous": 
            if method == "random": # have to write
                last_seq = self._data[-1]
                self.add_mutated_sequence(last_seq)
            elif method == "mut_type": #have to write
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

                    idx_target = _get_idx_from_hgvs_mrna(self._hgvs_mrnas[-1])

                    new_codon = _get_codon(curr_seq, idx_target)
                    old_codon = _get_codon(prev_seq, idx_target)
                    
                    if _is_codon_synonymous(new_codon, old_codon):
                        self._hgvs_aa.append(".")
                    else:
                        pos_aa = _get_pos_aa(idx_target)

                        new_aa = _get_amino_acid_from_codon(new_codon, three=True)
                        old_aa = _get_amino_acid_from_codon(old_codon, three=True)

                        self._hgvs_aa.append(f"p.{pos_aa}{old_aa}>{new_aa}")
                        stop_flag = True
            else:
                raise ValueError(
                    "condition must be either max_length or nonsynonymous"
                )
        else:
            raise ValueError(
                "condition must be either max_length or nonsynonymous"
            )

        return self._data

    def save(self, output_dir_path, prefix=None, suffix=None):
        records = []
        for i in range(self._length+1):
            description = f"{self._hgvs_mrnas[i]};{self._hgvs_aa[i]};{self._mut_types[i]}"
            if isinstance(self._data[i], str):
                seq = Seq(self._data[i])
            else:
                seq = self._data[i]
            record = SeqRecord(seq, id=self._id, description=description)
            records.append(record)

        if prefix:
            prefix = f"{prefix}_"
        else:
            prefix = ""

        if suffix:
            suffix = f"_{suffix}"
        else:
            suffix = ""

        file_name = f"{prefix}trajectory{suffix}.fa"
        output_path = os.path.join(output_dir_path, file_name)

        SeqIO.write(records, output_path, "fasta")
