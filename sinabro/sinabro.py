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


# A list of all trinucleotides. datatype: BioPtyhon Seq class
dna_bases = ["A", "C", "G", "T"]
trinucleotides = []
for nt1 in dna_bases:
    for nt2 in dna_bases:
        for nt3 in dna_bases:
            trinucleotides.append(Seq(f"{nt1}{nt2}{nt3}"))


def make_full_mutational_signature(
    input_file_path: str,
    column: str = None,
    strand_bias: float = 0.5
) -> np.ndarray:
    # read SBS file
    sbs = pd.read_csv(input_file_path, sep="\t", index_col=0, header=0)

    # get column
    if column is None:
        column = sbs.columns[0]

    # fill in full_mutational_signature numpy 2d array
    full_mutational_signature = np.zeros((4,64))
    mut_types = list(sbs.index)

    for mut_type in mut_types:
        trinucleotide_f = f"{mut_type[0]}{mut_type[2]}{mut_type[6]}"
        trinucleotide_r = str(Seq(trinucleotide_f).reverse_complement())
        changed_nt_f = mut_type[4]
        changed_nt_r = str(Seq(changed_nt_f).complement())

        mut_type_percentage = strand_bias*sbs.loc[mut_type, column]

        i1 = dna_bases.index(changed_nt_f)
        j1 = trinucleotides.index(Seq(trinucleotide_f))
        i2 = dna_bases.index(changed_nt_r)
        j2 = trinucleotides.index(Seq(trinucleotide_r))
        full_mutational_signature[i1,j1] = mut_type_percentage
        full_mutational_signature[i2,j2] = mut_type_percentage

        return full_mutational_signature


def compute_mut_prob_matrix(seq, full_mutational_signature):
    trint_vector = ["-"]
    trint_type_vector = [-1]
    for idx in range(1, len(seq)-1):
        trint_vector.append(seq[idx-1:idx+2])
        trint_type_vector.append(trinucleotides.index(seq[idx-1:idx+2]))
    trint_vector.append("-")
    trint_type_vector.append(-1)

    P_before_norm = np.zeros((4, len(seq)))

    for idx in range(1, len(seq)-1):
        P_before_norm[:, idx] = (
            full_mutational_signature[:, trint_type_vector[idx]]
        )

    S = np.sum(P_before_norm)
    P = P_before_norm/S

    return P


def mutate_seq_with_mutational_signature(
    seq: Union[str, Seq, MutableSeq],
    mutational_signature: str = "SBS1",
    cosmic_version: float = 3.3,
    genome_ref: str = "GRCh37",
    custom_signature_path: str = None,
    column: str = None,
    strand_bias: float = 0.5
) -> Tuple[Union[str, Seq, MutableSeq], int, str, int]:
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

    signatures = [
        SBS1, SBS2, SBS3, SBS4, SBS5, SBS6, SBS7a, SBS7b, SBS7c, SBS7d,  
    	SBS8, SBS9, SBS10a, SBS10b, SBS10c, SBS10d, SBS11, SBS12, SBS13,
        SBS14, SBS15, SBS16, SBS17a, SBS17b, SBS18, SBS19, SBS20, SBS21,
        SBS22, SBS23, SBS24, SBS25, SBS26, SBS27, SBS28, SBS29, SBS30,
        SBS31, SBS32, SBS33, SBS34, SBS35, SBS36, SBS37, SBS38, SBS39,
        SBS40, SBS41, SBS42, SBS43, SBS44, SBS45, SBS46, SBS47, SBS48,
        SBS49, SBS50, SBS51, SBS52, SBS53, SBS54, SBS55, SBS56, SBS57,
        SBS58, SBS59, SBS60, SBS84, SBS85, SBS86, SBS87, SBS88, SBS89,
        SBS90, SBS91, SBS92, SBS93, SBS94, SBS95
        ]

    if mutational_signature == "custom":
        if custom_signature_path is None:
            raise ValueError("custom_signature_path must be provided")
        else:
            full_mutational_signature = make_full_mutational_signature(
                input_file_path=custom_signature_path,
                column=column,
                strand_bias=strand_bias
                )
    elif mutational_signature in signatures:
        signature_file_path = (
            f"../cosmic_signatures/
            COSMIC_v{str(cosmic_version)}.1_SBS_{genome_ref}.txt"
        )
        full_mutational_signature = make_full_mutational_signature(
                input_file_path=signature_file_path,
                column=mutational_signature,
                strand_bias=strand_bias
                )
    else:
        raise ValueError("signature must be SBSXX or custom")

    P = compute_mut_prob_matrix(seq, full_mutational_signature)

    num_nonzero = np.count_nonzero(P)
    if not num_nonzero:
        return Seq(""), -1, "", 1
    
    idx_flat = np.random.choice(len(seq)*4, p=P.flatten("F"))[0]

    idx_base = idx_flat%4
    idx_target = math.floor(idx_flat/4)

    if str(seq[idx_target]) not in ["C", "T"]:
        description = f"{str(Seq(seq[idx_target+1]).complement())}[
            {str(Seq(seq[idx_target]).complement())}>
            {str(Seq(dna_bases[idx_base]).complement())}]
            {str(Seq(seq[idx_target-1]).complement())}"
    else:
        description = f"{str(seq[idx_target-1])}[
            {str(seq[idx_target])}>
            {dna_bases[idx_base]}]
            {str(seq[idx_target+1])}"

    seq[idx_target] = dna_bases[idx_base]

    if data_type == "str":
        return str(seq), idx_target, description, 0
    elif data_type == "Seq":
        return Seq(seq), idx_target, description, 0
    elif data_type == "MutableSeq":
        return seq, idx_target, description, 0
 

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
        return str(seq), idx, description
    elif data_type == "Seq":
        return Seq(seq), idx, description
    elif data_type == "MutableSeq":
        return seq, idx, description


def get_target_of_mut_type(mut_type):
    open_bracet_pos = re.search("\[", mut_type).start()
    change_symbol_pos = re.search("\>", mut_type).start()
    close_bracet_pos = re.search("\]", mut_type).start()

    before_bracet = mut_type[:open_bracet_pos]
    after_bracet = mut_type[close_bracet_pos+1:]
    seq_target = mut_type[open_bracet_pos+1:change_symbol_pos]

    seq_output = f"{before_bracet}{seq_target}{after_bracet}"

    return seq_output


def get_result_of_mut_type(mut_type):
    open_bracet_pos = re.search("\[", mut_type).start()
    change_symbol_pos = re.search("\>", mut_type).start()
    close_bracet_pos = re.search("\]", mut_type).start()

    before_bracet = mut_type[:open_bracet_pos]
    after_bracet = mut_type[close_bracet_pos+1:]
    seq_result = mut_type[change_symbol_pos+1:close_bracet_pos]

    seq_output = f"{before_bracet}{seq_result}{after_bracet}"

    return seq_output


def get_motif_indices(seq, mut_type):
    seq_target = get_target_of_mut_type(mut_type)
    idx_motifs_obj = re.finditer(pattern=seq_target, string=str(seq))
    idx_motifs = [idx.start() for idx in idx_motifs_obj]
    
    return idx_motifs


def get_idx_offset_of_mut_type(mut_type):
    idx_offset = re.search("\[", mut_type).start()

    return idx_offset


def get_idx_target_from_idx_motif(idx_motif, mut_type):
    idx_target = idx_motif+get_idx_offset_of_mut_type(mut_type)

    return idx_target


def get_idx_motif_from_idx_target(idx_target, mut_type):
    idx_motif = idx_target-get_idx_offset_of_mut_type(mut_type)

    return idx_motif
    
    
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

    idx_motifs = get_motif_indices(seq, mut_type)
    idx_motifs = [idx for idx in idx_motifs if idx >= start and idx <= end]

    if not idx_motifs:
        return Seq(""), -1, "", 1
    
    # select specific motif to mutate
    idx_motif = np.random.choice(idx_motifs, 1)[0] 
    idx_target = get_idx_target_from_idx_motif(idx_motif, mut_type)

    nt_before = seq[idx_target:idx_target+1]
    seq[idx_motif:idx_target+1] = get_result_of_mut_type(mut_type)
    nt_after = seq[idx_target:idx_target+1]
    
    description = f"c.{str(idx_target)}{nt_before}>{nt_after}"

    if data_type == "str":
        return str(seq), idx_target, description, 0
    elif data_type == "Seq":
        return Seq(seq), idx_target, description, 0
    elif data_type == "MutableSeq":
        return seq, idx_target, description, 0

    
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
        self._descriptions = ["original"]
        self._length = 0

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
        
    def add_mutated_sequence(self, method="random", mut_type=None, mutational_signature=None, cosmic_version=None, genome_ref=None, custom_signature_path=None, column=None, strand_bias=None):
        if method == "random":
            last_seq = self._data[-1]
            new_seq, idx_target, description = random_single_substitution(
                last_seq,
                start=1,
                end=self._seq_length-2
                )
            self.append(new_seq, description=description)
            self._length = self._length+1
        elif method == "mut_type":
            last_seq = self._data[-1]
            new_seq, idx_target, description, e = mutate_seq_with_mut_type(
                last_seq,
                mut_type=mut_type,
                start=1,
                end=self._seq_length-2
                )
            if not e:
                self.append(new_seq, description=description)
                self._length = self._length+1
            else:
                pass
        elif method == "signature":
            last_seq = self._data[-1]
            new_seq, idx_target, description, e = mutate_seq_with_mutational_signature(
                last_seq,
                mutational_signature=mutational_signature,
                cosmic_version=cosmic_version,
                genome_ref=genome_ref,
                custom_signature_path=custom_signature_path,
                column=column,
                strand_bias=strand_bias
            )
            if not e:
                self.append(new_seq, description=description)
                self._length = self._length+1
        else:
            raise ValueError(
                "mode must be either random, motif, or signature"
            )

        return self._data