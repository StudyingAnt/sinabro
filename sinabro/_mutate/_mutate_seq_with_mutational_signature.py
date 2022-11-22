import os
import math

import numpy as np
import pandas as pd

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
    start: int = None,
    end: int = None,
    mutational_signature: str = "SBS1",
    cosmic_version: float = 3.3,
    genome_ref: str = "GRCh37",
    custom_signature_path: str = None,
    column: str = None,
    strand_bias: float = 0.5
) -> tuple[Union[str, Seq, MutableSeq], int, str, int]:
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
        "SBS1", "SBS2", "SBS3", "SBS4", "SBS5", "SBS6", "SBS7a", "SBS7b", 
        "SBS7c", "SBS7d", "SBS8", "SBS9", "SBS10a", "SBS10b", "SBS10c", 
        "SBS10d", "SBS11", "SBS12", "SBS13", "SBS14", "SBS15", "SBS16", 
        "SBS17a", "SBS17b", "SBS18", "SBS19", "SBS20", "SBS21", "SBS22", 
        "SBS23", "SBS24", "SBS25", "SBS26", "SBS27", "SBS28", "SBS29", 
        "SBS30", "SBS31", "SBS32", "SBS33", "SBS34", "SBS35", "SBS36", 
        "SBS37", "SBS38", "SBS39", "SBS40", "SBS41", "SBS42", "SBS43", 
        "SBS44", "SBS45", "SBS46", "SBS47", "SBS48", "SBS49", "SBS50", 
        "SBS51", "SBS52", "SBS53", "SBS54", "SBS55", "SBS56", "SBS57",
        "SBS58", "SBS59", "SBS60", "SBS84", "SBS85", "SBS86", "SBS87", 
        "SBS88", "SBS89", "SBS90", "SBS91", "SBS92", "SBS93", "SBS94", 
        "SBS95"
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
        dir_name = os.path.dirname(__file__)
        tokens = [
            "../../cosmic_signatures/COSMIC_v",
            str(cosmic_version), ".1_SBS_",
            genome_ref, ".txt"
        ]
        signature_file = "".join(tokens)
        signature_file_path = os.path.join(dir_name, signature_file)
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
        return MutInfo(MutableSeq(""), -1, "", "", 1)
   
    idx_flat = np.random.choice(len(seq)*4, p=P.flatten("F"))

    idx_base = idx_flat%4
    idx_target = math.floor(idx_flat/4)

    if str(seq[idx_target]) not in ["C", "T"]:
        hgvs_mrna_tokens = [
            "c.", str(idx_target+1),
            str(seq[idx_target]), ">",
            dna_bases[idx_base]
        ]
        mut_type_tokens = [
            str(Seq(seq[idx_target+1]).complement()), "[", 
            str(Seq(seq[idx_target]).complement()), ">", 
            str(Seq(dna_bases[idx_base]).complement()), "]",
            str(Seq(seq[idx_target-1]).complement())
        ]
        hgvs_mrna = "".join(hgvs_mrna_tokens)
        mut_type = "".join(mut_type_tokens)
    else:
        hgvs_mrna_tokens = [
            "c.", str(idx_target+1),
            str(seq[idx_target]), ">",
            dna_bases[idx_base]
        ]
        mut_type_tokens = [
            str(seq[idx_target-1]), "[",
            str(seq[idx_target]), ">",
            dna_bases[idx_base], "]",
            str(seq[idx_target+1])    
        ]
        hgvs_mrna = "".join(hgvs_mrna_tokens)
        mut_type = "".join(mut_type_tokens)
        
    seq[idx_target] = dna_bases[idx_base]

    if data_type == "str":
        return MutInfo(str(seq), idx_target, hgvs_mrna, mut_type, 0)
    elif data_type == "Seq":
        return MutInfo(Seq(seq), idx_target, hgvs_mrna, mut_type, 0)
    elif data_type == "MutableSeq":
        return MutInfo(seq, idx_target, hgvs_mrna, mut_type, 0)