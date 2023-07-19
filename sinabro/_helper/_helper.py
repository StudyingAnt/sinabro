import re
import math

import numpy as np

from Bio.Seq import Seq, MutableSeq
from Bio.Data import CodonTable

from typing import Union

class MutInfo:
    def __init__(
        self, 
        new_seq: Union[str, Seq, MutableSeq], 
        idx_target: int, 
        hgvs_mrna: str,
        mut_type: str,
        e: int
        ) -> None:
        self.new_seq = new_seq
        self.idx_target = idx_target
        self.hgvs_mrna = hgvs_mrna
        self.mut_type = mut_type
        self.e = e

def _reverse_complement_mut_type(
        mut_type: str
    ) -> str:

    open_bracet_pos = re.search("\[", mut_type).start()
    change_symbol_pos = re.search("\>", mut_type).start()
    close_bracet_pos = re.search("\]", mut_type).start()

    before_bracet = mut_type[:open_bracet_pos]
    after_bracet = mut_type[close_bracet_pos+1:]

    seq_target = mut_type[open_bracet_pos+1:change_symbol_pos]
    seq_result = mut_type[change_symbol_pos+1:close_bracet_pos]

    # reverse complement
    new_before_bracet = str(Seq(after_bracet).reverse_complement())
    new_after_bracet = str(Seq(before_bracet).reverse_complement())
    new_seq_target = str(Seq(seq_target).reverse_complement())
    new_seq_result = str(Seq(seq_result).reverse_complement())

    revcomp_mut_type = f"{new_before_bracet}[{new_seq_target}>{new_seq_result}]{new_after_bracet}"

    return revcomp_mut_type


def _get_target_of_mut_type(
    mut_type: str,
    ) -> str:
    open_bracet_pos = re.search("\[", mut_type).start()
    change_symbol_pos = re.search("\>", mut_type).start()
    close_bracet_pos = re.search("\]", mut_type).start()

    before_bracet = mut_type[:open_bracet_pos]
    after_bracet = mut_type[close_bracet_pos+1:]
    seq_target = mut_type[open_bracet_pos+1:change_symbol_pos]

    seq_output = f"{before_bracet}{seq_target}{after_bracet}"

    return seq_output


def _get_result_of_mut_type(
    mut_type: str,
    ) -> str:
    open_bracet_pos = re.search("\[", mut_type).start()
    change_symbol_pos = re.search("\>", mut_type).start()
    close_bracet_pos = re.search("\]", mut_type).start()

    before_bracet = mut_type[:open_bracet_pos]
    after_bracet = mut_type[close_bracet_pos+1:]
    seq_result = mut_type[change_symbol_pos+1:close_bracet_pos]

    seq_output = f"{before_bracet}{seq_result}{after_bracet}"

    return seq_output


def _get_motif_indices(
    seq: Union[str, Seq, MutableSeq], 
    mut_type) -> list:
    if isinstance(seq, str):
        pass
    elif isinstance(seq, Seq):
        seq = str(seq)
    elif isinstance(seq, MutableSeq):
        seq = str(seq)
    else:
        raise TypeError(
            "seq should be a string, Seq, or MutableSeq object"
        )
    seq_target = _get_target_of_mut_type(mut_type)
    idx_motifs_obj = re.finditer(pattern=seq_target, string=seq)
    idx_motifs = [idx.start() for idx in idx_motifs_obj]
    
    return idx_motifs


def _get_idx_offset_of_mut_type(
    mut_type: str
    ) -> int:
    idx_offset = re.search("\[", mut_type).start()

    return idx_offset


def _get_idx_target_from_idx_motif(
    idx_motif: int,
    mut_type: str
    ) -> int:
    idx_target = idx_motif+_get_idx_offset_of_mut_type(mut_type)

    return idx_target


def _get_idx_motif_from_idx_target(
    idx_target: int,
    mut_type: str
    ) -> int:
    idx_motif = idx_target-_get_idx_offset_of_mut_type(mut_type)

    return idx_motif


def _is_codon_synonymous(codon1, codon2):
    codon_table = CodonTable.standard_dna_table.forward_table
    codon_table["TAA"] = "*"
    codon_table["TAG"] = "*"
    codon_table["TGA"] = "*"
    
    if codon_table[codon1] == codon_table[codon2]:
        return True
    else:
        return False


def _get_idx_codon(idx_target):
    idx_codon = int(math.floor((idx_target-1)/3)*3+1)

    return idx_codon


def _get_pos_aa(idx_target):
    pos_aa = int(math.floor((idx_target-1)/3)+1)

    return pos_aa


def _get_codon(seq, idx_target):
    if isinstance(seq, str):
        pass
    elif isinstance(seq, Seq):
        seq = str(seq)
    elif isinstance(seq, MutableSeq):
        seq = str(seq)
    else:
        raise TypeError(
            "seq should be a string, Seq, or MutableSeq object"
        )

    idx_codon = _get_idx_codon(idx_target)

    codon = seq[idx_codon:idx_codon+3]

    return codon


def _get_idx_from_hgvs_mrna(hgvs_mrna, offset=0):
    changed_symbol_pos = re.search(">", hgvs_mrna).start()

    return int(hgvs_mrna[2:changed_symbol_pos-1])-1+offset


def _get_amino_acid_from_codon(codon, three=False):
    codon_table = CodonTable.standard_dna_table.forward_table
    codon_table["TAA"] = "*"
    codon_table["TAG"] = "*"
    codon_table["TGA"] = "*"

    one_to_three = {
        "A": "Ala", "R": "Arg", "N": "Asn", "D": "Asp", "C": "Cys",
        "E": "Glu", "Q": "Gln", "G": "Gly", "H": "His", "I": "Ile",
        "L": "Leu", "K": "Lys", "M": "Met", "F": "Phe", "P": "Pro",
        "S": "Ser", "T": "Thr", "W": "Trp", "Y": "Tyr", "V": "Val",
        "*": "*"
    }

    if three:
        return one_to_three[codon_table[codon]] 
    else:
        return codon_table[codon]


def generate_random_sequence(length, gc=0.5):
    dna = ["A", "C", "G", "T"]
    seq = Seq("".join(
        np.random.choice(
            dna, 
            length, 
            p=0.5*np.array([1-gc, gc, gc, 1-gc])
            )
        )
    )
    
    return seq