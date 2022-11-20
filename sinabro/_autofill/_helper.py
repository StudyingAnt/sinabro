import re
import math

from Bio.Data import CodonTable
from Bio.Seq import Seq, MutableSeq


def _is_codon_synonymous(codon1, codon2):
    codon_table = CodonTable.standard_dna_table.forward_table
    codon_table["TAA"] = "-"
    codon_table["TAG"] = "-"
    codon_table["TGA"] = "-"
    
    if codon_table[codon1] == codon_table[codon2]:
        return True
    else:
        return False


def _get_idx_codon(idx_target):
    idx_codon = int(math.floor((idx_target-1)/3)*3+1)

    return idx_codon


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


def _get_idx_from_hgvs(hgvs):
    changed_symbol_pos = re.search(">", hgvs).start()

    return int(hgvs[2:changed_symbol_pos-1])-1
    