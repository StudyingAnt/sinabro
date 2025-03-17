from typing import Union
from Bio.Seq import Seq, MutableSeq

from ..types.types import MutInfo, MutationRecord
from ..utils import get_codon, get_amino_acid_from_codon

def eval_nonsym(records, mutinfo, **kwargs):
    old_seq = records[-2].sequence
    new_seq = records[-1].sequence

    old_aa = get_amino_acid_from_codon(get_codon(old_seq, mutinfo.idx_target))
    new_aa = get_amino_acid_from_codon(get_codon(new_seq, mutinfo.idx_target))

    if old_aa != new_aa:
        return True, None
    else:
        return False, None