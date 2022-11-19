import re

from Bio.Seq import Seq, MutableSeq
from typing import Union


class MutInfo:
    def __init__(
        self, 
        new_seq: Union[str, Seq, MutableSeq], 
        idx_target: int, 
        hgvs: str,
        mut_type: str,
        e: int
        ) -> None:
        self.new_seq = new_seq
        self.idx_target = idx_target
        self.hgvs = hgvs
        self.mut_type = mut_type
        self.e = e


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