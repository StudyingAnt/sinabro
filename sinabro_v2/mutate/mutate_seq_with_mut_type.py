import re
import numpy as np

from typing import Union
from Bio.Seq import Seq, MutableSeq

from ..types.types import MutInfo, MutationRecord
from ._helper import preserve_seq_type

def _unpack_complex_mut_type(
        mut_type: str
    ) -> str:


    non_single_iupac_nts = [
        "R", "Y", "S", "W", "K", "M", 
        "B", "D", "H", "V",
        "N"
        ]
    
    iupac_table = {
        "R": ["A", "G"],
        "Y": ["C", "T"],
        "S": ["G", "C"],
        "W": ["A", "T"],
        "K": ["G", "T"],
        "M": ["A", "C"],
        "B": ["C", "G", "T"],
        "D": ["A", "G", "T"],
        "H": ["A", "C", "T"],
        "V": ["A", "C", "G"],
        "N": ["A", "C", "G", "T"]
    }

    mut_types_pre = []
    for c in mut_type:
        if c in non_single_iupac_nts:
            mut_types_pre.append(iupac_table[c])
        else:
            mut_types_pre.append([c])

    mut_types_split = list(itertools.product(*mut_types_pre))

    mut_types = []
    for new_mut_type in mut_types_split:
        mut_types.append("".join(new_mut_type))       

    return mut_types

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

@preserve_seq_type
def _mutate_seq_with_mut_type_basic(seq, mut_type):
    start = 0
    end = len(seq) - 1

    idx_motifs = [idx for idx in _get_motif_indices(seq, mut_type) if start <= idx <= end]
    if not idx_motifs:
        return MutInfo(seq, -1, "", "", 1)

    # Create a random generator instance.
    rng = np.random.default_rng()
    
    idx_motif = rng.choice(idx_motifs)

    idx_target = _get_idx_target_from_idx_motif(idx_motif, mut_type)
    nt_before = seq[idx_target:idx_target+1]
    mutated_seq_str = _get_result_of_mut_type(mut_type)
    len_motif = len(mutated_seq_str)
    seq[idx_motif: idx_motif+len_motif] = mutated_seq_str
    nt_after = seq[idx_target:idx_target+1]
    hgvs_mrna = f"c.{idx_target-start}{nt_before}>{nt_after}"
    
    # Always work with MutableSeq type inside the function.
    return MutInfo(seq, idx_target, hgvs_mrna, mut_type, 0)

@preserve_seq_type
def _mutate_seq_with_mut_type_reverse(seq, mut_type):
    start = 0
    end = len(seq) - 1

    rev_mut_type = _reverse_complement_mut_type(mut_type)
    
    return _mutate_seq_with_mut_type_basic(seq, rev_mut_type)

@preserve_seq_type
def _mutate_seq_with_mut_type_both(seq, mut_type, strand_bias=0.5):
    start = 0
    end = len(seq) - 1

    # Generate a random number between 0 and 1.
    rand_val = np.random.uniform(0, 1)
    
    # Compute the reverse mutation type.
    rev_mut_type = _reverse_complement_mut_type(mut_type)
    
    # Decide the primary and secondary mutation types based on the bias.
    primary_mut = mut_type if rand_val < strand_bias else rev_mut_type
    secondary_mut = rev_mut_type if rand_val < strand_bias else mut_type
    
    # Attempt mutation using the primary mutation type.
    result = _mutate_seq_with_mut_type_basic(seq, primary_mut)
    if result.e == 0:
        return result
    # If the primary attempt fails, try the secondary mutation type.
    return _mutate_seq_with_mut_type_basic(seq, secondary_mut)

@preserve_seq_type
def _mutate_seq_with_mut_type_complex(seq, mut_type, mut_type_bias=None, strand="both", strand_bias=0.5):
    start = 0
    end = len(seq) - 1




@preserve_seq_type
def mutate_seq_with_mut_type(seq, **kwargs):
    # Set default parameter values
    defaults = {
        'mut_type': None,
        'strand': 0.5
    }
    # Update defaults with values provided in kwargs (if any)
    defaults.update(kwargs)
    
    # Extract parameters from the defaults dictionary
    mut_type = defaults['mut_type']
    strand = defaults['strand']

    # Check mut_type provided.
    if mut_type is None:
        raise ValueError("A value for 'mut_type' must be provided.")

    # Validate strand parameter.
    if strand not in ("both", "single"):
        raise ValueError("strand should be either 'single' or 'both', default is 'both'")

    start = 0
    end = len(seq) - 1

    # Determine if the mutation type is simple (only single nucleotide IUPAC symbols allowed).
    mut_type_nts = [nt for nt in mut_type if nt not in ["[", ">", "]"]]
    non_single_iupac_nts = ["R", "Y", "S", "W", "K", "M", "B", "D", "H", "V", "N"]
    simple_mut_type = not (set(mut_type_nts) & set(non_single_iupac_nts))
    
    # Create a random generator instance.
    rng = np.random.default_rng()

    

    

    # Always work with MutableSeq type inside the function.
    return MutInfo(seq, idx_target, hgvs_mrna, mut_type, 0)
    
