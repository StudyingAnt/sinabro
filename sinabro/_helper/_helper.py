"""Historical helper names backed by the shared coordinate model."""
import re
import numpy as np
from Bio.Seq import Seq
from ..core import MutInfo, expand_motif as _unpack_complex_mut_type, reverse_motif as _reverse_complement_mut_type, motif_target as _get_target_of_mut_type, motif_result as _get_result_of_mut_type, motif_indices as _get_motif_indices, synonymous as _is_codon_synonymous, codon_start as _get_idx_codon, get_codon as _get_codon, amino_acid as _get_amino_acid_from_codon, rng_from

def _get_idx_offset_of_mut_type(mut_type):
    from ..core import parse_motif
    return len(parse_motif(mut_type)[0])

def _get_idx_target_from_idx_motif(idx_motif, mut_type):
    return idx_motif + _get_idx_offset_of_mut_type(mut_type)

def _get_idx_motif_from_idx_target(idx_target, mut_type):
    return idx_target - _get_idx_offset_of_mut_type(mut_type)

def _get_idx_from_hgvs_mrna(hgvs_mrna, offset=1):
    match = re.fullmatch('c\\.(\\d+)[ACGT]>[ACGT]', hgvs_mrna)
    if not match or int(match[1]) < 1:
        raise ValueError('invalid HGVS substitution')
    return int(match[1]) - 1 + offset

def _get_pos_aa(idx_target, start=1):
    return (_get_idx_codon(idx_target, start) - start) // 3 + 1

def generate_random_sequence(length, gc=0.5, *, rng=None, seed=None):
    from ..core import integer
    integer(length, 'length', 1)
    if not np.isfinite(gc) or not 0 <= gc <= 1:
        raise ValueError('gc must be between zero and one')
    return Seq(''.join(rng_from(rng, seed).choice(list('ACGT'), length, p=np.array([1 - gc, gc, gc, 1 - gc]) / 2)))
