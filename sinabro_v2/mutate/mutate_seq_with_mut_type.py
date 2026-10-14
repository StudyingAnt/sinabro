from sinabro._helper._helper import _unpack_complex_mut_type, _reverse_complement_mut_type, _get_target_of_mut_type, _get_result_of_mut_type, _get_motif_indices, _get_idx_offset_of_mut_type, _get_idx_target_from_idx_motif, _get_idx_motif_from_idx_target
from sinabro._mutate._mutate_seq_with_mut_type import mutate_seq_with_mut_type, mutate_seq_with_mut_types

def _mutate_seq_with_mut_type_basic(seq, mut_type, **kwargs):
    return mutate_seq_with_mut_type(seq, mut_type, strand='single', **kwargs)

def _mutate_seq_with_mut_type_reverse(seq, mut_type, **kwargs):
    return _mutate_seq_with_mut_type_basic(seq, _reverse_complement_mut_type(mut_type), **kwargs)

def _mutate_seq_with_mut_type_both(seq, mut_type, strand_bias=0.5, **kwargs):
    return mutate_seq_with_mut_type(seq, mut_type, strand_bias=strand_bias, **kwargs)

def _mutate_seq_with_mut_type_complex(seq, mut_type, mut_type_bias=None, strand='both', strand_bias=0.5, **kwargs):
    if mut_type_bias is not None:
        motifs = _unpack_complex_mut_type(mut_type)
        return mutate_seq_with_mut_types(seq, motifs, mut_type_bias, strand=strand, strand_bias=strand_bias, **kwargs)
    return mutate_seq_with_mut_type(seq, mut_type, strand=strand, strand_bias=strand_bias, **kwargs)
