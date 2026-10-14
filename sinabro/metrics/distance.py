"""Distances are amino-acid distances after translating the specified CDS."""
from ..evaluation.policies import protein, alignment_distance
from ..core import coding_bounds
from Bio.Seq import Seq

def dna_to_aa(seq, table=1, to_stop=False):
    coding_bounds(seq, 0, len(seq))
    return str(Seq(str(seq)).translate(table=table, to_stop=to_stop))

def hamming(a, b):
    (a, b) = (protein(a), protein(b))
    if len(a) != len(b):
        raise ValueError('proteins must have equal length')
    return sum((x != y for (x, y) in zip(a, b)))

def compute_distance(records, **kwargs):
    if not records:
        raise ValueError('at least one record is required')
    include = kwargs.get('threshold_hit', False)
    if not include and len(records) < 2:
        raise ValueError('excluding terminal record requires at least two records')
    start = kwargs.get('cds_start', 1)
    end = kwargs.get('cds_end', len(records[0].sequence) - 1)
    (ref, _, _) = coding_bounds(records[0].sequence, start, end)
    (cur, _, _) = coding_bounds(records[-1 if include else -2].sequence, start, end)
    if len(ref) != len(cur):
        raise ValueError('trajectory sequences must have equal length')
    method = kwargs.get('method', 'hamming')
    if method == 'hamming':
        return hamming(ref[start:end], cur[start:end])
    if method == 'blosum':
        return alignment_distance(ref[start:end], cur[start:end])
    raise ValueError('method must be hamming or blosum')
